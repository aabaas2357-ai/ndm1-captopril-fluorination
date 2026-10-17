"""
charge_swap_test.py - does the carboxylate-Zn collapse during minimization depend on the ligand charges?
Starts from the crystal-pose structure md_5ZJ2/solvated.pdb and minimizes twice with the SAME system:
  A) original charges (Gasteiger)   B) AM1-BCC charges from amber_lig/cap_gaff2_xtal.mol2
Only ligand charges (and the matching 1-4 exceptions) change; types, LJ and bonded terms stay as they are.
Run in C:\\NDM1_Project\\md with the venv python (openmm + numpy only).
"""
import numpy as np
import openmm as mm
from openmm import app, unit

OUT = "md_5ZJ2"
MOL2 = "amber_lig/cap_gaff2_xtal.mol2"
SCALE14 = 1.0 / 1.2


def read_mol2(path):
    atoms, sect = [], None
    for line in open(path, encoding="utf-8-sig", errors="replace"):
        if line.startswith("@<TRIPOS>"):
            sect = line.strip()
            continue
        if sect == "@<TRIPOS>ATOM" and line.strip():
            p = line.split()
            atoms.append(dict(name=p[1], xyz=np.array([float(p[2]), float(p[3]), float(p[4])]),
                              elem=p[1][0] if p[1][0] != "H" else "H", q=float(p[8])))
    return atoms


pdb = app.PDBFile(OUT + "/solvated.pdb")
top = pdb.topology
pos = pdb.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
atoms = list(top.atoms())
lig = [a for a in atoms if a.residue.name == "CAP"]
bcc = read_mol2(MOL2)
print("ligand atoms in system: %d   in mol2: %d" % (len(lig), len(bcc)))

# --- map BCC charges onto the MD ligand atoms
bcc_heavy = {b["name"]: b for b in bcc if b["elem"] != "H"}
newq = {}
worst = 0.0
for a in lig:
    if a.element.symbol != "H":
        if a.name not in bcc_heavy:
            raise SystemExit("heavy atom name %s not in mol2 (names: %s)" % (a.name, sorted(bcc_heavy)))
        b = bcc_heavy[a.name]
        worst = max(worst, float(np.linalg.norm(pos[a.index] - b["xyz"])))
        newq[a.index] = b["q"]
print("max heavy-atom position difference crystal vs mol2: %.3f A (should be < 0.05)" % worst)
if worst > 0.05:
    raise SystemExit("atom names do not map onto the same atoms - stop and send me this output")


def parent(xyz, heavy):
    d = [(float(np.linalg.norm(xyz - h["xyz"])), h["name"]) for h in heavy]
    return min(d)[1]


hb = [b for b in bcc if b["elem"] == "H"]
heavy_b = list(bcc_heavy.values())
grp = {}
for b in hb:
    grp.setdefault(parent(b["xyz"], heavy_b), []).append(b["q"])
for a in lig:
    if a.element.symbol == "H":
        pn = parent(pos[a.index], [dict(name=n, xyz=pos[x.index]) for n, x in
                                    [(y.name, y) for y in lig if y.element.symbol != "H"]])
        if pn not in grp:
            raise SystemExit("hydrogen %s has no counterpart on %s" % (a.name, pn))
        newq[a.index] = float(np.mean(grp[pn]))
tot = sum(newq.values())
shift = (-2.0 - tot) / len(newq)
for k in newq:
    newq[k] += shift
print("AM1-BCC total %.4f -> shifted by %.5f per atom to exactly -2.0000" % (tot, shift))

system0 = mm.XmlSerializer.deserialize(open(OUT + "/system.xml").read())
nb0 = [f for f in system0.getForces() if isinstance(f, mm.NonbondedForce)][0]
oldq = {i: nb0.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in newq}
print("original ligand total charge: %.4f" % sum(oldq.values()))

systemB = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(system0))
nb = [f for f in systemB.getForces() if isinstance(f, mm.NonbondedForce)][0]
for i, q in newq.items():
    _, sg, ep = nb.getParticleParameters(i)
    nb.setParticleParameters(i, q, sg, ep)
qall = {i: nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(nb.getNumParticles())}
nfix = 0
for k in range(nb.getNumExceptions()):
    i, j, cp, sg, ep = nb.getExceptionParameters(k)
    if i in newq or j in newq:
        if cp.value_in_unit(unit.elementary_charge ** 2) != 0.0 or ep.value_in_unit(unit.kilojoule_per_mole) != 0.0:
            nb.setExceptionParameters(k, i, j, SCALE14 * qall[i] * qall[j], sg, ep)
            nfix += 1
print("1-4 exceptions rescaled for the new charges:", nfix)

zn = [a.index for a in atoms if a.element is not None and a.element.symbol == "Zn"]
byname = {a.name: a.index for a in lig}
lys = [a.index for a in atoms if a.residue.name == "LYS" and a.residue.id == "170" and a.name == "NZ"][0]
asn = [a.index for a in atoms if a.residue.name == "ASN" and a.residue.id == "179" and a.name == "N"][0]


def row(P):
    d = lambda i, j: float(np.linalg.norm(P[i] - P[j]))
    oz = [min(d(byname[o], z) for z in zn) for o in ("O1", "O2", "O3")]
    sz = [d(byname["S1"], z) for z in zn]
    lk = [d(byname[o], lys) for o in ("O1", "O2")]
    an = [d(byname[o], asn) for o in ("O1", "O2", "O3")]
    return "%6.2f %6.2f %6.2f | %5.2f %5.2f | %9.2f %9.2f | %8.2f %8.2f %8.2f" % tuple(oz + sz + lk + an)


def minimize(system):
    integ = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    try:
        sim = app.Simulation(top, system, integ, mm.Platform.getPlatformByName("OpenCL"),
                             {"OpenCLPlatformIndex": "0", "DeviceIndex": "0"})
    except Exception:
        integ = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        sim = app.Simulation(top, system, integ, mm.Platform.getPlatformByName("CPU"))
    sim.context.setPositions(pdb.positions)
    sim.minimizeEnergy()
    return sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)


print("\ncolumns (A): O1-Zn O2-Zn O3-Zn | S-Zn(a) S-Zn(b) | O1-Lys NZ O2-Lys NZ | O1-Asn N O2-Asn N O3-Asn N")
print("crystal (solvated.pdb)        ", row(pos))
try:
    pm = app.PDBFile(OUT + "/minimized.pdb").getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    print("minimized.pdb (old run)       ", row(pm))
except Exception as e:
    print("minimized.pdb not read:", e)
print("A) original charges, re-min   ", row(minimize(system0)))
print("B) AM1-BCC charges, min       ", row(minimize(systemB)))