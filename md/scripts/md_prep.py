"""
md_prep.py (v2, no OpenFF) - build the NDM-1 / D-captopril (5ZJ2) system for OpenMM,
minimize it on the GPU and time a short run.
Usage:  python md_prep.py 5ZJ2_clean.pdb mco_crystal.pdb

Ligand parameters: BASELINE. Atom types are borrowed from amber14 (PRO/CPRO/CYM/ALA),
partial charges are Gasteiger rescaled to the net charge (-2). Not validated.
"""
import sys, os, json, time
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import openmm as mm
from openmm import app, unit

PROT = sys.argv[1] if len(sys.argv) > 1 else "5ZJ2_clean.pdb"
LIG = sys.argv[2] if len(sys.argv) > 2 else "mco_crystal.pdb"
OUT = "md_5ZJ2"
K_ZN = 10000.0  # kJ/mol/nm^2, zinc coordination restraints
LIG_NAME = "CAP"
os.makedirs(OUT, exist_ok=True)

ZN_SITES = {
    "301": [("124", "OD2"), ("208", "SG"), ("250", "NE2")],
    "302": [("120", "NE2"), ("122", "ND1"), ("189", "NE2")],
}
HIS_STATES = {"120": "HID", "122": "HIE", "189": "HID", "250": "HID"}
CYS_STATES = {"208": "CYX"}

# captopril heavy atoms -> (amber14 template, atom) used only to borrow atom TYPES
SRC = {1: ("ALA", "CB"), 2: ("ALA", "CB"), 3: ("CYM", "CB"), 4: ("CYM", "SG"),
       5: ("PRO", "C"), 6: ("PRO", "O"), 7: ("PRO", "N"), 8: ("PRO", "CD"),
       9: ("PRO", "CG"), 10: ("PRO", "CB"), 11: ("PRO", "CA"),
       12: ("CPRO", "C"), 13: ("CPRO", "O"), 14: ("CPRO", "OXT")}
HSRC = {1: ("ALA", "HB1"), 2: ("ALA", "HB1"), 3: ("CYM", "HB2"), 8: ("PRO", "HD2"),
        9: ("PRO", "HG2"), 10: ("PRO", "HB2"), 11: ("PRO", "HA")}


def xyz(pos):
    return np.array(pos.value_in_unit(unit.nanometer))


def find_atom(top, resid, name, resname=None):
    for res in top.residues():
        if res.id == resid and (resname is None or res.name == resname):
            for a in res.atoms():
                if a.name == name:
                    return a.index
    raise SystemExit("atom not found: residue %s atom %s" % (resid, name))


def ptype(ff, res, atom):
    if res not in ff._templates:
        raise SystemExit("template %s not found in amber14" % res)
    for a in ff._templates[res].atoms:
        if a.name == atom:
            return a.type
    raise SystemExit("atom %s not in template %s (has: %s)" %
                     (atom, res, [a.name for a in ff._templates[res].atoms]))


def build_ligand(path):
    tmpl = Chem.MolFromSmiles("CC(CS)C(=O)N1CCCC1C(=O)O")
    raw = Chem.MolFromPDBFile(path, removeHs=True, proximityBonding=True)
    if raw is None or raw.GetNumHeavyAtoms() != 14:
        raise SystemExit("ligand PDB must contain the 14 heavy atoms of captopril")
    mol = AllChem.AssignBondOrdersFromTemplate(tmpl, raw)
    Chem.AssignStereochemistryFrom3D(mol)
    for a in mol.GetAtoms():
        if a.GetSymbol() == "S":
            a.SetFormalCharge(-1)  # thiolate
    patt = Chem.MolFromSmarts("[CX3](=[OX1])[OX2;D1]")
    for m in mol.GetSubstructMatches(patt):
        mol.GetAtomWithIdx(m[2]).SetFormalCharge(-1)  # carboxylate
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol, addCoords=True)
    print("ligand atoms:", mol.GetNumAtoms(), " net charge:", Chem.GetFormalCharge(mol))
    cc = Chem.FindMolChiralCenters(mol, useLegacyImplementation=False)
    print("chiral centers (index, label, in ring):",
          [(i, l, mol.GetAtomWithIdx(i).IsInRing()) for i, l in cc])
    ok = all((l == "R") if mol.GetAtomWithIdx(i).IsInRing() else (l == "S") for i, l in cc)
    print("stereo check vs D-captopril (ring R, chain S):", "PASS" if ok else "FAIL")
    return mol


def make_ligand_ffxml(mol, ff, xml_path):
    patt = Chem.MolFromSmarts(
        "[CX4:1][CX4:2]([CX4:3][SX1:4])[CX3:5](=[OX1:6])[NX3:7]1[CX4:8][CX4:9][CX4:10][CX4:11]1"
                "[CX3:12](~[OX1:13])~[OX1:14]")
    match = mol.GetSubstructMatch(patt)
    if not match:
        raise SystemExit("captopril pattern not found in the ligand")
    label_of = {idx: k + 1 for k, idx in enumerate(match)}
    types = {}
    for k, idx in label_of.items():
        pass
    for idx, k in label_of.items():
        types[idx] = ptype(ff, *SRC[k])
    for a in mol.GetAtoms():
        if a.GetSymbol() == "H":
            heavy = a.GetNeighbors()[0].GetIdx()
            k = label_of[heavy]
            res, nm = HSRC.get(k, ("PRO", "HB2"))
            types[a.GetIdx()] = ptype(ff, res, nm)
    AllChem.ComputeGasteigerCharges(mol)
    q = [a.GetDoubleProp("_GasteigerCharge") for a in mol.GetAtoms()]
    shift = (Chem.GetFormalCharge(mol) - sum(q)) / len(q)
    q = [x + shift for x in q]
    counts, names = {}, []
    for a in mol.GetAtoms():
        s = a.GetSymbol()
        counts[s] = counts.get(s, 0) + 1
        names.append("%s%d" % (s, counts[s]))
    lines = ["<ForceField>", "<Residues>", '<Residue name="%s">' % LIG_NAME]
    for i, a in enumerate(mol.GetAtoms()):
        lines.append('<Atom name="%s" type="%s" charge="%.5f"/>' % (names[i], types[i], q[i]))
    for b in mol.GetBonds():
        lines.append('<Bond atomName1="%s" atomName2="%s"/>' %
                     (names[b.GetBeginAtomIdx()], names[b.GetEndAtomIdx()]))
    lines += ["</Residue>", "</Residues>", "</ForceField>"]
    open(xml_path, "w").write("\n".join(lines))
    print("ligand atom types:", sorted(set(types.values())))
    print("ligand Gasteiger charge range: %.2f .. %.2f" % (min(q), max(q)))
    return names


# 1) protein + zinc
pdb = app.PDBFile(PROT)
print("zinc residues:", [(r.name, r.id) for r in pdb.topology.residues() if r.name.upper() == "ZN"])
ff_prot = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
variants = []
for r in pdb.topology.residues():
    if r.name in ("HIS", "HID", "HIE", "HIP"):
        variants.append(HIS_STATES.get(r.id))
    elif r.name in ("CYS", "CYM", "CYX"):
        variants.append(CYS_STATES.get(r.id))
    else:
        variants.append(None)
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(ff_prot, pH=7.0, variants=variants)

# 2) ligand
mol = build_ligand(LIG)
xml_path = os.path.join(OUT, "cap.xml")
names = make_ligand_ffxml(mol, ff_prot, xml_path)
lig_top = app.Topology()
chain = lig_top.addChain()
res = lig_top.addResidue(LIG_NAME, chain, "900")
atoms = [lig_top.addAtom(names[i], app.Element.getBySymbol(a.GetSymbol()), res)
         for i, a in enumerate(mol.GetAtoms())]
for b in mol.GetBonds():
    lig_top.addBond(atoms[b.GetBeginAtomIdx()], atoms[b.GetEndAtomIdx()])
coords = mol.GetConformer().GetPositions() / 10.0
lig_pos = unit.Quantity([mm.Vec3(float(p[0]), float(p[1]), float(p[2])) for p in coords], unit.nanometer)
modeller.add(lig_top, lig_pos)

# 3) crystal distances for the zinc restraints
top = modeller.topology
pos0 = xyz(modeller.positions)
s_idx = [a.index for a in top.atoms() if a.element is not None and a.element.symbol == "S"
         and a.residue.name == LIG_NAME]
if len(s_idx) != 1:
    raise SystemExit("could not find the single ligand sulfur")
pairs = []
for zn_id, partners in ZN_SITES.items():
    zn = find_atom(top, zn_id, "ZN")
    for resid, name in partners:
        pairs.append((zn, find_atom(top, resid, name), "Zn%s-%s%s" % (zn_id, resid, name)))
    pairs.append((zn, s_idx[0], "Zn%s-S(ligand)" % zn_id))
restr = []
print("\nCrystal coordination distances (Angstrom):")
for zn, j, label in pairs:
    d = float(np.linalg.norm(pos0[zn] - pos0[j]))
    print("  %-22s %.2f%s" % (label, d * 10, "   <-- CHECK (too long)" if d > 3.0 else ""))
    restr.append((zn, j, d, label))

# 4) solvate + build the system
ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml", xml_path)
modeller.addSolvent(ff, model="tip3p", padding=1.0 * unit.nanometer,
                    ionicStrength=0.15 * unit.molar, neutralize=True)
top = modeller.topology
system = ff.createSystem(top, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometer,
                         constraints=app.HBonds, rigidWater=True)
bf = mm.HarmonicBondForce()
for zn, j, d, label in restr:
    bf.addBond(zn, j, d, K_ZN)
system.addForce(bf)

nb = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
lig_atoms = [a.index for a in top.atoms() if a.residue.name == LIG_NAME]
q_lig = sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in lig_atoms)
q_all = sum(nb.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
            for i in range(system.getNumParticles()))
print("\natoms in system:", system.getNumParticles())
print("ligand total charge: %.2f   (expected -2.00)" % q_lig)
print("system total charge: %.2f   (expected 0.00)" % q_all)

app.PDBFile.writeFile(top, modeller.positions, open(os.path.join(OUT, "solvated.pdb"), "w"))
open(os.path.join(OUT, "system.xml"), "w").write(mm.XmlSerializer.serialize(system))
json.dump([[a, b, c, d] for a, b, c, d in restr], open(os.path.join(OUT, "restraints.json"), "w"))

# 5) minimize + speed test on the Quadro T2000 (OpenCL platform 0, device 0)
integ = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
plat = mm.Platform.getPlatformByName("OpenCL")
sim = app.Simulation(top, system, integ, plat, {"OpenCLPlatformIndex": "0", "DeviceIndex": "0"})
sim.context.setPositions(modeller.positions)
print("\nminimizing ...")
sim.minimizeEnergy(maxIterations=1000)
sim.context.setVelocitiesToTemperature(300 * unit.kelvin)
sim.step(500)
t0 = time.time()
sim.step(2000)
sim.context.getState(getEnergy=True)
dt = time.time() - t0
print("speed: %.1f ns/day  (2 fs steps, %d atoms)" % (2000 * 0.002 / 1000 / dt * 86400, system.getNumParticles()))

st = sim.context.getState(getPositions=True, getEnergy=True)
p = xyz(st.getPositions())
print("\nDistances after minimization + 1 ps (Angstrom):")
for zn, j, d, label in restr:
    print("  %-22s crystal %.2f -> now %.2f" % (label, d * 10, float(np.linalg.norm(p[zn] - p[j])) * 10))
app.PDBFile.writeFile(top, st.getPositions(), open(os.path.join(OUT, "minimized.pdb"), "w"))
print("\nSaved in folder", OUT, ": solvated.pdb, minimized.pdb, system.xml, restraints.json, cap.xml")