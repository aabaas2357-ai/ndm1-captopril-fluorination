"""
md_run.py - restrained equilibration + production MD for the NDM-1 / D-captopril (5ZJ2) system.

Usage:
  python md_run.py [production_ns] [K_ZN] [--rep NAME] [--seed N] [--resume] [--overwrite]
  python md_run.py 2.0 40000                       new run
  python md_run.py 2.0 40000 --rep r2 --seed 22    independent replicate (its own files)
  python md_run.py 2.0 40000 --resume              continue after a crash / power cut
                                                   (same arguments as the original run)

Needs md_5ZJ2/minimized.pdb, system.xml and restraints.json written by md_prep.py.
Production saves a checkpoint every 50 ps (ckpt_<tag>.chk + progress_<tag>.json).
If the run is interrupted, --resume continues from the last checkpoint and writes the new
frames to prod_<tag>_partN.dcd (progress_<tag>.json lists which time range of each part is valid).

BASELINE setup: ligand atom types borrowed from amber14, Gasteiger charges, non-bonded Zn ions
with harmonic coordination restraints. Not validated; results are qualitative only.
"""
import sys, os, json, time, csv, argparse
import numpy as np
import openmm as mm
from openmm import app, unit

ap = argparse.ArgumentParser()
ap.add_argument("prod_ns", nargs="?", type=float, default=1.0)
ap.add_argument("k_zn", nargs="?", type=float, default=10000.0)   # kJ/mol/nm^2
ap.add_argument("--rep", default="")
ap.add_argument("--seed", type=int, default=None)
ap.add_argument("--resume", action="store_true")
ap.add_argument("--overwrite", action="store_true")
ap.add_argument("--dir", default="md_5ZJ2")
args = ap.parse_args()

OUT = args.dir
PROD_NS = args.prod_ns
K_ZN = args.k_zn
TAG = "K%d" % int(K_ZN) + ("_" + args.rep if args.rep else "")
CKPT = "%s/ckpt_%s.chk" % (OUT, TAG)
PROG = "%s/progress_%s.json" % (OUT, TAG)
SOLV = {"HOH", "WAT", "NA", "CL", "Na+", "Cl-"}
DT_FS = 2.0
STEPS_PER_PS = int(round(1000.0 / DT_FS))
CHUNK_PS = 5
CKPT_EVERY_PS = 50
EQ_STAGES = [(1000.0, 50), (300.0, 50), (100.0, 50), (30.0, 50), (0.0, 50)]
EQ_PS = float(sum(ps for _, ps in EQ_STAGES))

if not args.resume and os.path.exists(PROG) and not args.overwrite:
    raise SystemExit("files for tag %s already exist. Use --resume to continue it, --rep NAME for a new "
                     "replicate, or --overwrite to start over." % TAG)


def val(x, u):
    return float(x.value_in_unit(u)) if hasattr(x, "value_in_unit") else float(x)


def kabsch(P, Q):
    """rotation R and centroids so that (X - Pc) @ R.T + Qc maps frame P onto frame Q"""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, Pc, Qc


pdb = app.PDBFile(OUT + "/minimized.pdb")
top = pdb.topology
pos0 = pdb.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
system = mm.XmlSerializer.deserialize(open(OUT + "/system.xml").read())
restr = json.load(open(OUT + "/restraints.json"))   # [zn_index, partner_index, d_nm, label]

# 1) Zn coordination restraints: set the force constant requested on the command line
bfs = [f for f in system.getForces()
       if isinstance(f, mm.HarmonicBondForce) and f.getNumBonds() == len(restr)]
if len(bfs) != 1:
    raise SystemExit("could not identify the zinc restraint force (found %d candidates)" % len(bfs))
bf = bfs[0]
a0, b0 = bf.getBondParameters(0)[0], bf.getBondParameters(0)[1]
if (a0, b0) != (restr[0][0], restr[0][1]):
    raise SystemExit("restraint force does not match restraints.json")
for k in range(bf.getNumBonds()):
    a, b, length, _ = bf.getBondParameters(k)
    bf.setBondParameters(k, a, b, val(length, unit.nanometer), K_ZN)
print("zinc restraint force constant: %.0f kJ/mol/nm^2" % K_ZN)

# 2) report the Zn non-bonded parameters actually used
nb = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
for zi in sorted({r[0] for r in restr}):
    q, sg, ep = nb.getParticleParameters(zi)
    print("Zn particle %d: charge %.2f, sigma %.4f nm, epsilon %.5f kJ/mol"
          % (zi, val(q, unit.elementary_charge), val(sg, unit.nanometer), val(ep, unit.kilojoule_per_mole)))

# 3) heavy-atom position restraints (protein + ligand + Zn), released stepwise
pr = mm.CustomExternalForce("0.5*kpos*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
pr.addGlobalParameter("kpos", 1000.0)
for nm in ("x0", "y0", "z0"):
    pr.addPerParticleParameter(nm)
nres = 0
for a in top.atoms():
    if a.element is not None and a.element.symbol != "H" and a.residue.name not in SOLV:
        p = pos0[a.index]
        pr.addParticle(a.index, [float(p[0]), float(p[1]), float(p[2])])
        nres += 1
system.addForce(pr)
print("position-restrained heavy atoms:", nres)

# 4) NPT
baro = mm.MonteCarloBarostat(1.0 * unit.bar, 300.0 * unit.kelvin, 25)
integ = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, DT_FS * unit.femtoseconds)
if args.seed is not None:
    baro.setRandomNumberSeed(args.seed)
    integ.setRandomNumberSeed(args.seed)
system.addForce(baro)
plat = mm.Platform.getPlatformByName("OpenCL")
sim = app.Simulation(top, system, integ, plat, {"OpenCLPlatformIndex": "0", "DeviceIndex": "0"})
sim.context.setPositions(pdb.positions)
if args.seed is not None:
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, args.seed)
    print("random seed:", args.seed)
else:
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin)

ca_idx = [a.index for a in top.atoms() if a.name == "CA" and a.residue.name != "CAP"]
lig_idx = [a.index for a in top.atoms() if a.residue.name == "CAP" and a.element.symbol != "H"]
print("CA atoms for alignment: %d, ligand heavy atoms: %d" % (len(ca_idx), len(lig_idx)))

clock = [0.0]
resuming = args.resume
prog = None
if resuming:
    if not (os.path.exists(CKPT) and os.path.exists(PROG)):
        raise SystemExit("nothing to resume for tag %s (no checkpoint saved yet)" % TAG)
    prog = json.load(open(PROG))
    if prog.get("done"):
        raise SystemExit("run %s already finished" % TAG)
    sim.loadCheckpoint(CKPT)
    sim.context.setParameter("kpos", 0.0)
    clock[0] = float(prog["clock_ps"])
    print("RESUMING %s from t = %.1f ps" % (TAG, clock[0]))

header = ["stage", "time_ps", "ligand_RMSD_A"] + [r[3] + "_A" for r in restr]
trk = "%s/track_%s.csv" % (OUT, TAG)
keep = []
if resuming and os.path.exists(trk):
    with open(trk, newline="") as f:
        rd = csv.reader(f)
        next(rd, None)
        for r in rd:
            try:
                if len(r) == len(header) and float(r[1]) <= clock[0] + 1e-6:
                    keep.append(r)
            except ValueError:
                pass
fh = open(trk, "w", newline="")
w = csv.writer(fh)
w.writerow(header)
w.writerows(keep)
fh.flush()

S_IDX = [i for i, r in enumerate(restr) if "S" in r[3] and "lig" in r[3].lower()]
if len(S_IDX) != 2:
    S_IDX = [3, 7]
print("S-Zn restraint indices:", S_IDX)

PARTS = []


def save_progress(done=False):
    """flush trajectory/log files, then write checkpoint + progress atomically"""
    for rp in sim.reporters:
        o = getattr(rp, "_out", None)
        if o is not None:
            try:
                o.flush()
                os.fsync(o.fileno())
            except Exception:
                pass
    sim.saveCheckpoint(CKPT + ".tmp")
    os.replace(CKPT + ".tmp", CKPT)
    if PARTS:
        PARTS[-1]["end_ps"] = clock[0]
    with open(PROG + ".tmp", "w") as f:
        json.dump({"clock_ps": clock[0], "parts": PARTS, "prod_ns": PROD_NS, "eq_ps": EQ_PS,
                   "k_zn": K_ZN, "done": done}, f, indent=1)
    os.replace(PROG + ".tmp", PROG)


def observe(stage):
    st = sim.context.getState(getPositions=True)
    p = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    R, Pc, Qc = kabsch(p[ca_idx], pos0[ca_idx])
    L = (p[lig_idx] - Pc) @ R.T + Qc
    rmsd = float(np.sqrt(((L - pos0[lig_idx]) ** 2).sum(1).mean()) * 10)
    d = [float(np.linalg.norm(p[a] - p[b]) * 10) for a, b, _, _ in restr]
    row = [stage, round(clock[0], 1), round(rmsd, 2)] + [round(x, 2) for x in d]
    w.writerow(row)
    fh.flush()
    return row


def run(stage, ps, verbose_every=100, checkpoint=False):
    row = None
    ck_chunks = max(1, int(CKPT_EVERY_PS // CHUNK_PS))
    for i in range(int(round(ps / CHUNK_PS))):
        sim.step(CHUNK_PS * STEPS_PER_PS)
        clock[0] += CHUNK_PS
        row = observe(stage)
        if (i + 1) % max(1, verbose_every // CHUNK_PS) == 0:
            print("  t=%7.1f ps  ligand RMSD %.2f A  S-Zn: %.2f %.2f"
                  % (clock[0], row[2], row[3 + S_IDX[0]], row[3 + S_IDX[1]]), flush=True)
        if checkpoint and (i + 1) % ck_chunks == 0:
            save_progress()
    return row


t0 = time.time()
if not resuming:
    print("\nequilibration (restraints released stepwise)")
    for k, ps in EQ_STAGES:
        sim.context.setParameter("kpos", k)
        row = run("equil_k%d" % int(k), ps, verbose_every=50)
        print("stage k=%g done: ligand RMSD %.2f A" % (k, row[2]), flush=True)
    state = sim.context.getState(getPositions=True)
    app.PDBFile.writeFile(top, state.getPositions(), open("%s/prod_start_%s.pdb" % (OUT, TAG), "w"))
    part_n = 1
    remaining_ps = PROD_NS * 1000.0
else:
    PARTS[:] = prog["parts"]
    part_n = len(PARTS) + 1
    remaining_ps = PROD_NS * 1000.0 - (clock[0] - EQ_PS)

suffix = "" if part_n == 1 else "_part%d" % part_n
dcd_name = "prod_%s%s.dcd" % (TAG, suffix)
log_name = "prod_%s%s.log" % (TAG, suffix)
PARTS.append({"dcd": dcd_name, "log": log_name, "start_ps": clock[0], "end_ps": clock[0]})
sim.reporters.append(app.DCDReporter("%s/%s" % (OUT, dcd_name), 5000))
sim.reporters.append(app.StateDataReporter("%s/%s" % (OUT, log_name), 5000, step=True, time=True,
                                          potentialEnergy=True, temperature=True, density=True,
                                          speed=True))
save_progress()
print("\nproduction: %.2f ns total, %.1f ps still to run, frame every 10 ps, checkpoint every %d ps"
      % (PROD_NS, remaining_ps, CKPT_EVERY_PS))
if remaining_ps > 0:
    run("prod", remaining_ps, verbose_every=100, checkpoint=True)
save_progress(done=True)
app.PDBFile.writeFile(top, sim.context.getState(getPositions=True).getPositions(),
                      open("%s/final_%s.pdb" % (OUT, TAG), "w"))
fh.close()
print("\nwall time: %.1f min" % ((time.time() - t0) / 60.0))
print("files in %s: track_%s.csv, %s, %s, prod_start_%s.pdb, final_%s.pdb, ckpt_%s.chk, progress_%s.json"
      % (OUT, TAG, dcd_name, log_name, TAG, TAG, TAG, TAG))