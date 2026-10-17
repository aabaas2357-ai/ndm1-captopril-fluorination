"""
build_ligand.py - build the D-captopril dianion in the CRYSTAL pose (mco_crystal.pdb) and parametrize it
with GAFF2 + AM1-BCC charges (antechamber / sqm). Run inside WSL with the 'amber' conda env python.
Outputs go to C:\\NDM1_Project\\md\\amber_lig
"""
import os, sys, shutil, subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
import parmed

SRC = "/mnt/c/NDM1_Project/md/mco_crystal.pdb"
OUTC = "/mnt/c/NDM1_Project/md/amber_lig"
WORK = os.path.expanduser("~/amber_work/lig")
os.makedirs(WORK, exist_ok=True)
os.makedirs(OUTC, exist_ok=True)
os.chdir(WORK)
os.environ["AMBERHOME"] = sys.prefix
os.environ["PATH"] = os.path.join(sys.prefix, "bin") + os.pathsep + os.environ.get("PATH", "")


def sh(cmd):
    print("$", " ".join(cmd), flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    out = (r.stdout + r.stderr).strip()
    if out:
        print("\n".join(out.splitlines()[-15:]))
    if r.returncode != 0:
        raise SystemExit("command failed (exit code %d)" % r.returncode)


# 1) heavy atoms from the crystal + bond orders from a neutral template
tmpl = Chem.MolFromSmiles("CC(CS)C(=O)N1CCCC1C(=O)O")
xt = Chem.MolFromPDBFile(SRC, removeHs=True, proximityBonding=True)
if xt is None:
    raise SystemExit("could not read " + SRC)
print("crystal ligand heavy atoms:", xt.GetNumAtoms(), "(expected 14)")
mol = AllChem.AssignBondOrdersFromTemplate(tmpl, xt)

# 2) thiolate and carboxylate (net charge -2)
for a in mol.GetAtoms():
    if a.GetSymbol() == "S":
        a.SetFormalCharge(-1)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
    elif a.GetSymbol() == "O" and a.GetDegree() == 1 and a.GetBonds()[0].GetBondType() == Chem.BondType.SINGLE:
        a.SetFormalCharge(-1)
        a.SetNoImplicit(True)
        a.SetNumExplicitHs(0)
Chem.SanitizeMol(mol)

# 3) stereo straight from the crystal coordinates
Chem.AssignStereochemistryFrom3D(mol)
cc = Chem.FindMolChiralCenters(mol, useLegacyImplementation=False)
print("chiral centers (index, label, in ring):", [(i, l, mol.GetAtomWithIdx(i).IsInRing()) for i, l in cc])
ok = len(cc) == 2 and all((l == "R") if mol.GetAtomWithIdx(i).IsInRing() else (l == "S") for i, l in cc)
print("stereo check vs D-captopril (ring R, chain S):", "PASS" if ok else "FAIL")
if not ok:
    raise SystemExit("stereo does not match D-captopril - stop here and send me this output")

# 4) hydrogens, relaxed with heavy atoms fixed at the crystal positions
mh = Chem.AddHs(mol, addCoords=True)
print("formula:", rdMolDescriptors.CalcMolFormula(mh), " formal charge:", Chem.GetFormalCharge(mh))
try:
    if AllChem.UFFHasAllMoleculeParams(mh):
        ff = AllChem.UFFGetMoleculeForceField(mh)
        for a in mh.GetAtoms():
            if a.GetAtomicNum() > 1:
                ff.AddFixedPoint(a.GetIdx())
        ff.Initialize()
        ff.Minimize(maxIts=5000)
        print("hydrogens relaxed (UFF, heavy atoms fixed)")
    else:
        print("UFF parameters missing: hydrogens left as built")
except Exception as e:
    print("hydrogen relaxation skipped:", e)
Chem.MolToMolFile(mh, "lig_xtal.mol")

# 5) GAFF2 types + AM1-BCC charges
sh(["antechamber", "-i", "lig_xtal.mol", "-fi", "mdl", "-o", "cap_gaff2.mol2", "-fo", "mol2",
    "-c", "bcc", "-nc", "-2", "-at", "gaff2", "-rn", "CAP", "-s", "2", "-pf", "n"])
sh(["parmchk2", "-i", "cap_gaff2.mol2", "-f", "mol2", "-o", "cap_gaff2.frcmod", "-s", "gaff2"])

# 6) put the crystal-pose coordinates back (antechamber optimizes the geometry for the charges)
s = parmed.load_file("cap_gaff2.mol2", structure=True)
xyz = mh.GetConformer().GetPositions()
if len(s.atoms) == xyz.shape[0]:
    s.coordinates = xyz
    s.save("cap_gaff2_xtal.mol2", overwrite=True)
else:
    print("WARNING: atom count differs (%d vs %d), coordinates not restored" % (len(s.atoms), xyz.shape[0]))

# 7) report
AllChem.ComputeGasteigerCharges(mh)
print("\nidx el  name  gaff2   AM1-BCC   Gasteiger(old)")
tot = 0.0
for i, at in enumerate(s.atoms):
    g = mh.GetAtomWithIdx(i).GetDoubleProp("_GasteigerCharge") if mh.GetAtomWithIdx(i).HasProp("_GasteigerCharge") else float("nan")
    tot += at.charge
    print("%3d %-3s %-5s %-6s %8.4f %10.4f" % (i, mh.GetAtomWithIdx(i).GetSymbol(), at.name, at.type, at.charge, g))
print("sum of AM1-BCC charges: %.4f (expected -2.0000)" % tot)

for f in ("lig_xtal.mol", "cap_gaff2.mol2", "cap_gaff2.frcmod", "cap_gaff2_xtal.mol2", "sqm.out"):
    if os.path.exists(f):
        shutil.copy(f, OUTC)
print("\nfiles copied to", OUTC)