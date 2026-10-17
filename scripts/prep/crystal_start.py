import os, subprocess, sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

os.makedirs("xtal", exist_ok=True)
tmpl = Chem.MolFromSmiles("C[C@H](CS)C(=O)N1CCC[C@@H]1C(=O)O")

def base(thiolate):
    ref = Chem.MolFromPDBFile("mco_crystal.pdb")
    m = AllChem.AssignBondOrdersFromTemplate(tmpl, ref)
    rw = Chem.RWMol(m)
    for hit in rw.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[OX2H1]")):
        a = rw.GetAtomWithIdx(hit[2]); a.SetFormalCharge(-1); a.SetNumExplicitHs(0); a.SetNoImplicit(True)
    if thiolate:
        for hit in rw.GetSubstructMatches(Chem.MolFromSmarts("[SX2H1]")):
            a = rw.GetAtomWithIdx(hit[0]); a.SetFormalCharge(-1); a.SetNumExplicitHs(0); a.SetNoImplicit(True)
    m = rw.GetMol(); Chem.SanitizeMol(m)
    m = Chem.AddHs(m, addCoords=True)
    Chem.AssignStereochemistryFrom3D(m)
    return m

def to_f2(m):
    m = Chem.RWMol(m)
    hit = m.GetSubstructMatch(Chem.MolFromSmarts("N1[CH2][CH2][CH2][CH]1C(=O)O"))
    c4 = m.GetAtomWithIdx(hit[2])
    conf = m.GetConformer()
    pc = np.array(conf.GetAtomPosition(c4.GetIdx()))
    hs = [n for n in c4.GetNeighbors() if n.GetAtomicNum() == 1]
    assert len(hs) == 2, "C4 does not carry two H"
    for h in hs:
        ph = np.array(conf.GetAtomPosition(h.GetIdx()))
        v = (ph - pc) / np.linalg.norm(ph - pc)
        h.SetAtomicNum(9)
        conf.SetAtomPosition(h.GetIdx(), (pc + 1.36 * v).tolist())
    m = m.GetMol(); Chem.SanitizeMol(m)
    return m

def meeko(sdf, pdbqt):
    subprocess.run([sys.executable, "-m", "meeko.cli.mk_prepare_ligand", "-i", sdf, "-o", pdbqt], check=True)
    lines = []
    for ln in open(pdbqt):
        ln = ln.rstrip("\n")
        if ln.startswith(("ATOM", "HETATM")) and ln[77:79].strip() == "S": ln = ln[:77] + "SA"
        lines.append(ln)
    open(pdbqt, "w").write("\n".join(lines) + "\n")

for tag, thiolate in (("SH", False), ("Sm", True)):
    m = base(thiolate)
    print("crystal ligand CIP (expect methyl C S, proline C R):", Chem.FindMolChiralCenters(m, useLegacyImplementation=False))
    for name, mol in (("captopril", m), ("F2captopril", to_f2(m))):
        sdf = f"xtal\\{name}_xtal_{tag}.sdf"
        Chem.MolToMolFile(mol, sdf)
        meeko(sdf, sdf.replace(".sdf", ".pdbqt"))
        print(name, "heavy atoms:", mol.GetNumHeavyAtoms(), "net charge:", Chem.GetFormalCharge(mol))
