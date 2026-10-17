from rdkit import Chem

sdf = r".\xtal\F2captopril_xtal_SH.sdf"
mol = Chem.MolFromMolFile(sdf, removeHs=False)
print("SMILES:", Chem.MolToSmiles(Chem.RemoveHs(mol)))

# proline ring: N1, C2 (carries COOH), C3, C4, C5
pat = Chem.MolFromSmarts("[NX3;R]1[CX4;R]([CX3](~O)~O)[CX4;R][CX4;R][CX4;R]1")
m = mol.GetSubstructMatch(pat)
if not m:
    print("proline ring pattern not found")
else:
    names = {m[0]: "N1", m[1]: "C2 (alpha, COOH)", m[5]: "C3", m[6]: "C4", m[7]: "C5"}
    for a in mol.GetAtoms():
        if a.GetSymbol() == "F":
            c = [n for n in a.GetNeighbors()][0]
            print(f"F idx {a.GetIdx()} is on carbon idx {c.GetIdx()} ->", names.get(c.GetIdx(), "NOT a proline ring carbon"))
