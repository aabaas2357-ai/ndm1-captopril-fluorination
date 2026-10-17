from rdkit import Chem
m = Chem.MolFromMolFile(r"xtal\captopril_xtal_SH.sdf")
for idx, lab in Chem.FindMolChiralCenters(m, useLegacyImplementation=False):
    print(idx, lab, "ring (proline alpha)" if m.GetAtomWithIdx(idx).IsInRing() else "chain (methyl-bearing)")
