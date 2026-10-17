import csv
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, rdCIPLabeler
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

LIB = [
 ("ref_Dcap", "D-captopril", "C[C@H](CS)C(=O)N1CCC[C@@H]1C(=O)O", "ref"),
 ("ref_CF3chain", "CF3-chain (5aC)", "FC(F)(F)[C@H](CS)C(=O)N1CCC[C@@H]1C(=O)O", "ref"),
 ("ref_Lcap", "L-captopril", "C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O", "ref"),
 ("ref_epiDcap", "epi-D-captopril", "C[C@@H](CS)C(=O)N1CCC[C@@H]1C(=O)O", "ref"),
 ("t1_44diF", "4,4-diF", "C[C@H](CS)C(=O)N1CC(F)(F)C[C@@H]1C(=O)O", "1"),
 ("t1_CF3_44diF", "CF3-chain + 4,4-diF", "FC(F)(F)[C@H](CS)C(=O)N1CC(F)(F)C[C@@H]1C(=O)O", "1"),
 ("t2_33diF", "3,3-diF", "C[C@H](CS)C(=O)N1CCC(F)(F)[C@@H]1C(=O)O", "2"),
 ("t2_4F", "4-F", "C[C@H](CS)C(=O)N1CC(F)C[C@@H]1C(=O)O", "2"),
 ("t2_3F", "3-F", "C[C@H](CS)C(=O)N1CCC(F)[C@@H]1C(=O)O", "2"),
 ("t2_CHF2chain", "CHF2-chain", "FC(F)[C@H](CS)C(=O)N1CCC[C@@H]1C(=O)O", "2"),
 ("t3_pip", "pipecolic D", "C[C@H](CS)C(=O)N1CCCC[C@@H]1C(=O)O", "3"),
 ("t3_isonip", "isonipecotic", "C[C@H](CS)C(=O)N1CCC(CC1)C(=O)O", "3"),
 ("t3_CF3_isonip", "CF3 + isonipecotic", "FC(F)(F)[C@H](CS)C(=O)N1CCC(CC1)C(=O)O", "3"),
 ("t3_azet", "azetidine D", "C[C@H](CS)C(=O)N1CC[C@@H]1C(=O)O", "3"),
]

rows, w3d = [], Chem.SDWriter("library_3d.sdf")
for cid, name, smi, tier in LIB:
    m = Chem.MolFromSmiles(smi)
    isos = list(EnumerateStereoisomers(m))
    for k, iso in enumerate(isos):
        rdCIPLabeler.AssignCIPLabels(iso)
        cip = ",".join(a.GetProp("_CIPCode") for a in iso.GetAtoms() if a.HasProp("_CIPCode"))
        uid = cid if len(isos) == 1 else f"{cid}_iso{k+1}"
        rows.append({"id": uid, "name": name, "tier": tier,
                     "smiles": Chem.MolToSmiles(iso), "CIP(atom order)": cip,
                     "MW": round(Descriptors.MolWt(iso), 1),
                     "cLogP": round(Crippen.MolLogP(iso), 2)})
        mh = Chem.AddHs(iso)
        if AllChem.EmbedMolecule(mh, randomSeed=7) == 0:
            AllChem.MMFFOptimizeMolecule(mh)
            mh.SetProp("_Name", uid)
            w3d.write(mh)
w3d.close()
for r in rows:
    print(r["id"], r["CIP(atom order)"], r["MW"], r["cLogP"])
with open("library.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=rows[0].keys())
    wr.writeheader()
    wr.writerows(rows)
print("saved library.csv and library_3d.sdf:", len(rows), "compounds")
