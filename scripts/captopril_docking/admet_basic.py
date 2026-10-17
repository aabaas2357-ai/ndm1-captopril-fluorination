import csv
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

cmpds = {
    "D-captopril": "C[C@H](CS)C(=O)N1CCC[C@@H]1C(=O)O",
    "4,4-difluoro-D-cap": "C[C@H](CS)C(=O)N1CC(F)(F)C[C@@H]1C(=O)O",
}

def cat(c):
    p = FilterCatalogParams()
    p.AddCatalog(c)
    return FilterCatalog(p)

pains = cat(FilterCatalogParams.FilterCatalogs.PAINS)
brenk = cat(FilterCatalogParams.FilterCatalogs.BRENK)

def alerts(c, m):
    return "; ".join(e.GetDescription() for e in c.GetMatches(m)) or "none"

rows = []
for name, smi in cmpds.items():
    m = Chem.MolFromSmiles(smi)
    mw, lp = Descriptors.MolWt(m), Crippen.MolLogP(m)
    hbd, hba = Lipinski.NumHDonors(m), Lipinski.NumHAcceptors(m)
    r = {
        "name": name,
        "MW": round(mw, 2),
        "cLogP": round(lp, 2),
        "TPSA": round(rdMolDescriptors.CalcTPSA(m), 1),
        "HBD": hbd,
        "HBA": hba,
        "RotBonds": Lipinski.NumRotatableBonds(m),
        "Fsp3": round(rdMolDescriptors.CalcFractionCSP3(m), 2),
        "QED": round(QED.qed(m), 2),
        "Lipinski_violations": int(mw > 500) + int(lp > 5) + int(hbd > 5) + int(hba > 10),
        "PAINS": alerts(pains, m),
        "Brenk": alerts(brenk, m),
    }
    rows.append(r)
    print(r)

with open("admet_basic.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)
print("saved admet_basic.csv")
