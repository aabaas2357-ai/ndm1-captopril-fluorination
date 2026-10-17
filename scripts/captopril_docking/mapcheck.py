import numpy as np
from rdkit import Chem
cen = np.array([30.5, 40.5, 26.7]); sp = 0.375; n = 67
vals = np.array([float(x) for x in open("protein_tz.SA.map").read().split("\n")[6:] if x.strip()])
g = vals.reshape((n, n, n))  # x fastest in AutoGrid files
m = Chem.MolFromMolFile("captopril_xtal_SH.sdf", removeHs=False)
p = [m.GetConformer().GetAtomPosition(a.GetIdx()) for a in m.GetAtoms() if a.GetSymbol() == "S"][0]
idx = np.rint((np.array([p.x, p.y, p.z]) - cen) / sp).astype(int) + n // 2
print("SA map energy at crystal S:", g[idx[2], idx[1], idx[0]])
print("map min / median / max:", vals.min(), np.median(vals), vals.max())
