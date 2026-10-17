import numpy as np
a = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in open("mco_crystal.pdb") if l.startswith(("ATOM", "HETATM"))])
print("%.3f,%.3f,%.3f" % tuple(a.mean(0)))
