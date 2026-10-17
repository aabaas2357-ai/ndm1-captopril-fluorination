import collections
c = collections.Counter(l[16:20].strip() for l in open("5ZJ2_H.pdb") if l.startswith(("ATOM", "HETATM")))
print(sorted(c.items()))
