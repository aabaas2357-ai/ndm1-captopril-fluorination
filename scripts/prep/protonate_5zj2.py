import math, subprocess, sys
lines = open("5ZJ2_clean.pdb").read().splitlines()
zn = [l for l in lines if l.startswith("HETATM") and l[17:20].strip() == "ZN"]
prot = [l for l in lines if l.startswith("ATOM")]
xyz = lambda l: (float(l[30:38]), float(l[38:46]), float(l[46:54]))
his, cys = {}, set()
for z in zn:
    for a in prot:
        d = math.dist(xyz(z), xyz(a))
        if a[76:78].strip() in ("N", "O", "S") and d < 2.7:
            res, num, name = a[17:20], int(a[22:26]), a[12:16].strip()
            print(f"Zn{z[22:26].strip()}: {res}{num} {name} {d:.2f} A")
            if res == "HIS": his[num] = "HID" if name == "NE2" else "HIE"
            if res == "CYS" and name == "SG": cys.add(num)
print("His states:", his, "| Cys thiolate:", sorted(cys))
out = []
for l in prot:
    if l[17:20] == "HIS" and int(l[22:26]) in his:
        l = l[:17] + his[int(l[22:26])] + l[20:]
    out.append(l)
open("5ZJ2_noZn.pdb", "w").write("\n".join(out) + "\nEND\n")
subprocess.run([sys.executable, "-m", "pdb2pqr", "--ff", "AMBER", "--ffout", "AMBER", "--keep-chain",
                "--noopt", "--pdb-output", "5ZJ2_Hraw.pdb", "5ZJ2_noZn.pdb", "5ZJ2_H.pqr"], check=True)
final = []
for l in open("5ZJ2_Hraw.pdb"):
    if l.startswith(("ATOM", "HETATM")):
        if l[12:16].strip() == "HG" and int(l[22:26]) in cys: continue
        if l[12:16].strip() in ("HD1", "HE2") and int(l[22:26]) in his:
            print("His", int(l[22:26]), l[17:20], "has", l[12:16].strip())
        final.append(l.rstrip("\n"))
open("5ZJ2_H.pdb", "w").write("\n".join(final + zn) + "\nEND\n")
print("written 5ZJ2_H.pdb:", len(final), "atoms + 2 Zn")
