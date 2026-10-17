import sys
src, LIG, CH = "5ZJ2.pdb", (sys.argv[1] if len(sys.argv) > 1 else "MCO"), (sys.argv[2] if len(sys.argv) > 2 else "A")
prot, lig, chains, hets = [], [], set(), set()
for ln in open(src):
    if not ln.startswith(("ATOM", "HETATM")): continue
    if ln[16] not in (" ", "A"): continue
    ln = ln[:16] + " " + ln[17:]
    res, ch = ln[17:20].strip(), ln[21]
    chains.add(ch)
    if ln.startswith("HETATM"): hets.add(res)
    if ch != CH: continue
    if ln.startswith("ATOM") or res == "ZN": prot.append(ln)
    elif res == LIG: lig.append(ln)
open("5ZJ2_clean.pdb", "w").write("".join(prot) + "END\n")
open("mco_crystal.pdb", "w").write("".join(lig) + "END\n")
print("chains:", sorted(chains), "| hetero residues:", sorted(hets))
print("protein+Zn lines:", len(prot), "| ligand lines:", len(lig))
print("Zn lines:", sum(1 for l in prot if l[17:20].strip() == "ZN"))
