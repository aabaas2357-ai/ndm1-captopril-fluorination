import csv
rows = list(csv.DictReader(open("library_scores3.csv")))
rows = [r for r in rows if r["local_only"]]
names = [r["id"] for r in rows]
cols = {"ref_Dcap": "green", "ref_CF2_44diF": "magenta", "ref_CF3chain": "cyan", "ref_CF3_44diF": "teal"}
extra = ["yellow", "salmon", "orange", "slate", "violet", "lime", "hotpink", "wheat"]
L = ["reinitialize", "load target_5ZJ2.pdb, prot", "hide everything, prot", "show cartoon, prot",
     "color gray80, prot", "set cartoon_transparency, 0.5, prot",
     "select zn, prot and elem Zn", "show spheres, zn", "color orange, zn", "set sphere_scale, 0.45, zn",
     "select key, prot and resi 67+70+73+120+122+124+189+208+211+220+250",
     "show sticks, key", "color palecyan, key and elem C", "set stick_radius, 0.14",
     "label key and name CA, resn+resi", "set label_size, -0.6", "bg_color white"]
ci = 0
for r in rows:
    n = r["id"]
    c = cols.get(n)
    if c is None:
        c = extra[ci % len(extra)]
        ci += 1
    L += [f"load {n}_local.sdf, {n}", f"show sticks, {n}", f"color {c}, {n} and elem C",
          f"set stick_radius, 0.2, {n}", f"dist {n}_SZn, ({n} and elem S), zn, 3.5", f"hide labels, {n}_SZn"]
L.append("disable all")
L += ["enable prot", "enable zn", "enable key"]
for r in rows:
    n = r["id"]
    L += ["disable " + " ".join(f"{x} {x}_SZn" for x in names),
          "enable ref_Dcap", "enable ref_Dcap_SZn", f"enable {n}", f"enable {n}_SZn",
          f"zoom ({n} or zn), 9",
          f"scene {n}, store, \"{n} | local {r['local_only']} | S-Zn {r['S_Zn_min']} | contact {r['min_contact']}\""]
L += ["set scene_buttons, 1", "scene ref_Dcap, recall"]
open("view_all.pml", "w").write("\n".join(L) + "\n")
print("wrote view_all.pml for", len(rows), "ligands")
