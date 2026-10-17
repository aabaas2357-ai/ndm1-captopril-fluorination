from pymol import cmd
import math, os

D = r"."
LOG = []
def out(s=""):
    print(s)
    LOG.append(s)

cmd.reinitialize()
cmd.load(D + r"\5ZJ2_H.pdb", "prot")
cmd.load(D + r"\mco_crystal.pdb", "lig")
cmd.load(D + r"\xtal\F2captopril_xtal_Sm.sdf", "F2cap")
orig_path = D + r"\5ZJ2.pdb"
if os.path.exists(orig_path):
    cmd.load(orig_path, "orig")
else:
    cmd.fetch("5ZJ2", name="orig", type="pdb")

# scene
cmd.bg_color("white")
cmd.hide("everything")
cmd.show("cartoon", "prot")
cmd.color("gray80", "prot")
cmd.set("cartoon_transparency", 0.5)
cmd.select("zns", "prot and elem Zn")
cmd.show("spheres", "zns")
cmd.set("sphere_scale", 0.5, "zns")
cmd.color("purpleblue", "zns")
cmd.select("coord", "prot and chain A and resi 120+122+124+189+208+250")
cmd.show("sticks", "coord")
cmd.color("marine", "coord and elem C")
cmd.show("sticks", "lig")
cmd.color("orange", "lig and elem C")
cmd.show("sticks", "F2cap")
cmd.color("green", "F2cap and elem C")
cmd.set("stick_radius", 0.2)
cmd.select("Fatoms", "F2cap and elem F")
cmd.show("spheres", "Fatoms")
cmd.set("sphere_scale", 0.3, "Fatoms")
cmd.color("cyan", "Fatoms")
cmd.disable("orig")
cmd.distance("zn_S", "zns", "lig and elem S", 3.2)
cmd.distance("zn_coord", "zns", "coord and (elem N+O+S)", 2.6)
cmd.set("dash_color", "black")
cmd.set("dash_width", 3)
cmd.hide("labels")
cmd.select("pocket", "byres (prot within 5 of lig)")
cmd.show("lines", "pocket")
cmd.orient("lig")
cmd.zoom("lig", 8)

# sanity: counts and frame check
out("F atoms in F2cap: %d" % cmd.count_atoms("Fatoms"))
out("Zn in prot: %d | Zn in orig: %d" % (cmd.count_atoms("prot and elem Zn"), cmd.count_atoms("orig and elem Zn")))
out("waters in orig: %d" % cmd.count_atoms("orig and resn HOH"))
Zp = cmd.get_model("prot and elem Zn").atom
Zo = cmd.get_model("orig and elem Zn").atom
out("frame check (orig Zn -> nearest prot Zn, should be ~0):")
for a in Zo:
    out("   %s%s  %.2f A" % (a.resn, a.resi, min(math.dist(a.coord, b.coord) for b in Zp)))
Fm = cmd.get_model("Fatoms").atom

# crystal ligand C6 vs Val73 CG1
try:
    out("lig C6 - Val73 CG1: %.2f A" % cmd.get_distance("lig and name C6", "prot and resn VAL and resi 73 and name CG1"))
except Exception as e:
    out("C6-Val73 distance not available: %s" % e)

# F2cap thiolate S - Zn
out("")
for s in cmd.get_model("F2cap and elem S").atom:
    out("F2cap S-Zn: " + ", ".join("%.2f" % math.dist(s.coord, z.coord) for z in Zp) + " A")

# protein atoms near F
cmd.select("Fnear", "(prot and not elem H) within 4.5 of Fatoms")
N = cmd.get_model("Fnear").atom
for f in Fm:
    out("")
    out("F atom index %d: closest non-H protein/Zn atoms" % f.index)
    rows = sorted((math.dist(f.coord, n.coord), n.resn, n.resi, n.chain, n.name) for n in N)
    for r in rows[:12]:
        flag = "  <-- CLOSE (<2.9)" if r[0] < 2.9 else ""
        out("   %.2f A  %s%s:%s %s%s" % (r + (flag,)))

# waters in the crystal structure near F
cmd.select("Wnear", "(orig and resn HOH) within 5.0 of Fatoms")
out("")
out("waters within 5 A of F (orig): %d" % cmd.count_atoms("Wnear"))
for w in cmd.get_model("Wnear").atom:
    dz = ["%.2f" % math.dist(w.coord, z.coord) for z in Zp]
    df = ["%.2f" % math.dist(w.coord, f.coord) for f in Fm]
    out("water %s B=%.1f occ=%.2f | to Zn: %s | to F: %s" % (w.resi, w.b, w.q, dz, df))

# residues within 4 A of Zn
seen = set()
cmd.iterate("(prot and polymer.protein) within 4.0 of (prot and elem Zn)", "seen.add((resn, resi, chain))", space={"seen": seen})
out("")
out("protein residues within 4 A of Zn: " + str(sorted(seen, key=lambda x: int(x[1]))))

cmd.save(D + r"\scene_F2.pse")
open(D + r"\pymol_report.txt", "w", encoding="utf-8").write("\n".join(LOG))
out("")
out("saved scene_F2.pse and pymol_report.txt")
