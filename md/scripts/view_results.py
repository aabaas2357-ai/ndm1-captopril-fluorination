# view_results.py - PyMOL scenes comparing the crystal pose with the end of each 2 ns MD run.
# In PyMOL:  run C:\NDM1_Project\md\view_results.py      then press F1 ... F5 (scenes)
#   F1 crystal   F2 Gasteiger seed 11   F3 Gasteiger seed 22   F4 AM1-BCC seed 11   F5 all ligands overlaid
# Residue numbers in the MD files are PDB numbers minus 41; labels below show the PDB numbers.
from pymol import cmd
import os

BASE = r"C:\NDM1_Project\md"
RUNS = [("crystal", "md_5ZJ2/solvated.pdb", "green"),
        ("gast_s11", "md_5ZJ2/final_K40000.pdb", "cyan"),
        ("gast_s22", "md_5ZJ2/final_K40000_r2.pdb", "magenta"),
        ("bcc_s11", "md_5ZJ2_bcc/final_K40000.pdb", "orange")]
RES = "79+81+83+148+167+209+170+179"      # 5 Zn-binding residues + Lys211 (170) + Asn220 (179), MD numbering
LYS, ASN = "170", "179"

os.chdir(BASE)
cmd.reinitialize()
cmd.bg_color("white")
for name, path, col in RUNS:
    cmd.load(path, name)
    cmd.remove("%s and (solvent or resn NA+CL+HOH+WAT)" % name)
for name, path, col in RUNS[1:]:
    cmd.fit("%s and name CA" % name, "crystal and name CA")

cmd.hide("everything")
cmd.show("cartoon", "crystal and polymer")
cmd.color("gray70", "crystal and polymer")
cmd.set("cartoon_transparency", 0.5)
cmd.set("stick_radius", 0.15)
cmd.set("sphere_scale", 0.45)
members = {}
for name, path, col in RUNS:
    sel = "%s and (resn CAP or (polymer and resi %s)) and not hydro" % (name, RES)
    cmd.show("sticks", sel)
    cmd.color(col, sel + " and elem C")
    cmd.set("stick_radius", 0.28, "%s and resn CAP" % name)
    cmd.show("spheres", "%s and resn ZN" % name)
    zA = "first (%s and resn ZN)" % name
    zB = "last (%s and resn ZN)" % name
    lig = "%s and resn CAP" % name
    pairs = {"O2_ZnA": (lig + " and name O2", zA),
             "S_ZnA": (lig + " and name S1", zA),
             "S_ZnB": (lig + " and name S1", zB),
             "O1_LysNZ": (lig + " and name O1", "%s and resi %s and name NZ" % (name, LYS)),
             "O1_AsnN": (lig + " and name O1", "%s and resi %s and name N" % (name, ASN))}
    dn = []
    for k, (s1, s2) in pairs.items():
        cmd.distance("%s_%s" % (name, k), s1, s2)
        dn.append("%s_%s" % (name, k))
    members[name] = [name] + dn
cmd.color("slate", "resn ZN")
cmd.hide("everything", "hydro")
cmd.set("dash_gap", 0.3)
cmd.set("dash_radius", 0.06)
cmd.color("black", "*_O2_ZnA or *_S_ZnA or *_S_ZnB or *_O1_LysNZ or *_O1_AsnN")
cmd.label("crystal and name CA and resi %s" % RES, "'%s%d' % (resn, int(resi)+41)")
cmd.set("label_size", 14)
cmd.set("label_color", "black")


def make_scene(key, names, with_dist=True):
    cmd.disable("all")
    for n in names:
        for m in (members[n] if with_dist else members[n][:1]):
            cmd.enable(m)
    if "crystal" in names:
        cmd.enable("crystal")
    cmd.orient("crystal and resn CAP")
    cmd.zoom("crystal and resn CAP", 9)
    cmd.scene(key, "store")


for i, (name, path, col) in enumerate(RUNS, start=1):
    make_scene("F%d" % i, [name])
make_scene("F5", [n for n, _, _ in RUNS], with_dist=False)
cmd.scene("F1", "recall")
print("scenes ready: F1 crystal, F2 Gasteiger s11, F3 Gasteiger s22, F4 AM1-BCC s11, F5 overlay")