from pymol import cmd
import os

D = r"."
FIG = D + r"\figures"
os.makedirs(FIG, exist_ok=True)
RAY = 1   # غيّرها إلى 0 لمعاينة سريعة بدون ray tracing

cmd.reinitialize()
cmd.load(D + r"\5ZJ2_H.pdb", "prot")
cmd.load(D + r"\mco_crystal.pdb", "capto")
cmd.load(D + r"\xtal\F2captopril_xtal_Sm.sdf", "F2cap")
if os.path.exists(D + r"\5ZJ2.pdb"):
    cmd.load(D + r"\5ZJ2.pdb", "orig")
else:
    cmd.fetch("5ZJ2", name="orig", type="pdb")

# ---------- style ----------
cmd.bg_color("white")
for k, v in [("ray_trace_mode", 1), ("ray_shadow", 0), ("antialias", 2), ("depth_cue", 0),
             ("cartoon_fancy_helices", 1), ("cartoon_transparency", 0.35),
             ("stick_radius", 0.17), ("sphere_quality", 3), ("dash_gap", 0.25),
             ("dash_width", 3), ("dash_radius", 0.06), ("dash_color", "black"),
             ("label_size", 18), ("label_color", "black"), ("label_font_id", 7),
             ("ray_opaque_background", 1), ("transparency", 0.6)]:
    cmd.set(k, v)

cmd.select("zns", "prot and elem Zn")
cmd.select("coord", "prot and chain A and resi 120+122+124+189+208+250 and not elem H")
cmd.select("Fatoms", "F2cap and elem F")
cmd.select("pocket", "(byres (prot within 5 of capto)) and not elem H")
LBL = '"%s%s" % ({"HID":"His","HIE":"His","HIS":"His"}.get(resn, resn.title()), resi)'

def base():
    cmd.hide("everything")
    cmd.delete("d_*")
    cmd.show("cartoon", "prot")
    cmd.color("gray85", "prot")

def finish(name, target, buf):
    cmd.hide("everything", "elem H")
    cmd.orient(target)
    cmd.zoom(target, buf)
    cmd.scene(name, "store")
    cmd.png(FIG + "\\" + name + ".png", width=1800, height=1350, dpi=300, ray=RAY)

# ---------- F1: overview ----------
base()
cmd.show("surface", "pocket")
cmd.set("transparency", 0.6, "prot")
cmd.show("spheres", "zns"); cmd.set("sphere_scale", 0.6, "zns"); cmd.color("purpleblue", "zns")
cmd.show("sticks", "capto"); cmd.color("orange", "capto and elem C"); cmd.util.cnc("capto")
finish("F1", "capto or zns", 14)

# ---------- F2: zinc site ----------
base()
cmd.show("sticks", "coord"); cmd.color("marine", "coord and elem C"); cmd.util.cnc("coord")
cmd.show("spheres", "zns"); cmd.set("sphere_scale", 0.5, "zns"); cmd.color("purpleblue", "zns")
cmd.show("sticks", "capto"); cmd.color("orange", "capto and elem C"); cmd.util.cnc("capto")
cmd.distance("d_ZnS", "zns", "capto and elem S", 3.2)
cmd.color("yellow", "d_ZnS")
cmd.distance("d_coord", "zns", "coord and (elem N+O+S)", 2.6)
cmd.hide("labels", "d_coord")
cmd.label("coord and name CA", LBL)
finish("F2", "zns or capto", 5)

# ---------- F3: captopril vs F2 in the pocket ----------
base()
cmd.show("sticks", "capto"); cmd.color("orange", "capto and elem C"); cmd.util.cnc("capto")
cmd.show("sticks", "F2cap"); cmd.color("limegreen", "F2cap and elem C"); cmd.util.cnc("F2cap")
cmd.show("spheres", "Fatoms"); cmd.set("sphere_scale", 0.3, "Fatoms"); cmd.color("cyan", "Fatoms")
cmd.select("site", "prot and chain A and resi 67+70+73+250 and not elem H")
cmd.show("sticks", "site"); cmd.color("wheat", "site and elem C"); cmd.util.cnc("site")
cmd.show("sticks", "coord and resi 250")
cmd.show("spheres", "zns"); cmd.set("sphere_scale", 0.4, "zns"); cmd.color("purpleblue", "zns")
cmd.select("wat", "orig and resn HOH and resi 661+673")
cmd.show("spheres", "wat"); cmd.set("sphere_scale", 0.25, "wat"); cmd.color("red", "wat")
cmd.distance("d_F_Val73", "Fatoms", "prot and resn VAL and resi 73 and name CG1", 3.5)
cmd.color("red", "d_F_Val73")
cmd.distance("d_F_wat", "Fatoms", "wat", 3.0)
cmd.color("magenta", "d_F_wat")
cmd.label("site and name CA", LBL)
finish("F3", "Fatoms", 8)

cmd.scene("F1", "recall")
cmd.save(D + r"\presentation.pse")
print("done: scenes F1,F2,F3 | images in", FIG)
