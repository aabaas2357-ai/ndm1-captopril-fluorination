"""
bcc_system_export.py - writes md_5ZJ2_bcc/ (system.xml with AM1-BCC ligand charges, minimized.pdb from the
crystal pose, restraints.json) so that md_run.py can run it with  --dir md_5ZJ2_bcc
Needs charge_swap_test.py in the same folder. Run in C:\\NDM1_Project\\md with the venv python.
"""
import os, shutil
exec(compile(open("charge_swap_test.py", encoding="utf-8-sig").read(), "charge_swap_test.py", "exec"))

NEW = "md_5ZJ2_bcc"
os.makedirs(NEW, exist_ok=True)
Pb = minimize(systemB)
with open(NEW + "/minimized.pdb", "w") as fh:
    app.PDBFile.writeFile(top, unit.Quantity(Pb, unit.angstrom), fh)
with open(NEW + "/system.xml", "w") as fh:
    fh.write(mm.XmlSerializer.serialize(systemB))
shutil.copy(OUT + "/restraints.json", NEW + "/restraints.json")
chk = app.PDBFile(NEW + "/minimized.pdb").getPositions(asNumpy=True).value_in_unit(unit.angstrom)
print("\nexported structure row:        ", row(chk))
print("written to", NEW, ": system.xml, minimized.pdb, restraints.json")