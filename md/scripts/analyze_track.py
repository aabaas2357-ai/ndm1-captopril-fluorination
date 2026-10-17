import csv, glob, statistics as st, os

OUT = "md_5ZJ2"
for f in sorted(glob.glob(OUT + "/track_K*.csv")):
    rows = list(csv.DictReader(open(f)))
    prod = [r for r in rows if r["stage"] == "prod"]
    print("\n==", os.path.basename(f), " production frames:", len(prod))
    for col in [c for c in rows[0] if c not in ("stage", "time_ps")]:
        v = [float(r[col]) for r in prod]
        print("  %-22s mean %.2f  sd %.2f  min %.2f  max %.2f" % (col, st.mean(v), st.pstdev(v), min(v), max(v)))

for f in sorted(glob.glob(OUT + "/prod_K*.log")):
    lines = [l.strip() for l in open(f) if l.strip()]
    head = [h.strip('"# ') for h in lines[0].split(",")]
    data = []
    for l in lines[1:]:
        try:
            data.append([float(x) for x in l.split(",")])
        except ValueError:
            pass
    print("\n==", os.path.basename(f))
    for i, name in enumerate(head):
        if any(k in name for k in ("Temperature", "Density", "Speed")):
            v = [d[i] for d in data]
            print("  %-22s mean %.2f  sd %.2f" % (name, st.mean(v), st.pstdev(v)))