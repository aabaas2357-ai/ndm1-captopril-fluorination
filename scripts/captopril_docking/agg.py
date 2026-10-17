import csv, glob, numpy as np
D = {}
for f in sorted(glob.glob("library_scores_s*.csv")):
    for r in csv.DictReader(open(f)):
        for k in ("local_only", "score_only", "S_Zn_min"):
            if r[k]:
                D.setdefault(r["id"], {}).setdefault(k, []).append(float(r[k]))
base = np.array(D["ref_Dcap"]["local_only"])
print("id | local mean+-sd | d vs D-cap | score_only mean | S-Zn mean | n")
for i, d in sorted(D.items(), key=lambda t: np.mean(t[1]["local_only"])):
    l = np.array(d["local_only"])
    print(f"{i}: {l.mean():.2f}+-{l.std():.2f} | {l.mean()-base.mean():+.2f} | "
          f"{np.mean(d['score_only']):.2f} | {np.mean(d.get('S_Zn_min', [float('nan')])):.2f} | {len(l)}")
