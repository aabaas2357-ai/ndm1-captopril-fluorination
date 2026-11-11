# NDM-1 computational study: fluorinated captopril analogues and zinc-aware docking

> **Status: work in progress (computational only).** Nothing here demonstrates biological activity.
> Several results are negative or inconclusive and are reported as such. Items marked `TODO` must be
> completed or removed before publishing.

## 1. Background

New Delhi metallo-β-lactamase-1 (NDM-1) is a di-zinc B1 metallo-β-lactamase that hydrolyses carbapenems.
Thiols such as captopril inhibit it by coordinating the zinc ions (PDB 5ZJ2: D-captopril, 1.10 Å).
Kondratieva et al. (*Eur J Med Chem* 2024, 266:116140) reported CF3-captopril analogues with activity
comparable to the parent compounds; the present work asks whether **4,4-difluorination of the proline
ring** (CF2 at C4) changes conformation or binding in NDM-1.

Two related questions are studied:

| | Question | Status |
|---|---|---|
| A | Can approved drugs (e-Drug3D, ~2,100 compounds) be screened against NDM-1 (PDB 3S0Z) for zinc-binding poses? | Incomplete (`TODO`: final counts) |
| B | Does 4,4-difluoro-D-captopril differ from D-captopril in strain, contacts or docking in the 5ZJ2 pocket? | Done at the level described below |

## 2. Methods (B)

* **Structures:** PDB 5ZJ2 (D-captopril in NDM-1). Ligand hydrogens and the CF2 fluorines were built by script
  on the crystal geometry, they are *not* experimentally observed.
* **Conformational strain:** GFN2-xTB with ALPB(water); bound geometry (heavy atoms frozen, L-BFGS) versus
  ETKDG-sampled free conformers.
* **Zinc cluster model:** 115 atoms (two Zn, eight residues, ligand), GFN2-xTB.
* **Docking:** AutoDock4Zn (AD4Zn) and plain AutoDock Vina 1.2.x; re-docking validation against the crystal pose.
* **Ligand library:** 4 reference compounds, CF3/CF2 variants and 5 heterocyclic variants, placed on the
  D-captopril crystal pose by scaffold alignment, scored with Vina (15 seeds).
* **Properties:** RDKit descriptors only (not true ADMET).

## 3. Results (B)

| Test | Result | Reliability |
|---|---|---|
| Conformational strain (thiolate forms) | 0.92 vs 1.08 kcal/mol (captopril vs difluoro) | Difference below method accuracy; no effect detected |
| Strain (neutral SH forms) | 16.8 vs 16.7 kcal/mol | Relaxed pose sits in a different trans basin; treated as an upper bound only |
| F contacts in the pocket | No clash after relaxation (F···Val73 2.36 Å in the rigid model, 3.57 Å after cluster relaxation); F atoms >4.5 Å from both Zn | Model-based (built fluorines) |
| Zinc cluster interaction energy | -91.1 vs -96.7 kcal/mol in the first run, but it follows the Zn–Zn geometry (differences 1.3 vs 5 kcal/mol depending on basin) | Not reliable |
| AD4Zn re-docking of captopril | Failed: best pose 3.7 Å RMSD, sulfur 4.6–5.4 Å from Zn (crystal 2.3 Å) | Protocol rejected |
| Vina re-docking of captopril | Best pose 2.09 Å RMSD (another 1.94 Å), but sulfur still 4.4–5.4 Å from Zn | Marginal; scores shape/contacts only |
| Crystal-pose control, same scoring | captopril vs difluoro differ by about 0.1 kcal/mol | No difference |

Library scoring (Vina `local_only`, mean ± SD over 15 seeds; Δ vs D-captopril, kcal/mol):

| Compound | Mean ± SD | Δ |
|---|---|---|
| D-captopril | -3.81 ± 0.20 | 0 |
| 4,4-CF2 | -3.61 ± 0.27 | +0.21 |
| CF3 on chain | -4.21 ± 0.14 | -0.40 |
| CF3 + 4,4-CF2 | -4.53 ± 0.35 | -0.72 (within one SD of CF3 alone) |

The five heterocyclic variants (tetrazole, isonipecotic-tetrazole, thioproline, indoline, thiazole) were not
better than D-captopril. A validation set of five CF3/L-captopril stereoisomers with published IC50
(Kondratieva 2024) gave Spearman ρ between 0.9 and 1.0 for `score_only` and between -0.1 and 0.8 for
`local_only` (n = 5), so ranking within close analogues is not reliable.

**Take-home:** no evidence that 4,4-difluorination changes conformation or binding; the only trend that
survives repetition is a small advantage of a chain CF3 group, consistent with the cited paper.


## 4. Limitations

* Fluorine and hydrogen positions are built, not measured; the pocket is rigid and waters are removed.
* AD4Zn is parametrised for tetrahedral Zn and did not reproduce a thiolate bridging two zinc ions.
* xTB/GFN2 is approximate for dinuclear zinc sites; cluster energies depend on the optimisation basin.
* `$fix` behaviour in the xTB runs (that frozen atoms stayed frozen) was not independently verified.
* No molecular dynamics, no DFT thiol pKa, no experimental data.

## 5. Reproduce

```bash
conda env create -f environment.yml        # docking / RDKit environment
conda activate ndm1-docking
# xTB (separate): conda create -n xtb -c conda-forge xtb=6.7.1
```

1. Download PDB 5ZJ2 and prepare it as described in `docs/` (`TODO`: add step list).
2. `scripts/captopril_docking/build_library.py` → `score_library2.py` / `hetero_run.py` → `agg.py`.
3. `scripts/captopril_xtb/` for strain and cluster calculations (paths at the top of each script must be edited).

## 6. Layout

```
scripts/screening/          FDA-library docking and Zn-distance filtering
scripts/captopril_xtb/      strain, cluster, CIP and fluorine-position checks
scripts/captopril_docking/  Vina / AD4Zn validation and library scoring
data/inputs/                ligand SDF files built from the 5ZJ2 geometry
results/                    CSV summaries
figures/                    PyMOL session script and images
```

## MD baseline (folder md/)

A separate OpenMM molecular-dynamics baseline of NDM-1 with D-captopril (PDB 5ZJ2) lives in [md/](md/README.md): set-up, restrained equilibration, 2 ns production runs with checkpoint/resume, and analysis scripts. Main finding: with a non-bonded Zn(2+) model the ligand carboxylate moves from 5.1 A to about 1.8 A from Zn in every run, so these simulations are a documented baseline, not evidence for the binding mode.

## 7. References

* Kondratieva A. et al. *Eur J Med Chem* 2024, 266, 116140.
* Brem J. et al. *Antimicrob Agents Chemother* 2016, 60, 142–150.
* Eberhardt J. et al. *J Chem Inf Model* 2021, 61, 3891 (AutoDock Vina 1.2).
* Santos-Martins D. et al. *J Chem Theory Comput* 2014, 10, 1255 (AutoDock4Zn).
* Bannwarth C. et al. *J Chem Theory Comput* 2019, 15, 1652 (GFN2-xTB).

## 8. License and acknowledgements

Code: MIT (see `LICENSE`). `TODO`: supervisor and institution, contact.
