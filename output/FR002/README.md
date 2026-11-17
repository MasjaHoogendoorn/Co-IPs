# FR002 — PfSUN1 Co-immunoprecipitation

Bait: PfSUN1, `PF3D7_1215100` (UniProt `Q8I5Q5`), an inner-nuclear-membrane SUN-domain
protein. Three independent biological preparations per condition: control runs
001/003/005, PfSUN1 pull-downs 002/004/006.

Input: `input/FR002/list_of_identified_proteins_FR002_proteoflux_table.xlsx`, SHA-256
`fa4a2ab4c0273ee32e0e765f5eafa897af7500465689ba24ad8b4009a9a07e08`.

## How to reproduce

`notebooks/FR002_SUN.ipynb` is the whole analysis. Open it in Jupyter and run all
cells; it needs nothing from any other notebook. It requires R with limma and proDA,
and downloads SAINTexpress-spc into `.tools/` on first use. Everything in
`figures/`, `tables/` and `saint/` comes from that run.

## What the experiment shows

The pull-down worked. PfSUN1 is quantified in all three pull-downs (14, 15 and 14
spectra, about 1% of each run's signal) and in one control at a single spectrum.

The dataset is shallow and contaminated. 567 protein groups, of which 223 are parasite.
Three quarters of the intensity is human, including an epidermal and sweat signature
that varies between runs.

Differential recovery is exploratory. Thirteen parasite groups are enriched at BH ≤ 0.05
(limma, 3 against 3, residual df 4), but the same contrast fitted under six analyses
gives counts from 0 to 29, and none of the thirteen is significant under proDA. The
strongest, `PF3D7_1413200`, reaches BH 0.07 there and SAINTexpress scores it 0.00.

Two protein groups reach SAINT BFDR ≤ 0.01: the bait, and dipeptidyl aminopeptidase 2
(`PF3D7_1247800`, ten spectra across the pull-downs, none in any control). A DnaJ
protein follows at BFDR 0.11.

Seven of the ten highest-scoring preys overall are human, scoring 0.23 to 0.33, and all
of them have their pull-down spectra in run 004 alone. A score in that range reflects
very little replicate evidence in this dataset.

## Limitations

1. The 89 testable parasite groups are a bright subset, and the bait is not among them
   because it is absent from two controls and nothing is imputed.
2. Outside the bait, candidates rest on 2 to 10 spectra. SAINT's BFDR at this depth is a
   ranking, not a calibrated error rate.
3. Nothing here is evidence of direct binding.

## Outputs

`figures/` holds eleven figures as PDF, SVG and 400 dpi PNG:

```
01_measurement_depth      06_bait_recovery
02_missing_measurements   06b_top10_recovery
03_run_similarity         07_normalisation
04_composition_across_runs  08_differential_recovery
05_top_contributors       10_candidate_evidence_with_human
                          11_candidate_evidence_parasite
```

`tables/` holds the source data behind every figure, including `differential_recovery`,
`candidate_proteins`, `candidate_evidence`, `saint_scores` and the six `diff_*` fits
that figure 08 is drawn from.

`saint/` holds the three files SAINTexpress was given and the `list.txt` it returned:
`bait.dat` (IP, bait, T or C), `prey.dat` (prey id, length, gene name) and `inter.dat`
(IP, bait, prey, spectral count), tab-separated without headers, so the run can be
repeated by hand.
