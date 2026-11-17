# CK001 — PfRON4 Co-immunoprecipitation

Bait: RON4, `Q8IIK5`, part of the rhoptry neck complex. One biological preparation per
condition, each injected twice: DMSO runs 001/003, rapamycin 002/004, untagged NF54
control 005/006. The two runs within a condition are technical replicates, so the
comparisons are descriptive and carry no p-values.

Input: `input/CK001/list_of_identified_proteins.xlsx`, SHA-256
`efba217897061154baadabb97489910ced987d1349e9af1735feb2eff7e573a7`.

## How to reproduce

`notebooks/CK001_RON4.ipynb` is the whole analysis. Open it in Jupyter and run all
cells; everything in `figures/` and `tables/` comes from that run. CK001 needs only
the Python packages, no R and no SAINTexpress.

## What the experiment shows

The pull-down worked. RON4 is the first or second most intense protein group in every
tagged run, with 79 to 82 assigned spectra and 6 to 9% of the run's signal. It is also
present in both NF54 runs, at 5 and 21 spectra, because the parental line is untagged
rather than RON4-free; the gap between the two is 7.2 log2.

Ranked against the untagged line within each condition, the invasion-associated
proteins and mCherry fill most of both top tens. 1,021 protein groups could be ranked
in DMSO and 992 in RAPA; a further 538 and 434 were quantified in both injections of
their condition but in neither NF54 run, and carry no difference.

Between conditions the comparison is descriptive. 1,311 of 1,802 groups have a
complete-case fold change, the median is −0.229 log2, and RON4 itself is at −0.049:
the same amount of bait came down either way. RON2 and RON5 were recovered about
fourfold less per unit of bait under RAPA, AMA1 and RON3 about twofold.

## Limitations

1. No statistical inference is possible. One preparation per condition means the two
   runs measure instrument and injection reproducibility, not biological variability.
   No p-values, no FDR, no volcano plot.
2. The NF54 runs are the shallowest in the set, with 2,945 and 4,736 assigned spectra
   against 9,855 to 13,657 in the tagged runs, so a difference against NF54 reflects
   measurement depth as well as recovery.
3. NF54 is untagged, not RON4-free, so it is a control for the tag rather than for the
   protein.
4. Nothing is imputed; a missing measurement stays missing and is never drawn as zero.
5. Median equalisation covers the four tagged runs only. It removes a systematic offset
   between runs and does not make the measurements less variable.
6. Recovery relative to the bait is not a binding affinity and is not evidence that a
   complex came apart.

## Outputs

`figures/` holds nine figures as PDF, SVG and 400 dpi PNG:

```
01_measurement_depth      05b_top10_recovery
02_missing_measurements   06_normalisation
03_run_similarity         07_differential_recovery
04_technical_replicates   08_selected_proteins
05_bait_recovery
```

`tables/` holds the source data behind every figure, including
`qc_ranked_difference_vs_control` with its companion `qc_recovery_no_control` for the
groups that cannot be ranked, `normalisation_offsets`, `normalisation_sensitivity` and
`focus_protein_summary`.

