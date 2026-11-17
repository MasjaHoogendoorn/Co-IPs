"""Read and validate a ProteoFlux protein workbook, and describe its quality.

A ProteoFlux export has a 'Summary' sheet with one row per protein group and
one intensity column per MS run, plus optional per-run evidence sheets
('Spectral Counts', 'Identification Qvalue', 'Peptides (raw)').

Nothing in this module assumes a particular experiment: the caller passes the
sample identifiers it expects to find.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

def match_sample_columns(columns, sample, prefix=""):
    """Columns belonging to one sample identifier.

    A ProteoFlux column is <prefix><sample>_<experiment>_<batch>, e.g.
    'Raw_Intensities001_CK001_G26'. Searching for '001' anywhere in that string
    also hits the 'CK001' in every other column, so the identifier is required
    to sit directly after the prefix.
    """
    pattern = re.compile(rf"^{re.escape(prefix)}{re.escape(str(sample))}(?!\d)")
    return [c for c in columns if pattern.match(str(c))]


ANNOTATION_COLUMNS = [
    "UNIPROT_AC", "GENE_NAMES", "PROTEIN_DESCRIPTIONS", "FASTA_HEADERS",
    "NUM_PRECURSORS", "NUM_UNIQUE_PEPTIDES", "PROTEIN_WEIGHT",
]


def sha256(path) -> str:
    """Checksum of the input file, so a result can be traced to its data."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass
class Workbook:
    """One ProteoFlux workbook, read and checked.

    annotation  one row per protein group: accessions, gene names, description,
                peptide and precursor counts
    intensity   the same rows, one column per sample identifier, raw intensities
                with missing measurements left as NaN
    """

    path: Path
    sha256: str
    annotation: pd.DataFrame
    intensity: pd.DataFrame
    column_map: dict           # sample id -> workbook column name
    sheets: list
    checks: list = field(default_factory=list)

    @property
    def samples(self) -> list:
        return list(self.intensity.columns)

    def row_of(self, accession: str) -> int:
        """Index of the protein group containing an accession.

        Protein groups hold several accessions separated by ';', so a bait is
        matched against the members of the group, not against the whole string.
        """
        hits = [i for i, group in enumerate(self.annotation["UNIPROT_AC"])
                if accession in str(group).split(";")]
        if len(hits) != 1:
            raise KeyError(f"{accession} matches {len(hits)} protein groups, expected 1")
        return hits[0]

    def evidence(self, sheet: str) -> pd.DataFrame:
        """A per-run evidence sheet, aligned to the annotation rows.

        Used for spectral counts and identification confidence, which say how
        much peptide evidence stands behind an intensity.
        """
        if sheet not in self.sheets:
            raise KeyError(f"the workbook has no '{sheet}' sheet")
        table = pd.read_excel(self.path, sheet_name=sheet).set_index("UNIPROT_AC")
        table = table.reindex(self.annotation["UNIPROT_AC"])
        keep = {}
        for sample in self.column_map:
            hits = match_sample_columns(table.columns, sample)
            if len(hits) != 1:
                raise ValueError(
                    f"sample {sample} matches {hits} on sheet '{sheet}'")
            keep[hits[0]] = sample
        out = table[list(keep)].rename(columns=keep)
        out.index = self.annotation.index
        return out


def load_workbook(path, samples, sheet="Summary", intensity_prefix="Raw_Intensities"):
    """Read the workbook and check it before anything is calculated on it.

    `samples` is the list of sample identifiers expected in the file, e.g.
    ['001', '002', ...]. Each one must match exactly one intensity column.
    Raises on anything that would quietly corrupt the analysis.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"no workbook at {path}")

    checks = []

    def check(name, passed, detail=""):
        checks.append({"check": name, "result": "pass" if passed else "FAIL",
                       "detail": detail})
        if not passed:
            raise ValueError(f"{name} , {detail}")

    sheets = pd.ExcelFile(path).sheet_names
    check("the workbook has the expected sheet", sheet in sheets,
          f"'{sheet}' not among {sheets}")

    raw = pd.read_excel(path, sheet_name=sheet)
    check("protein identifiers are present", "UNIPROT_AC" in raw.columns,
          "no UNIPROT_AC column")

    intensity_columns = [c for c in raw.columns if str(c).startswith(intensity_prefix)]
    check("intensity columns are present", bool(intensity_columns),
          f"no column starts with '{intensity_prefix}'")

    column_map = {}
    for sample in samples:
        hits = match_sample_columns(intensity_columns, sample, intensity_prefix)
        check(f"sample {sample} maps to one column", len(hits) == 1,
              f"matched {hits}")
        column_map[str(sample)] = hits[0]
    check("no column serves two samples",
          len(set(column_map.values())) == len(column_map))
    unused = sorted(set(intensity_columns) - set(column_map.values()))
    checks.append({"check": "all intensity columns are accounted for",
                   "result": "pass" if not unused else "note",
                   "detail": f"not mapped: {unused}" if unused else ""})

    annotation = raw[[c for c in ANNOTATION_COLUMNS if c in raw.columns]].copy()
    check("protein groups are unique", not annotation["UNIPROT_AC"].duplicated().any(),
          "a protein group appears twice")

    intensity = raw[list(column_map.values())].rename(
        columns={v: k for k, v in column_map.items()}).astype(float)
    check("intensities are numeric and non-negative",
          bool((intensity.stack(future_stack=True).dropna() >= 0).all()))
    check("no intensity is exactly zero",
          not bool((intensity == 0).any().any()),
          "zeros would become -inf on the log2 scale")
    check("no protein group is missing from every sample",
          not bool(intensity.isna().all(axis=1).any()))

    return Workbook(path=path, sha256=sha256(path), annotation=annotation,
                    intensity=intensity, column_map=column_map, sheets=sheets,
                    checks=checks)


# --- quality control ---------------------------------------------------------

def detection_summary(intensity, condition_of):
    """Per sample: how many protein groups were quantified, and how many not."""
    observed = intensity.notna()
    return pd.DataFrame({
        "sample": intensity.columns,
        "condition": [condition_of[s] for s in intensity.columns],
        "quantified": observed.sum().to_numpy(),
        "missing": (~observed).sum().to_numpy(),
        "missing_fraction": (~observed).mean().round(4).to_numpy(),
    })


def detection_pattern(intensity, conditions):
    """How many protein groups fall into each pattern of presence per condition.

    'all' means quantified in every sample of that condition, 'some' in at
    least one but not all, 'none' in none. This is the honest alternative to
    imputing a value for something that was never measured.
    """
    observed = intensity.notna()
    labels = {}
    for condition, samples in conditions.items():
        n = observed[list(samples)].sum(axis=1)
        labels[condition] = np.select(
            [n == len(samples), n > 0], ["all", "some"], default="none")
    pattern = pd.DataFrame(labels, index=intensity.index)
    counts = (pattern.value_counts().rename("protein_groups")
              .reset_index().sort_values("protein_groups", ascending=False))
    return pattern, counts.reset_index(drop=True)


def pair_disagreement(log2_intensity, samples=None):
    """Agreement between every pair of samples, on the log2 scale.

    Pearson correlation is dominated by the wide dynamic range of a pull-down
    and looks excellent even when replicates disagree, so this reports the
    spread of the pairwise log2 difference as well: its SD, its median absolute
    value, and the systematic offset between the two samples.
    """
    samples = list(samples or log2_intensity.columns)
    rows = []
    for i, a in enumerate(samples):
        for b in samples[i + 1:]:
            diff = (log2_intensity[a] - log2_intensity[b]).dropna()
            if len(diff) < 3:
                continue
            both = log2_intensity[[a, b]].dropna()
            rows.append({
                "sample_a": a, "sample_b": b,
                "shared_protein_groups": len(diff),
                "sd_log2_difference": round(float(diff.std()), 3),
                "median_abs_log2_difference": round(float(diff.abs().median()), 3),
                "median_log2_offset": round(float(diff.median()), 3),
                "pearson_r": round(float(both[a].corr(both[b])), 4),
            })
    return pd.DataFrame(rows).sort_values("sd_log2_difference").reset_index(drop=True)


def bait_recovery(log2_intensity, bait_row, condition_of):
    """Did the pull-down work? Where does the bait sit in each sample.

    Rank and share of total signal say more than the raw intensity: a bait that
    is not near the top of its own pull-down is a warning about the experiment,
    not about the analysis.
    """
    rows = []
    for sample in log2_intensity.columns:
        column = log2_intensity[sample]
        observed = column.dropna()
        linear = np.power(2.0, observed)
        bait = column.loc[bait_row]
        rows.append({
            "sample": sample,
            "condition": condition_of[sample],
            "bait_log2_intensity": round(float(bait), 3) if pd.notna(bait) else np.nan,
            "rank_by_intensity": (int((observed > bait).sum()) + 1
                                  if pd.notna(bait) else np.nan),
            "of_n_quantified": len(observed),
            "percent_of_sample_signal": (round(float(2.0 ** bait / linear.sum() * 100), 2)
                                         if pd.notna(bait) else np.nan),
        })
    return pd.DataFrame(rows)


def missingness_by_detection(log2_intensity, samples):
    """How often a protein was seen, against how intense it was when seen.

    Mass spectrometry misses faint things, so missing values are usually not
    missing at random. This table shows whether that is the case here, which
    decides how the proteins excluded from a complete-case comparison should
    be read.
    """
    samples = list(samples)
    observed = log2_intensity[samples].notna().sum(axis=1)
    rows = []
    for k in range(len(samples), 0, -1):
        group = log2_intensity.loc[observed == k, samples]
        if group.empty:
            continue
        rows.append({
            "detected_in_runs": f"{k} of {len(samples)}",
            "protein_groups": int(len(group)),
            "mean_log2_of_observed": round(float(group.stack().mean()), 2),
            "max_log2_of_observed": round(float(group.stack().max()), 2)})
    table = pd.DataFrame(rows)

    mean_intensity = log2_intensity[samples].mean(axis=1)
    complete = observed == len(samples)
    deciles = []
    for name, q in (("lowest decile", 0.1), ("highest decile", 0.9)):
        mask = (mean_intensity <= mean_intensity.quantile(q) if q < 0.5
                else mean_intensity >= mean_intensity.quantile(q))
        deciles.append({"abundance": name,
                        "fraction_detected_in_all_runs":
                            f"{complete[mask].mean():.0%}"})
    return table, pd.DataFrame(deciles)


def detection_versus_abundance(log2_intensity, samples, bins=8):
    """Fraction of protein groups seen in every run, by abundance band.

    If detection were unrelated to abundance this would be a flat line. It
    never is, and how steeply it falls says how much of the missingness is
    simply the detection limit.
    """
    samples = list(samples)
    mean_intensity = log2_intensity[samples].mean(axis=1)
    complete = log2_intensity[samples].notna().sum(axis=1) == len(samples)
    edges = np.nanquantile(mean_intensity, np.linspace(0, 1, bins + 1))
    edges[-1] += 1e-9
    rows = []
    for low, high in zip(edges[:-1], edges[1:]):
        band = (mean_intensity >= low) & (mean_intensity < high)
        if not band.any():
            continue
        rows.append({"mean_log2_intensity": float(mean_intensity[band].mean()),
                     "protein_groups": int(band.sum()),
                     "fraction_complete": float(complete[band].mean())})
    return pd.DataFrame(rows)
