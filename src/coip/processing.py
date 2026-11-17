"""Transformations applied to intensities: log2, median equalisation, summaries.

The original intensities are never modified. Every function returns a new
table, and the caller decides which one to carry forward.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def log2_intensities(intensity):
    """log2 of the raw intensities. Missing measurements stay missing.

    No offset and no imputation: a protein that was not measured has no value,
    and giving it one would invent data.
    """
    return np.log2(intensity.astype(float))


def median_equalise(log2_intensity, samples, rows=None):
    """Shift the named samples so each has the same median log2 intensity.

    The common target is the mean of those samples' medians, so the data are
    recentred rather than moved onto one arbitrary sample. Samples not named
    are returned unchanged.

    `rows` restricts the median to an anchor set of protein groups, for when
    the bulk of the table is not a fair reference , a pull-down dominated by
    contamination, for instance. The whole table is still shifted; only the
    size of the shift is decided by the anchor. Whatever the anchor is, this
    assumes the median of that set is the same in every sample, which is an
    assumption about the biology and belongs in the notebook, not here.

    Returns the shifted table and the offset added to each sample.
    """
    samples = [str(s) for s in samples]
    missing = [s for s in samples if s not in log2_intensity.columns]
    if missing:
        raise KeyError(f"not in the table: {missing}")

    anchor = log2_intensity if rows is None else log2_intensity.loc[rows]
    if anchor.empty:
        raise ValueError("the anchor set is empty")
    medians = anchor[samples].median()
    target = float(medians.mean())
    offsets = pd.Series(0.0, index=log2_intensity.columns, dtype=float)
    offsets[samples] = target - medians

    out = log2_intensity.copy()
    out[samples] = log2_intensity[samples] + offsets[samples]
    return out, offsets


def condition_summary(log2_intensity, conditions):
    """Mean log2 intensity and number of quantified samples, per condition.

    The mean is taken over the samples that have a value, so it is only
    comparable between proteins with the same number of observations , which is
    why the counts are returned alongside it.
    """
    mean = pd.DataFrame({c: log2_intensity[list(s)].mean(axis=1)
                         for c, s in conditions.items()})
    n = pd.DataFrame({c: log2_intensity[list(s)].notna().sum(axis=1)
                      for c, s in conditions.items()})
    return mean, n


def complete_cases(log2_intensity, conditions, min_per_condition=2):
    """Which protein groups have enough observations to be modelled.

    Returns a boolean mask and, for the rest, the reason they were left out.
    """
    counts = pd.DataFrame({c: log2_intensity[list(s)].notna().sum(axis=1)
                           for c, s in conditions.items()})
    eligible = (counts >= min_per_condition).all(axis=1)

    total = counts.sum(axis=1)
    reason = pd.Series("tested", index=log2_intensity.index, dtype=object)
    reason[~eligible] = (f"fewer than {min_per_condition} observations "
                         "in both conditions")
    for condition in conditions:
        short = ~eligible & (counts[condition] < min_per_condition)
        other = [c for c in conditions if c != condition]
        only_this = short & (counts[other] >= min_per_condition).all(axis=1)
        reason[only_this] = f"fewer than {min_per_condition} observations in {condition}"
    reason[total == 0] = "not quantified in either condition"

    return eligible, reason


def descriptive_comparison(log2_intensity, conditions, treatment, reference,
                           min_observations=2):
    """Mean log2 intensity per condition and the fold change between two of them.

    Descriptive only: the fold change is the difference of two means, with no
    test attached. A protein gets one when it was quantified `min_observations`
    times in each of the two compared conditions; otherwise the fold change is
    left empty and the reason recorded.
    """
    mean, observed = condition_summary(log2_intensity, conditions)
    eligible, reason = complete_cases(
        log2_intensity, {c: conditions[c] for c in (reference, treatment)},
        min_observations)

    table = pd.DataFrame(index=log2_intensity.index)
    for condition in conditions:
        table[f"n_observed_{condition}"] = observed[condition]
        table[f"mean_log2_{condition}"] = mean[condition].round(4)
    table["avg_log2"] = mean[[reference, treatment]].mean(axis=1).round(4).where(eligible)
    table["log2FC"] = (mean[treatment] - mean[reference]).round(4).where(eligible)
    table["quantified_in_both"] = eligible
    table["exclusion_reason"] = reason
    return table
