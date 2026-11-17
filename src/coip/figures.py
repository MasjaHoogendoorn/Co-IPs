"""The figure style from figure_style.txt, and the figures the notebooks use.

Colour roles (bait, highlighted proteins, background) are fixed here so a
colour means the same thing in every figure. Which condition gets which colour
is passed in by the notebook, because that belongs to the experiment.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

# Okabe-Ito: colour-blind safe and still distinguishable in greyscale.
OKABE_ITO = {
    "black": "#000000", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73",
    "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00",
    "purple": "#CC79A7",
}
# Colour carries meaning on two levels, and never crosses between them.
#
# Sample level , the role a run plays in the design. Set once here so that a
# condition keeps its colour in every figure of both experiments.
CONTROL_COLOUR = OKABE_ITO["green"]         # negative control IP (NF54, FR002 control)
PULLDOWN_COLOUR = OKABE_ITO["blue"]         # bait pull-down, untreated (DMSO, SUN)
TREATED_COLOUR = OKABE_ITO["vermillion"]    # bait pull-down, treated (RAPA)

# Protein level , what a point or a bar segment stands for. Orange is the
# proteins the figure is about, grey is everything it is measured against;
# that is also why the parasite class is orange and the human classes grey in
# a composition plot.
BAIT_COLOUR = OKABE_ITO["black"]
HIGHLIGHT_COLOUR = OKABE_ITO["orange"]
BACKGROUND_GREY = "#D9D9D9"
# Dark to light, for ranked subdivisions of the background.
BACKGROUND_RAMP = ("#595959", "#7F7F7F", "#A6A6A6", "#C6C6C6", "#D9D9D9")

ZERO_LINE_GREY = "#8C8C8C"
REFERENCE_RED = "#B02418"     # the twofold guides; a threshold, not a test

# One sequential map for every matrix in the project. Single hue, so it cannot
# be read as a categorical scale, and dark means "most alike" because the
# similar pairs are what the figure is asking about.
SIMILARITY_MAP = matplotlib.colors.LinearSegmentedColormap.from_list(
    "coip_similarity", ["#16324F", "#4A7BA7", "#9FC0DA", "#E8EEF4"])


def shorten_name(name, limit=26):
    """Cut a protein name at a word boundary, never mid-word.

    A name cut mid-word reads as a different protein, which on a figure is
    worse than a name that is simply shorter.
    """
    name = str(name)
    if len(name) <= limit:
        return name
    cut = name[:limit].rsplit(" ", 1)[0]
    return (cut or name[:limit]) + "\u2026"


def _readable_on(colour):
    """Black or white text, whichever the background can carry."""
    r, g, b = matplotlib.colors.to_rgb(colour)
    return "#FFFFFF" if (0.299 * r + 0.587 * g + 0.114 * b) < 0.55 else "#000000"

SINGLE_COLUMN = 85 / 25.4     # inches, for an 85 mm journal column
DOUBLE_COLUMN = 180 / 25.4    # inches, for a 180 mm journal column


def use_style():
    """Apply the project figure style. Call once, at the top of a notebook."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                            "Nimbus Sans", "DejaVu Sans"],
        "font.size": 8,
        "axes.labelsize": 8.5, "axes.titlesize": 8.5,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#000000", "axes.labelcolor": "#000000",
        "text.color": "#000000", "xtick.color": "#000000", "ytick.color": "#000000",
        "axes.linewidth": 0.75,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.width": 0.75, "ytick.major.width": 0.75,
        "axes.grid": False,
        "legend.frameon": False,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "figure.dpi": 120, "savefig.dpi": 400,
        # No tight bbox: the saved width must stay equal to the column width.
        "savefig.bbox": None, "savefig.facecolor": "white",
        "figure.constrained_layout.use": True,
        "figure.constrained_layout.h_pad": 0.03,
        "figure.constrained_layout.w_pad": 0.03,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    })


# matplotlib renamed the boxplot label argument in 3.9; Colab may ship either.
_BOXPLOT_LABELS = ("tick_labels"
                   if tuple(int(x) for x in
                            re.findall(r"\d+", matplotlib.__version__)[:2]) >= (3, 9)
                   else "labels")


def save(fig, name, directory, panel_axes=None):
    """Write the figure as PDF, SVG and 400 dpi PNG under `directory`.

    Panel letters are added only when there is genuinely more than one panel,
    in bold just outside the top-left corner of each.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    if panel_axes is not None and len(panel_axes) > 1:
        for letter, ax in zip("ABCDEFGH", np.ravel(panel_axes)):
            ax.text(-0.02, 1.04, letter, transform=ax.transAxes, fontsize=9,
                    fontweight="bold", va="bottom", ha="right")
    for extension in ("pdf", "svg", "png"):
        fig.savefig(directory / f"{name}.{extension}")
    return {extension: f"{name}.{extension}" for extension in ("pdf", "svg", "png")}


# Where a protein name may go relative to its point, as (sideways, vertical)
# offsets in points. Tried in order, away from the nearest axis edge first.
# Kept deliberately short: a name that has to travel 70 points to find room
# arrives with a leader line longer than the data it is describing.
_LABEL_ROWS = [0, 11, -11, 22, -22, 33, -33, 44, -44]
_LABEL_PAD = 5.0          # clear space to keep around a label, in points


def label_points(ax, points, labels, colours, passes=4, offsets=None,
                 avoid_y=(), avoid_points=True):
    """Name a handful of points without letting the names run into each other.

    Each label is tried in a few positions around its point and kept in the
    first one that leaves clear space around every label already placed. A thin
    leader line is drawn only when the label had to move off its own point.

    Placing a label can nudge the layout, which moves the points it was
    measured against, so the whole placement is repeated until it stops
    changing. `label_overlaps` checks the result.

    `offsets` overrides the search for individual points, as (dx, dy) in
    points or None to keep searching. `avoid_y` lists data-space y values the
    labels must not sit on, such as a zero line or a fold-change guide, and
    `avoid_points` keeps them off the marks they name. Letting the search
    handle those beats fixed offsets, which stop working the moment the
    figure changes size.
    """
    if offsets is not None and len(offsets) != len(labels):
        raise ValueError("offsets must have one entry per label")
    figure = ax.figure
    chosen, annotations = None, []
    for _ in range(passes):
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        axes_box = ax.get_window_extent(renderer)
        for annotation in annotations:
            annotation.remove()
        annotations, placed, picked = [], [], {}

        for i in sorted(range(len(labels)), key=lambda j: -points[j][1]):
            x, y = points[i]
            fixed = None if offsets is None else offsets[i]
            if fixed is not None:
                candidates = [tuple(fixed)]
            else:
                at_right = (ax.transData.transform((x, y))[0]
                            > axes_box.x0 + axes_box.width * 0.6)
                near, far = (-1, 1) if at_right else (1, -1)
                # Nearest placements first, then further out sideways: points
                # a few pixels apart need their names pushed apart laterally,
                # not just stacked.
                sides = [near * 7, far * 7, near * 30, far * 30]
                candidates = [(side, row) for row in _LABEL_ROWS
                              for side in sides]
                candidates.sort(key=lambda c: (abs(c[1]), abs(c[0])))
            forbidden = [ax.transData.transform((0, v))[1] for v in avoid_y]
            marks = ([ax.transData.transform(q) for q in points]
                     if avoid_points else [])
            for k, (dx, dy) in enumerate(candidates):
                # A label moved off its point needs a thread back to it, or
                # the reader has to guess which dot it belongs to. Sideways
                # counts as well as upwards.
                leader = abs(dy) >= 8 or abs(dx) > 12
                text = ax.annotate(
                    labels[i], (x, y), textcoords="offset points", xytext=(dx, dy),
                    fontsize=7, color=colours[i], zorder=5,
                    va="bottom" if dy > 0 and dx == 0 else
                       "top" if dy < 0 and dx == 0 else "center",
                    ha="center" if dx == 0 else "left" if dx > 0 else "right",
                    arrowprops=(dict(arrowstyle="-", lw=0.4, color=colours[i],
                                     shrinkA=0, shrinkB=3) if leader else None))
                tight = text.get_window_extent(renderer=renderer)
                box = tight.padded(_LABEL_PAD)
                # Padding buys space between labels; a reference line or a
                # mark only has to miss the glyphs themselves.
                outside = (tight.x0 < axes_box.x0 - 1
                           or tight.x1 > axes_box.x1 + 1
                           or tight.y0 < axes_box.y0 - 1
                           or tight.y1 > axes_box.y1 + 1)
                clash = (outside
                         or any(box.overlaps(other) for other in placed)
                         or any(tight.y0 - 1 <= y <= tight.y1 + 1 for y in forbidden)
                         or any(tight.x0 - 2 <= mx <= tight.x1 + 2
                                and tight.y0 - 2 <= my <= tight.y1 + 2
                                for mx, my in marks))
                if not clash or k == len(candidates) - 1:
                    placed.append(box)
                    annotations.append(text)
                    picked[i] = (dx, dy)
                    break
                text.remove()

        if picked == chosen:
            break
        chosen = picked
    return annotations


def replicate_difference(log2_intensity, replicate_pairs, colours, ax=None):
    """How far the two injections of one preparation differ, across abundance.

    One point per protein group quantified in both runs of a pair: their mean
    on the x axis, their difference on the y. A pair that sits as a band away
    from zero differs by a constant, which is what median equalisation can
    correct; a pair that fans out at low intensity differs by more where there
    is least signal to measure, which it cannot.

    This is the difference between two injections of one digest. It is not an
    estimate of biological variability, which this design does not measure.
    """
    fig = ax.figure if ax is not None else plt.figure(
        figsize=(SINGLE_COLUMN * 1.25, 2.6))
    ax = ax if ax is not None else fig.add_subplot()
    for condition, (a, b) in replicate_pairs:
        both = log2_intensity[[a, b]].dropna()
        ax.scatter(both.mean(axis=1), both[a] - both[b], s=4, alpha=0.35,
                   linewidths=0, color=colours[condition], rasterized=True,
                   label=f"{condition} ({a} \u2212 {b}), n = {len(both):,}")
    ax.axhline(0, color=OKABE_ITO["black"], lw=0.75)
    ax.set_xlabel("Mean log$_2$ intensity")
    ax.set_ylabel("log$_2$ difference between injections")
    ax.legend(loc="lower left", markerscale=2.8, handletextpad=0.3,
              borderaxespad=0.2)
    return fig, ax
def normalisation(before, after, samples, condition_of, colours,
                  replicate_pairs=None):
    """What median equalisation did, in three panels.

    A and B are the intensity distributions before and after. C shows the
    quantity the correction is actually aimed at: the difference between the
    two runs of one preparation, which should sit on zero if the two
    injections are on the same scale.

    Median equalisation adds one constant per run, so in C it slides each
    distribution along the axis without changing its width. That is the
    expected behaviour, not a shortcoming, and the panel is there to show how
    far each pair was moved towards zero.
    """
    widths = [1, 1] if replicate_pairs is None else [1, 1, 1.15]
    fig, axes = plt.subplots(1, len(widths), figsize=(DOUBLE_COLUMN, 2.5),
                             gridspec_kw=dict(width_ratios=widths))
    for ax, table in zip(axes[:2], (before, after)):
        values = [table[s].dropna().to_numpy() for s in samples]
        box = ax.boxplot(values, **{_BOXPLOT_LABELS: samples}, widths=0.6,
                         showfliers=False, patch_artist=True,
                         medianprops=dict(color="#FFFFFF", lw=1.1),
                         whiskerprops=dict(lw=0.75), capprops=dict(lw=0.75))
        for patch, sample in zip(box["boxes"], samples):
            patch.set(facecolor=colours[condition_of[sample]], alpha=1.0,
                      edgecolor=OKABE_ITO["black"], linewidth=0.6)
        ax.set_xlabel("MS run")
    axes[1].sharey(axes[0])
    axes[1].tick_params(labelleft=False)
    axes[0].set_ylabel("log$_2$ intensity")

    if replicate_pairs is not None:
        ax = axes[2]
        values, labels, faces = [], [], []
        for condition, (first, second) in replicate_pairs:
            for stage, table in (("before", before), ("after", after)):
                pair = table[[first, second]].dropna()
                values.append((pair[first] - pair[second]).to_numpy())
                labels.append(f"{condition}\n{stage}")
                faces.append((colours[condition], stage == "after"))
        positions = [i + 0.35 * (i // 2) for i in range(len(values))]
        box = ax.boxplot(values, positions=positions,
                         **{_BOXPLOT_LABELS: labels}, widths=0.62,
                         showfliers=False, patch_artist=True,
                         medianprops=dict(color=OKABE_ITO["black"], lw=1.0),
                         whiskerprops=dict(lw=0.75), capprops=dict(lw=0.75))
        for patch, (colour, filled) in zip(box["boxes"], faces):
            patch.set(facecolor=colour if filled else "white", alpha=1.0,
                      edgecolor=colour, linewidth=1.0)
        ax.axhline(0, color=ZERO_LINE_GREY, lw=0.9, zorder=0)
        ax.set_ylabel("log$_2$ difference within pair")
        ax.tick_params(axis="x", length=0)
    return fig, axes


def ma_plot(results, highlight_rows=(), bait_row=None, fold_change_line=1.0,
            x="avg_log2", y="log2FC", x_label=None, y_label=None, ax=None,
            label_offsets=None, label_rows=None):
    """Fold change against mean abundance, for every protein with both.

    The main global view of a descriptive comparison: it shows the size of
    every change and whether changes depend on how abundant a protein is.
    Rows in `highlight_rows` are drawn on top, the bait in black, and named
    unless `label_rows` restricts naming to a subset.

    The dashed guides mark a twofold change. They are a magnitude, not a test,
    and the notebook text has to say so.
    """
    data = results.dropna(subset=[x, y]).copy()
    fig = ax.figure if ax is not None else plt.figure(figsize=(DOUBLE_COLUMN, 3.4))
    ax = ax if ax is not None else fig.add_subplot()
    ax.scatter(data[x], data[y], s=5, linewidths=0, color=BACKGROUND_GREY,
               rasterized=True)
    ax.axhline(0, color=ZERO_LINE_GREY, lw=0.75)
    if fold_change_line:
        for cut in (-fold_change_line, fold_change_line):
            ax.axhline(cut, color=BACKGROUND_RAMP[2], lw=0.6, ls=(0, (4, 3)))

    low, high = data[x].min(), data[x].max()
    ax.set_xlim(low - 0.03 * (high - low), high + 0.05 * (high - low))
    # Show every point, but take the limits from the data rather than forcing
    # symmetry: one extreme value would otherwise leave half the panel empty.
    bottom, top = float(data[y].min()), float(data[y].max())
    pad = 0.07 * (top - bottom)
    ax.set_ylim(bottom - pad, top + pad)
    ax.set_xlabel(x_label or "Mean log$_2$ intensity")
    ax.set_ylabel(y_label or "log$_2$ fold change")

    rows = [r for r in highlight_rows if r in data.index]
    if rows:
        # The proteins of interest are circled rather than filled: the ring
        # marks them out by shape, so no second colour has to carry a meaning
        # here. The bait keeps its own black dot, which the ring sits around.
        ringed = [r for r in rows if r != bait_row]
        if ringed:
            ax.scatter(data.loc[ringed, x], data.loc[ringed, y], s=52, zorder=3,
                       facecolor="none", edgecolor=OKABE_ITO["black"],
                       linewidths=1.0)
        if bait_row in rows:
            ax.scatter(data.loc[[bait_row], x], data.loc[[bait_row], y], s=26,
                       zorder=4, color=BAIT_COLOUR, edgecolor="white",
                       linewidths=0.5)
        named = rows if label_rows is None else [r for r in rows
                                                 if r in set(label_rows)]
        label_points(ax, list(zip(data.loc[named, x], data.loc[named, y])),
                     list(data.loc[named, "label"]),
                     [OKABE_ITO["black"] for _ in named],
                     offsets=(None if label_offsets is None else
                              [label_offsets.get(r) for r in named]),
                     # Only the solid zero line is kept clear; the dashed
                     # guides are faint enough to read a name across.
                     avoid_y=(0.0,))
    return fig, ax


def protein_recovery(log2_intensity, rows, labels, condition_of, colours,
                     conditions, bait_row=None, ax=None,
                     condition_names=None):
    """Every individual measurement for a short list of proteins.

    One row per protein and one marker per sample, so a reader sees the spread
    behind each mean rather than a summary of it. Samples sit in one lane per
    condition, not scattered over the row: the vertical position then carries
    no meaning of its own and the two conditions read as two groups.

    A protein not quantified in a run is a missing measurement, not a low one,
    so it cannot be placed on the intensity scale. It is drawn as an open ring
    in a separate strip at the left, behind a rule, with the count of missing
    runs. Nothing is imputed and nothing is silently dropped.
    """
    samples = [s for condition in conditions for s in conditions[condition]]
    order = list(conditions)
    lane = {c: (i - (len(order) - 1) / 2) * 0.26 for i, c in enumerate(order)}
    fig = (ax.figure if ax is not None else
           plt.figure(figsize=(DOUBLE_COLUMN, 0.30 * len(rows) + 1.0)))
    ax = ax if ax is not None else fig.add_subplot()

    values = log2_intensity.loc[list(rows), samples].to_numpy(dtype=float)
    low, high = float(np.nanmin(values)), float(np.nanmax(values))
    span = high - low
    strip = low - 0.14 * span          # where the unquantified markers sit
    edge = low - 0.07 * span           # the rule between strip and data

    for y, row in enumerate(rows):
        if y:
            ax.axhline(y - 0.5, color="#EEEEEE", lw=0.6, zorder=0)
        for condition in order:
            missing = 0
            for sample in conditions[condition]:
                value = log2_intensity.at[row, sample]
                if pd.isna(value):
                    missing += 1
                    continue
                ax.plot(value, y + lane[condition], "o", ms=4.0,
                        mfc=colours[condition], mec="white", mew=0.4, zorder=3)
            if missing:
                ax.plot(strip, y + lane[condition], "o", ms=3.4, mfc="none",
                        mec=colours[condition], mew=0.9, zorder=3)
                ax.annotate(f"{missing}", (strip, y + lane[condition]),
                            textcoords="offset points", xytext=(6, 0),
                            fontsize=5.5, va="center", ha="left",
                            color=colours[condition])

    ax.axvline(edge, color=ZERO_LINE_GREY, lw=0.6, zorder=1)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{l} (bait)" if r == bait_row else l
                        for r, l in zip(rows, labels)])
    for tick, row in zip(ax.get_yticklabels(), rows):
        if row == bait_row:
            tick.set_fontweight("bold")
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("log$_2$ intensity")
    ax.set_xlim(strip - 0.04 * span, high + 0.09 * span)
    ticks = [v for v in ax.get_xticks() if low - 0.02 * span <= v <= high]
    ax.set_xticks([strip] + ticks)
    ax.set_xticklabels(["not\nquantified"] + [f"{v:g}" for v in ticks])
    for tick in ax.get_xticklabels()[:1]:
        tick.set_fontsize(6)
        tick.set_color(ZERO_LINE_GREY)

    names = condition_names or {}
    handles = [plt.Line2D([], [], marker="o", ls="", ms=4.0, mfc=colours[c],
                          mec="white", mew=0.4, label=names.get(c, c))
               for c in conditions]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.0),
              ncol=len(conditions), handletextpad=0.2, columnspacing=1.2)
    return fig, ax


def bait_relative_order(results, rows, bait_row, column="log2FC"):
    """Rows ordered by their distance from the bait, largest first."""
    distance = results.loc[list(rows), column] - results.at[bait_row, column]
    return list(distance.sort_values(ascending=False).index)


def bait_relative(results, rows, bait_row, highlight_rows=(), ax=None,
                  sort=True, highlight_colour=None):
    """Fold change of each protein beside the fold change of the bait itself.

    If the bait came down less efficiently in one condition, everything it
    carries moves with it. The distance from the bait's own line is what is
    left after that shared component.

    `sort` orders the rows by that distance. Pass False to keep the order
    given, which is what a neighbouring panel needs if the two are to be read
    across. `highlight_colour` overrides the project highlight, for a notebook
    whose palette differs.
    """
    data = results.loc[list(rows)].copy()
    bait_fc = float(results.at[bait_row, "log2FC"])
    data["relative"] = data["log2FC"] - bait_fc
    if sort:
        data = data.sort_values("relative", ascending=False)

    fig = (ax.figure if ax is not None else
           plt.figure(figsize=(SINGLE_COLUMN, 0.34 * len(data) + 1.1)))
    ax = ax if ax is not None else fig.add_subplot()
    positions = np.arange(len(data))
    for position, (row, record) in zip(positions, data.iterrows()):
        colour = (BAIT_COLOUR if row == bait_row
                  else (highlight_colour or HIGHLIGHT_COLOUR)
                  if row in set(highlight_rows) else "#808080")
        ax.hlines(position, 0, record["relative"], color=colour, lw=0.9)
        ax.plot(record["relative"], position, "o", ms=4.0, mfc=colour,
                mec="white", mew=0.4, zorder=3)
    ax.axvline(0, color=OKABE_ITO["black"], lw=0.75)
    ax.set_yticks(positions)
    ax.set_yticklabels(data["label"])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(len(data) - 0.4, -0.6)
    ax.set_xlabel("log$_2$ fold change relative to the bait")
    return fig, ax


# --- interaction confidence --------------------------------------------------

def candidate_dot_map(counts, rows, labels, sample_order, condition_of, colours,
                      scores, score_colour, bait_row=None, group_labels=None,
                      size_key=(1, 5, 10, 20), score_label="SAINT score",
                      max_area=165.0):
    """Per-run spectral counts beside the interaction score, for named proteins.

    A: one circle per protein and run, its *area* proportional to the spectral
    count, so a reader sees how many spectra there were and how evenly they
    fell across the runs. A measured zero is drawn as a short grey dash, which
    says nothing was assigned to that protein in that run; a missing cell is
    left blank, because missing and zero are not the same observation.

    B: the same proteins in the same rows, with the score the scoring program
    gave each one on a common 0 to 1 scale. The bait is drawn in the bait
    colour: it is the positive control for the method, not a discovery.

    No threshold is drawn. A cutoff line would read as a validated boundary,
    and the quantity on the axis is a model probability, not a test result.

    `rows` fixes the order top to bottom and is shared by both panels.
    """
    rows = list(rows)
    conditions = list(dict.fromkeys(condition_of[s] for s in sample_order))
    group_labels = group_labels or {c: c for c in conditions}

    # Columns sit one unit apart, with a gap where the condition changes, so
    # the two groups read as groups without a box or a bracket around them.
    position, previous, x_of = 0.0, None, {}
    for sample in sample_order:
        if previous is not None and condition_of[sample] != previous:
            position += 0.9
        x_of[sample] = position
        position += 1.0
        previous = condition_of[sample]

    shown = counts.loc[rows, sample_order].astype(float)
    largest = float(np.nanmax(shown.to_numpy())) or 1.0
    area = max_area / largest          # points^2 per spectrum: area, not radius

    fig = plt.figure(figsize=(DOUBLE_COLUMN, 0.285 * len(rows) + 1.05))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.0, 0.58],
                            height_ratios=[1.0, 0.15])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1], sharey=ax_a)
    ax_key = fig.add_subplot(grid[1, 0])
    ax_key.set_axis_off()

    for y, row in enumerate(rows):
        for sample in sample_order:
            value = shown.at[row, sample]
            if np.isnan(value):
                continue               # not measured: draw nothing at all
            if value == 0:
                ax_a.plot(x_of[sample], y, "_", ms=4.0, mew=1.0,
                          color=ZERO_LINE_GREY, zorder=2)
            else:
                ax_a.scatter(x_of[sample], y, s=value * area, zorder=3,
                             color=colours[condition_of[sample]],
                             edgecolors="white", linewidths=0.5)

    ax_a.set_xticks([x_of[s] for s in sample_order])
    ax_a.set_xticklabels(sample_order)
    ax_a.set_yticks(range(len(rows)))
    ax_a.set_yticklabels([labels[r] for r in rows])
    ax_a.set_xlim(min(x_of.values()) - 0.5, max(x_of.values()) + 0.5)
    ax_a.set_ylim(len(rows) - 0.5, -0.5)
    ax_a.tick_params(axis="both", length=0)
    for side in ("left", "bottom", "top", "right"):
        ax_a.spines[side].set_visible(False)

    blend = matplotlib.transforms.blended_transform_factory(
        ax_a.transData, ax_a.transAxes)
    for condition in conditions:
        members = [x_of[s] for s in sample_order if condition_of[s] == condition]
        ax_a.text(sum(members) / len(members), 1.02, group_labels[condition],
                  transform=blend, ha="center", va="bottom", fontsize=7.5,
                  color=colours[condition])

    for y, row in enumerate(rows):
        value = float(scores[row])
        tone = BAIT_COLOUR if row == bait_row else score_colour
        ax_b.hlines(y, 0, value, color="#BFBFBF", lw=0.9, zorder=2)
        ax_b.plot(value, y, "o", ms=5.0, mfc=tone, mec="white", mew=0.6,
                  zorder=3)
        ax_b.text(value + 0.045, y, f"{value:.2f}", va="center", ha="left",
                  fontsize=7, color=tone)

    ax_b.set_xlim(-0.02, 1.26)
    ax_b.set_xticks([0, 0.5, 1.0])
    ax_b.set_xlabel(score_label)
    ax_b.tick_params(axis="y", length=0, labelleft=False)
    ax_b.spines["left"].set_visible(False)

    # The bait is the row above the line, the candidates are below it.
    if bait_row in rows:
        for ax in (ax_a, ax_b):
            ax.axhline(rows.index(bait_row) + 0.5, color="#D0D0D0", lw=0.6,
                       zorder=1)
    for tick, row in zip(ax_a.get_yticklabels(), rows):
        if row == bait_row:
            tick.set_fontweight("bold")

    # Size key: radius goes as the square root of the count, so the spacing
    # between the circles has to as well or the big ones collide.
    radii = np.sqrt(np.asarray(size_key, dtype=float))
    radii = radii / radii.max()
    steps = np.cumsum(np.r_[0.0, radii[:-1] + radii[1:]])
    start, span = 0.30, 0.34
    at = start + steps * (span / steps[-1] if steps[-1] else 0.0)
    ax_key.set_xlim(0, 1)
    ax_key.set_ylim(0, 1)
    ax_key.text(0.0, 0.78, "Spectral count", va="center", ha="left", fontsize=7.5)
    for x, value in zip(at, size_key):
        ax_key.scatter(x, 0.78, s=value * area, color="#6E6E6E",
                       edgecolors="white", linewidths=0.5)
        ax_key.text(x, 0.06, f"{value:g}", va="bottom", ha="center", fontsize=7)
    ax_key.plot(at[-1] + 0.13, 0.78, "_", ms=4.0, mew=1.0, color=ZERO_LINE_GREY)
    ax_key.text(at[-1] + 0.17, 0.78, "0, no spectra assigned",
                va="center", ha="left", fontsize=7.5)

    # Panel letters go at the left edge of the figure, clear of the protein
    # names, so save() is called without panel_axes for this one.
    corner = matplotlib.transforms.blended_transform_factory(
        fig.transFigure, ax_a.transAxes)
    fig.text(0.004, 1.02, "A", transform=corner, fontsize=9,
             fontweight="bold", va="bottom", ha="left")
    ax_b.text(-0.09, 1.02, "B", transform=ax_b.transAxes, fontsize=9,
              fontweight="bold", va="bottom", ha="left")
    return fig, (ax_a, ax_b)


def enrichment_ranking(series, labels, emphasise=None, colours=None,
                       x_label="", ax=None, sort_by=None,
                       highlight_only=None, annotations=None,
                       annotation_header=None, series_label=True,
                       annotation_size=6.5):
    """Ranked horizontal dot plot of one value per protein, per condition.

    `highlight_only` draws the rows it marks in the series colour and the rest
    in grey, for a single-series plot where colour is free to mark the
    proteins of interest rather than a condition.

    `series` is an ordered mapping of condition name to values. With one
    condition this is a plain dot plot; with two the pair shares a row, joined
    by a thin connector so the shift between them reads as a length. Rows in
    `emphasise` get a bold label, which is how the proteins of interest are
    marked once colour is being used for the conditions.

    `annotations` is one short string per protein, in the order given, printed
    down the right of the panel: how many measurements each value rests on,
    which a ranked difference otherwise hides. `series_label` names the single
    series in the corner, and is worth switching off when the panel stands
    alone and the axis already says what the value is.
    """
    names = list(series)
    table = pd.DataFrame({name: np.asarray(series[name], dtype=float)
                          for name in names})
    reference = table[sort_by] if sort_by in table else table.mean(axis=1)
    order = np.argsort(reference.to_numpy())[::-1]
    table = table.iloc[order].reset_index(drop=True)
    labels = [labels[i] for i in order]
    # Index positionally: these may arrive as Series carrying the row labels
    # of the table they were selected from.
    flags = (np.asarray(emphasise, dtype=bool) if emphasise is not None
             else np.zeros(len(labels), dtype=bool))
    bold = [bool(flags[i]) for i in order]
    chosen = (bold if highlight_only is None
              else [bool(np.asarray(highlight_only, dtype=bool)[i]) for i in order])
    colours = colours or {name: OKABE_ITO["blue"] for name in names}

    fig = (ax.figure if ax is not None else
           plt.figure(figsize=(DOUBLE_COLUMN, 0.23 * len(labels) + 1.0)))
    ax = ax if ax is not None else fig.add_subplot()

    # Two conditions sit slightly above and below the row, so that a protein
    # whose values coincide still shows both dots instead of one hiding the
    # other.
    nudge = 0.13 if len(names) > 1 else 0.0
    rows_at = {name: (i - (len(names) - 1) / 2) * 2 * nudge
               for i, name in enumerate(names)}
    for position in range(len(table)):
        values = [table.at[position, name] for name in names]
        if len(values) > 1:
            low, high = names[int(np.argmin(values))], names[int(np.argmax(values))]
            ax.plot([min(values), max(values)],
                    [position + rows_at[low], position + rows_at[high]],
                    color="#BFBFBF", lw=1.1, zorder=2, solid_capstyle="round")
        else:
            tone = (colours[names[0]] if highlight_only is None
                    or chosen[position] else BACKGROUND_GREY)
            ax.hlines(position, 0, values[0], color=tone, lw=0.9, zorder=2)
        for name, value in zip(names, values):
            tone = colours[name]
            if highlight_only is not None and not chosen[position]:
                tone = "#9A9A9A"
            ax.plot(value, position + rows_at[name], "o", ms=5.0, mfc=tone,
                    mec="white", mew=0.6, zorder=3)

    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    for tick, is_bold in zip(ax.get_yticklabels(), bold):
        if is_bold:
            tick.set_fontweight("bold")
            if highlight_only is not None:
                tick.set_color(colours[names[0]])
    ax.set_ylim(len(labels) - 0.5, -0.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(x_label)
    # A protein can sit below the control, so the axis has to reach there.
    low = min(0.0, float(table.to_numpy().min()))
    high = float(table.to_numpy().max())
    # Room on the right for the per-row annotation, where there is one.
    right = 0.22 if annotations is not None else 0.05
    ax.set_xlim(low - 0.03 * (high - low), high + right * (high - low))

    if annotations is not None:
        for position, index in enumerate(order):
            ax.text(0.995, position, annotations[index],
                    transform=ax.get_yaxis_transform(), ha="right",
                    va="center", fontsize=annotation_size, color="#6E6E6E")
        if annotation_header:
            ax.text(0.995, -0.9, annotation_header,
                    transform=ax.get_yaxis_transform(), ha="right",
                    va="center", fontsize=annotation_size, color="#6E6E6E")
    if low < 0:
        ax.axvline(0, color=ZERO_LINE_GREY, lw=0.75, zorder=1)
    if len(names) > 1:
        handles = [plt.Line2D([], [], marker="o", ls="", ms=5.0,
                              mfc=colours[name], mec="white", mew=0.6,
                              label=name) for name in names]
        ax.legend(handles=handles, loc="lower right", handletextpad=0.2,
                  columnspacing=1.0)
    elif highlight_only is not None and series_label:
        # One condition per panel: name it on the panel, in its own colour,
        # so the panel does not depend on the caption to be read.
        ax.text(0.99, 1.01, names[0], transform=ax.transAxes, ha="right",
                va="bottom", fontsize=7.5, fontweight="bold",
                color=colours[names[0]])
    return fig, ax


def top_contributors(intensity, class_of, labels, colours, top_n=12,
                     class_order=None, italic_rows=(), class_labels=None,
                     edge_colour=None, edge_width=0.4):
    """The individual protein groups that carry most of the signal.

    A composition plot says how much of a run is contamination; this says
    which entries it actually is, which is what tells you where it came from.
    Share is of the total intensity over all runs, so the bars of one class sum
    to that class's share of the dataset.

    `class_labels` renames the classes in the legend, for a notebook whose
    composition figure already names them a particular way; `edge_colour` and
    `edge_width` match that figure's bar outline. Both default to the plain
    class name and a thin black edge, which is what CK001 uses.
    """
    edge_colour = edge_colour or OKABE_ITO["black"]
    share = intensity.sum(axis=1) / intensity.sum().sum() * 100
    top = share.nlargest(top_n).iloc[::-1]
    y = np.arange(len(top))
    # Double column: a protein named in words needs the room, and the names
    # are what makes this figure worth looking at.
    fig, ax = plt.subplots(figsize=(DOUBLE_COLUMN, 0.21 * len(top) + 0.95))
    ax.barh(y, top.to_numpy(), height=0.7, linewidth=edge_width,
            edgecolor=edge_colour,
            color=[colours[class_of[i]] for i in top.index])
    ax.set_yticks(y)
    ax.set_yticklabels([labels[i] for i in top.index])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Share of total intensity (%)")
    ax.set_ylim(-0.7, len(top) - 0.3)

    seen = [c for c in (class_order or pd.unique(class_of))
            if c in {class_of[i] for i in top.index}]
    legend = ax.legend(handles=[matplotlib.patches.Patch(
        facecolor=colours[c], edgecolor=edge_colour, linewidth=edge_width,
        label=(class_labels or {}).get(c, c)) for c in seen],
        loc="lower right", handlelength=1.2, handletextpad=0.5)
    for text, name in zip(legend.get_texts(), seen):
        if name in italic_rows:
            text.set_style("italic")
    return fig, ax


def volcano_grid(panels, rows, columns, fdr_cutoff=0.05, fold_change_line=1.0,
                 unavailable=None, labels=None, label_rows=(),
                 direction=("higher in pull-down", "higher in control"),
                 direction_colours=None, panel_letters=True):
    """One volcano per strategy, on shared axes, so the strategies can be read
    against each other.

    `panels` maps (row, column) to a limma result table, or to None for a
    combination that cannot be computed; `unavailable` then gives the reason,
    which is printed in the empty panel rather than left to the caption. A
    strategy that is impossible is a result, not a gap.
    """
    higher_colour, lower_colour = (direction_colours
                                   or (PULLDOWN_COLOUR, CONTROL_COLOUR))
    present = [f for f in panels.values() if f is not None]
    span = max(float(f["log2FC"].abs().max()) for f in present) * 1.1
    top = max(float((-np.log10(f["p_adjusted"].clip(lower=1e-300))).max())
              for f in present) * 1.2

    # A column with nothing to draw gets less room than one carrying data.
    widths = [0.55 if all(panels.get((row, column)) is None for row in rows)
              else 1.0 for column in columns]
    fig, axes = plt.subplots(len(rows), len(columns), squeeze=False,
                             figsize=(DOUBLE_COLUMN, 2.5 * len(rows)),
                             sharex=True, sharey=True, layout="constrained",
                             gridspec_kw=dict(width_ratios=widths))
    for r, row in enumerate(rows):
        for c, column in enumerate(columns):
            ax = axes[r][c]
            fit = panels.get((row, column))
            if fit is None:
                ax.text(0.5, 0.5, (unavailable or {}).get((row, column), "not computable"),
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=6.5, color="#777777", linespacing=1.5)
                for side in ("left", "bottom"):
                    ax.spines[side].set_visible(False)
                # Not set_xticks([]): these axes are shared, so clearing the
                # locator here would strip the ticks from every panel.
                ax.tick_params(labelbottom=False, labelleft=False, length=0)
            else:
                y = -np.log10(fit["p_adjusted"].clip(lower=1e-300))
                passing = fit["p_adjusted"] <= fdr_cutoff
                for mask, colour in ((~passing, BACKGROUND_GREY),
                                     (passing & (fit["log2FC"] < 0), lower_colour),
                                     (passing & (fit["log2FC"] > 0), higher_colour)):
                    ax.scatter(fit.loc[mask, "log2FC"], y[mask], s=5, linewidths=0,
                               color=colour, rasterized=True)
                ax.axhline(-np.log10(fdr_cutoff), color=REFERENCE_RED, lw=0.7,
                           ls=(0, (4, 2)), zorder=1)
                for guide in (-fold_change_line, fold_change_line):
                    ax.axvline(guide, color=ZERO_LINE_GREY, lw=0.7, ls=(0, (1, 2)),
                               zorder=1)
                ax.text(0.02, 0.97, f"{int(passing.sum())} significant\n"
                                    f"{len(fit)} tested",
                        transform=ax.transAxes, va="top", fontsize=6.5)
                # The same few proteins are named wherever they appear. A
                # protein missing from a panel was not estimable there, which
                # is not the same as an effect of zero, so nothing is drawn.
                for stagger, named in enumerate(label_rows):
                    if named not in fit.index or labels is None:
                        continue
                    x = float(fit.at[named, "log2FC"])
                    low, high = ax.get_xlim()
                    edge = x > low + 0.65 * (high - low)
                    # A tracked protein keeps the colour its own result earns
                    # and is marked by a ring, so colour never says anything
                    # except significance and direction.
                    own = fit.at[named, "p_adjusted"] <= fdr_cutoff
                    fill = (BACKGROUND_GREY if not own else
                            higher_colour if fit.at[named, "log2FC"] > 0
                            else lower_colour)
                    ax.scatter([x], [y[named]], s=20, facecolor=fill,
                               edgecolor=BAIT_COLOUR, linewidths=0.7, zorder=4)
                    ax.annotate(labels[named], (x, y[named]),
                                textcoords="offset points",
                                xytext=(-6 if edge else 6,
                                        5 + 9 * (stagger % 4)), fontsize=6,
                                ha="right" if edge else "left",
                                color=BAIT_COLOUR, zorder=5)
            if r == 0:
                ax.set_title(column, fontsize=7.5, pad=4)
            if c == 0:
                ax.set_ylabel(f"{row}\n$-$log$_{{10}}$ BH-adjusted p", fontsize=7.5)
            if r == len(rows) - 1 and fit is not None:
                ax.set_xlabel("log$_2$ fold change (pull-down / control)")
    axes[0][0].set_xlim(-span, span)
    axes[0][0].set_ylim(-0.2, top)

    # The dotted guides are a twofold magnitude. They are explained in the
    # caption rather than written on the canvas, which kept a stray label in
    # the corner of the bottom-left panel.
    if panel_letters:
        for letter, ax in zip("ABCDEFGHIJKL", np.ravel(axes)):
            ax.text(-0.01, 1.03, letter, transform=ax.transAxes, fontsize=9,
                    fontweight="bold", va="bottom", ha="right")

    handles = [plt.Line2D([], [], marker="o", ls="", ms=4, mfc=c, mec="none",
                          label=l) for l, c in
               ((direction[0], higher_colour),
                (direction[1], lower_colour),
                ("not significant", BACKGROUND_GREY))]
    fig.legend(handles=handles, loc="outside lower center", ncols=3, fontsize=7,
               handletextpad=0.3, columnspacing=1.4)
    return fig, axes


def _run_positions(samples, condition_of, gap=0.6):
    """x positions for one bar per run, with a gap between condition blocks.

    The design is then carried by position as well as by tick colour, so the
    grouping survives printing in black and white.
    """
    positions, cursor, previous = [], 0.0, None
    for sample in samples:
        condition = condition_of.get(sample) if condition_of else None
        if previous is not None and condition != previous:
            cursor += gap
        positions.append(cursor)
        cursor += 1.0
        previous = condition
    return np.asarray(positions)


def measurement_depth(depth, condition_of, colours, n_identified,
                      run_order=None, condition_names=None,
                      intensity_scale=1e7, intensity_unit=r"10$^7$"):
    """How much evidence each run produced, on the three scales that differ.

    `depth` is one row per run with `total_spectra`, `groups_with_spectra`,
    `groups_quantified` and `total_intensity`.

    The three panels are deliberately not interchangeable. **A** is sequencing
    effort, **B** is how many protein groups that effort reached, and **C** is
    how much signal came back. A run can be high on one and low on another,
    and that mismatch is the thing worth seeing.

    **B** is stacked rather than nested, because quantification is a subset of
    identification: the coloured part is what was quantified, the grey cap is
    what carried spectra without a usable intensity, and a run where every
    identified group was also quantified shows no grey at all. Those capped
    groups are what section 05 can score and section 04 cannot.
    """
    table = depth.set_index("run")
    samples = list(run_order or table.index)
    table = table.loc[samples]
    x = _run_positions(samples, condition_of)
    bar = [colours[condition_of[s]] for s in samples]

    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE_COLUMN, 2.6))

    ax = axes[0]
    ax.bar(x, table["total_spectra"], width=0.65, color=bar,
           edgecolor=OKABE_ITO["black"], linewidth=0.75)
    ax.set_ylabel("Assigned MS/MS spectra")

    ax = axes[1]
    quantified = table["groups_quantified"].to_numpy()
    spectra_only = table["groups_with_spectra"].to_numpy() - quantified
    ax.bar(x, quantified, width=0.65, color=bar,
           edgecolor=OKABE_ITO["black"], linewidth=0.6)
    ax.bar(x, spectra_only, bottom=quantified, width=0.65,
           color=BACKGROUND_RAMP[2], edgecolor=OKABE_ITO["black"], linewidth=0.6)
    ax.axhline(n_identified, color=OKABE_ITO["black"], lw=0.75, ls=(0, (3, 2)))
    ax.annotate(f"{n_identified:,} identified", (x[-1] + 0.4, n_identified),
                textcoords="offset points", xytext=(0, 3), fontsize=6.5,
                ha="right", va="bottom")
    ax.set_ylabel("Protein groups")
    ax.set_ylim(0, n_identified * 1.16)

    ax = axes[2]
    ax.bar(x, table["total_intensity"] / intensity_scale, width=0.65,
           color=bar, edgecolor=OKABE_ITO["black"], linewidth=0.6)
    ax.set_ylabel(f"Summed raw intensity ({intensity_unit})")

    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(samples)
        ax.set_xlim(x[0] - 0.75, x[-1] + 0.75)
        ax.set_xlabel("MS run")

    # One legend for the whole figure: the conditions and the grey cap of
    # panel B. The run labels are left black; the legend carries the design.
    names = condition_names or {c: c for c in colours}
    handles = [Patch(facecolor=colours[condition], label=names[condition])
               for condition in colours]
    handles.append(Patch(facecolor=BACKGROUND_RAMP[2],
                         label="spectra only, no intensity (B)"))
    fig.legend(handles=handles, loc="outside upper center",
               ncols=len(handles), handlelength=1.0, handletextpad=0.4,
               columnspacing=1.8)
    return fig, axes


def missingness_overview(detection_frequency, detection_curve, pattern_matrix,
                         runs_label="runs", condition_labels=None,
                         cmap="viridis"):
    """How complete the measurements are, and what the gaps depend on.

    `detection_frequency` is one row per possible number of runs, including
    the counts that are zero, so a bar that is empty is visible as empty.
    `pattern_matrix` is a condition-by-condition count of how completely each
    protein group was quantified.

    **B** is the panel that changes how everything else is read: detection
    depends on intensity, so the protein groups quantified in every run are a
    bright subset and not a random one. Nothing is imputed anywhere; a group
    that was not measured counts in a lower bar, never as a zero intensity.
    """
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE_COLUMN, 2.5),
                             gridspec_kw=dict(width_ratios=[1.05, 1.05, 0.9]))

    ax = axes[0]
    runs = detection_frequency["runs_detected"].to_numpy()
    groups = detection_frequency["protein_groups"].to_numpy()
    ax.bar(runs, groups, width=0.65, color=BACKGROUND_GREY,
           edgecolor=OKABE_ITO["black"], linewidth=0.75)
    for position, value in zip(runs, groups):
        ax.text(position, value, f"{value:,}" if value else "0", ha="center",
                va="bottom", fontsize=7)
    ax.set_xticks(runs)
    ax.set_xlabel("MS runs with quantified intensity")
    ax.set_ylabel("Protein groups")
    ax.set_ylim(0, max(groups) * 1.18)

    ax = axes[1]
    ax.plot(detection_curve["mean_log2_intensity"],
            detection_curve["fraction_complete"] * 100, "o-", ms=3.5, lw=1.0,
            color=OKABE_ITO["black"], mfc="white", mew=0.9)
    ax.set_xlabel("Mean observed log$_2$ intensity")
    ax.set_ylabel(f"Quantified in all {runs_label} (%)")
    ax.set_ylim(-4, 104)

    # C: a nine-cell contingency table, every value printed. Viridis, the
    # same map as the run-agreement matrix, so the notebook carries one
    # sequential scale rather than two.
    ax = axes[2]
    values = pattern_matrix.to_numpy(dtype=float)
    image = ax.imshow(values, cmap=cmap, vmin=0, vmax=values.max(),
                      aspect="auto")
    for position in range(1, values.shape[1]):      # white gutters between cells
        ax.axvline(position - 0.5, color="white", lw=1.6)
    for position in range(1, values.shape[0]):
        ax.axhline(position - 0.5, color="white", lw=1.6)
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            count = int(values[row, column])
            ax.text(column, row, f"{count:,}", ha="center", va="center",
                    fontsize=7,
                    color=_readable_on(image.cmap(image.norm(count))))
    names = condition_labels or {}
    ax.set_xticks(range(values.shape[1]), list(pattern_matrix.columns))
    ax.set_yticks(range(values.shape[0]), list(pattern_matrix.index))
    ax.set_xlabel(f"{names.get('columns', 'Pull-down')} runs quantified")
    ax.set_ylabel(f"{names.get('rows', 'Control')} runs quantified")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return fig, axes


def run_similarity_grid(agreement, samples, value="sd_log2_difference",
                        label="SD of log$_2$ intensity difference",
                        cmap="viridis", figsize=None, block_of=None,
                        diagonal_colour="#F0F0F0", annotate=True, ax=None):
    """Disagreement between every pair of MS runs, as a full square matrix.

    The companion to `run_similarity`, which draws the lower triangle in the
    project's single-hue map. This one draws all 36 cells with a blank
    diagonal, because reading a run across a full row is what shows that the
    runs fall into two blocks rather than into pairs. Both exist on purpose;
    a notebook uses one or the other, never both.

    Only meaningful if every pair was scored on the same protein groups, so
    pass an agreement table computed on a common set. `block_of` maps a run to
    its group and draws a thin white gutter between the blocks, so the design
    is visible in the layout and the run labels can stay black.
    """
    samples = list(samples)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize or (SINGLE_COLUMN * 1.25,
                                                   SINGLE_COLUMN * 1.05))
    else:
        fig = ax.figure

    grid = pd.DataFrame(np.nan, index=samples, columns=samples, dtype=float)
    for _, row in agreement.iterrows():
        if row["sample_a"] in samples and row["sample_b"] in samples:
            grid.at[row["sample_a"], row["sample_b"]] = row[value]
            grid.at[row["sample_b"], row["sample_a"]] = row[value]
    values = grid.to_numpy(dtype=float)
    if np.isnan(values[~np.eye(len(samples), dtype=bool)]).any():
        missing = [(samples[i], samples[j])
                   for i in range(len(samples)) for j in range(len(samples))
                   if i != j and np.isnan(values[i, j])]
        raise ValueError(f"no {value} for these run pairs: {missing}")

    image = ax.imshow(np.ma.masked_invalid(values),
                      cmap=plt.get_cmap(cmap).with_extremes(bad=diagonal_colour),
                      interpolation="nearest")
    if annotate:
        for i in range(len(samples)):
            for j in range(len(samples)):
                if i == j:
                    continue
                ax.text(j, i, f"{values[i, j]:.2f}", ha="center", va="center",
                        fontsize=6.5,
                        color=_readable_on(image.cmap(image.norm(values[i, j]))))

    # A white gutter wherever the group changes, on both axes.
    if block_of is not None:
        for position in range(1, len(samples)):
            if block_of[samples[position]] != block_of[samples[position - 1]]:
                for line in (ax.axvline, ax.axhline):
                    line(position - 0.5, color="white", lw=2.0)

    ax.set_xticks(range(len(samples)), samples)
    ax.set_yticks(range(len(samples)), samples)
    ax.tick_params(length=0)
    for side in ("left", "bottom", "top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("MS run")
    ax.set_ylabel("MS run")
    bar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    bar.set_label(label, fontsize=7.5)
    bar.ax.tick_params(labelsize=7)
    bar.outline.set_linewidth(0.75)
    return fig, ax


def composition_across_runs(source_counts, class_counts, source_percent,
                            class_percent, colours, labels, condition_of=None,
                            condition_names=None, figsize=None,
                            bar_width=0.72, block_gap=0.7,
                            edge_colour="#333333", edge_width=0.35,
                            count_label="Protein groups quantified",
                            percent_label="Share of raw intensity (%)"):
    """What each MS run is made of, at two levels of detail and two denominators.

    Each of the four tables is one row per run, in the order they are drawn,
    and one column per category. The left column of panels splits the run two
    ways, parasite against human; the right column splits the human side into
    its classes. The top row counts protein groups, the bottom row weighs them
    by intensity.

    Both denominators are shown because they answer different questions and
    disagree: a class can be one per cent of the protein groups and a third of
    the signal. Neither axis is left unnamed for that reason.

    Categories whose column is empty everywhere are dropped from the legend
    rather than shown as an absent category.
    """
    runs = list(source_counts.index)
    x = _run_positions(runs, condition_of, gap=block_gap)

    fig = plt.figure(figsize=figsize or (DOUBLE_COLUMN, 5.4), layout="none")
    # The legend row is a strip, not a third panel row: the keys belong just
    # under the column they describe.
    grid = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.18], left=0.09,
                            right=0.985, top=0.94, bottom=0.06, wspace=0.32,
                            hspace=0.13)
    ax_a, ax_b = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])
    ax_c, ax_d = fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])
    key_left, key_right = fig.add_subplot(grid[2, 0]), fig.add_subplot(grid[2, 1])
    key_left.axis("off")
    key_right.axis("off")

    names = condition_names or {}

    def stacked(ax, table, order):
        bottom = np.zeros(len(runs))
        for name in order:
            values = table[name].to_numpy(dtype=float)
            ax.bar(x, values, bottom=bottom, width=bar_width,
                   color=colours[name], edgecolor=edge_colour,
                   linewidth=edge_width)
            bottom += values
        ax.set_xticks(x, runs)
        ax.set_xlim(x[0] - 0.55, x[-1] + 0.55)
        ax.tick_params(axis="x", length=0, pad=4)
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        # Name each condition block above its bars, so the grouping is read
        # from the figure rather than from the caption.
        if condition_of is not None:
            for condition in dict.fromkeys(condition_of[r] for r in runs):
                at = [position for position, run in zip(x, runs)
                      if condition_of[run] == condition]
                ax.text(float(np.mean(at)), 1.045, names.get(condition, condition),
                        transform=ax.get_xaxis_transform(), ha="center",
                        va="bottom", fontsize=8)

    shown_source = [n for n in source_counts.columns if source_counts[n].sum() > 0]
    shown_class = [n for n in class_counts.columns if class_counts[n].sum() > 0]

    stacked(ax_a, source_counts, shown_source)
    stacked(ax_b, class_counts, shown_class)
    for ax in (ax_a, ax_b):
        ax.set_ylabel(count_label)
        ax.set_ylim(0, source_counts.sum(axis=1).max() * 1.08)

    stacked(ax_c, source_percent, shown_source)
    stacked(ax_d, class_percent, shown_class)
    for ax in (ax_c, ax_d):
        ax.set_ylabel(percent_label)
        ax.set_ylim(0, 100)
        ax.set_yticks(np.arange(0, 101, 20))

    def key(ax, order, ncol):
        ax.legend([Patch(facecolor=colours[name], edgecolor=edge_colour,
                         linewidth=edge_width) for name in order],
                  [labels.get(name, name) for name in order],
                  loc="upper center", ncol=ncol, frameon=False, fontsize=7,
                  handlelength=1.1, columnspacing=0.9, labelspacing=0.45,
                  borderaxespad=0)

    key(key_left, shown_source, ncol=2)
    key(key_right, shown_class, ncol=3)
    return fig, (ax_a, ax_b, ax_c, ax_d)


def bait_recovery(bait_table, colours, run_order=None, missing_note="not quantified"):
    """The bait's own evidence in every run, on the two scales it is measured on.

    `bait_table` is one row per run, as `data.bait_recovery` returns it, with a
    `spectral_counts` column added. **A** counts identified MS/MS spectra, where
    zero is a measured count: the run produced no spectrum the search assigned
    to the bait. **B** is the log2 protein intensity, which is missing rather
    than zero in a run with too little evidence to quantify. Those runs are left
    empty and named under the axis, so a gap is never read as a low value.
    """
    table = bait_table.set_index("sample")
    samples = list(run_order or table.index)
    table = table.loc[samples]
    condition_of = table["condition"].to_dict()
    x = _run_positions(samples, condition_of)
    colour = [colours[condition_of[s]] for s in samples]

    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COLUMN, 2.6),
                             layout="constrained")

    ax = axes[0]
    ax.bar(x, table["spectral_counts"], width=0.65, color=colour,
           edgecolor=OKABE_ITO["black"], linewidth=0.6)
    ax.set_ylabel("Assigned MS/MS spectra")

    ax = axes[1]
    intensity = table["bait_log2_intensity"]
    observed = intensity.notna().to_numpy()
    ax.scatter(x[observed], intensity[observed], s=26, zorder=3,
               color=[c for c, keep in zip(colour, observed) if keep],
               edgecolor=OKABE_ITO["black"], linewidth=0.5)
    span = intensity.max() - intensity.min()
    ax.set_ylim(intensity.min() - 0.5 * span, intensity.max() + 0.15 * span)
    for position, keep in zip(x, observed):
        if not keep:
            ax.annotate(missing_note, (position, ax.get_ylim()[0]),
                        textcoords="offset points", xytext=(0, 4), rotation=90,
                        fontsize=6.5, ha="center", va="bottom",
                        color=BACKGROUND_RAMP[0])
    ax.set_ylabel("Raw log$_2$ intensity")

    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(samples)
        ax.set_xlim(x[0] - 0.75, x[-1] + 0.75)
        ax.set_xlabel("MS run")
    return fig, axes
