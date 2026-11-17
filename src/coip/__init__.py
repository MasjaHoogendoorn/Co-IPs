"""Shared code for the Co-IP mass spectrometry analyses.

Each module has one job:

    data          read a ProteoFlux workbook, validate it, compute QC numbers
    processing    log2 transformation, median equalisation, condition summaries
    differential  run R limma and proDA on a log2 intensity matrix
    saint         spectral count summaries, and the SAINTexpress run for FR002
    figures       the figure style and the publication figures

Nothing here knows about a particular experiment. Sample identifiers, baits,
conditions and normalisation choices are arguments, set in the notebook.
"""

__all__ = ["data", "processing", "differential", "saint", "figures"]
