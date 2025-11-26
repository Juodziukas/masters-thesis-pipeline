"""
run_filter_contigs.py
====================

Filter contigs prior to binning.  This script reads a FASTA of
assembled contigs and a depth file produced by
``jgi_summarize_bam_contig_depths`` and writes a filtered FASTA
containing only those contigs that meet the following criteria:

  * Minimum contig length (default 2,500 bp)
  * Minimum mean coverage (default 5×)
  * GC content within a specified range (default 0.2–0.8)

The depth file is assumed to have whitespace‑separated columns with
the contig name in the first column and the mean coverage in the
third column.  If coverage is missing for a contig, zero is assumed.

Input (snakemake.input):
    fasta: path to contig FASTA
    depth: path to JGI depth file

Output (snakemake.output):
    filtered: path to filtered contig FASTA

Params (snakemake.params):
    min_length: minimum contig length
    min_coverage: minimum mean coverage
    min_gc: minimum GC fraction (0–1)
    max_gc: maximum GC fraction (0–1)
    outdir: directory to place the filtered FASTA

Dependencies: biopython and pandas must be installed in the conda environment.
"""

from pathlib import Path
import pandas as pd
from Bio import SeqIO

# Snakemake parameters
fasta_path: Path = Path(str(snakemake.input.fasta))          # type: ignore
depth_path: Path = Path(str(snakemake.input.depth))          # type: ignore
filtered_path: Path = Path(str(snakemake.output.filtered))    # type: ignore
min_length: int = int(snakemake.params.min_length)            # type: ignore
min_coverage: float = float(snakemake.params.min_coverage)    # type: ignore
min_gc: float = float(snakemake.params.min_gc)                # type: ignore
max_gc: float = float(snakemake.params.max_gc)                # type: ignore

# Create output directory if necessary
filtered_path.parent.mkdir(parents=True, exist_ok=True)

# Read depth file into a dictionary mapping contig IDs to coverage
coverage_dict = {}
if depth_path.exists():
    with open(depth_path) as df:
        header = df.readline()  # skip header
        for line in df:
            parts = line.strip().split()
            if not parts:
                continue
            contig = parts[0]
            # If there are at least 3 columns, take the third as mean coverage
            try:
                cov = float(parts[2]) if len(parts) > 2 else float(parts[1])
            except (IndexError, ValueError):
                cov = 0.0
            coverage_dict[contig] = cov

# Filter contigs
with open(filtered_path, "w") as out_fh:
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        length = len(record.seq)
        cov = coverage_dict.get(record.id, 0.0)
        # compute GC fraction
        seq = record.seq.upper()
        gc_count = seq.count("G") + seq.count("C")
        gc = (gc_count / length) if length else 0.0
        if length >= min_length and cov >= min_coverage and min_gc <= gc <= max_gc:
            SeqIO.write(record, out_fh, "fasta")