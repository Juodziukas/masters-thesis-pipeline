"""
run_maxbin2.py
==============

Run MaxBin2 on filtered contigs to generate an alternative set of metagenomic bins.
This script expects the `run_MaxBin.pl` wrapper to be available in the current
conda environment.  It uses the coverage file (depth) to weight contigs during
binning.  After MaxBin2 finishes, the script writes a manifest file listing
each bin and its constituent contigs.

Input (snakemake.input):
    fasta: path to filtered contigs FASTA
    depth: path to JGI depth file with contig abundances

Output (snakemake.output):
    directory: directory to which bins will be written
    manifest: tab‑separated file with two columns: bin filename and contig ID

Params (snakemake.params):
    outdir: path to output directory (used by MaxBin2)

Notes
-----
MaxBin2 produces bin FASTA files with names like ``bin.001.fasta`` in the
specified output directory.  This script collects all files with the suffix
``.fasta`` or ``.fa`` from that directory to build the manifest.
"""

import subprocess
from pathlib import Path

# Snakemake inputs and parameters
fasta: Path = Path(str(snakemake.input.fasta))                 # type: ignore
depth: Path = Path(str(snakemake.input.depth))                 # type: ignore
outdir: Path = Path(str(snakemake.params.outdir))              # type: ignore
threads: int = int(snakemake.threads)                          # type: ignore
manifest_path: Path = Path(str(snakemake.output[1]))           # type: ignore

# Create output directory
outdir.mkdir(parents=True, exist_ok=True)

# Run MaxBin2 via its Perl wrapper
subprocess.run([
    "run_MaxBin.pl",
    "-contig", str(fasta),
    "-abund", str(depth),
    "-out", str(outdir / "bin"),
    "-thread", str(threads)
], check=True)

# Collect bin files; MaxBin2 numbers bins starting at 001
bins = sorted(list(outdir.glob("bin.*.fasta")) + list(outdir.glob("bin.*.fa")))

# Build manifest
with open(manifest_path, "w") as mf:
    mf.write("bin\tcontigs\n")
    for bin_file in bins:
        with open(bin_file) as fh:
            for line in fh:
                if line.startswith(">"):
                    contig = line[1:].strip().split()[0]
                    mf.write(f"{bin_file.name}\t{contig}\n")