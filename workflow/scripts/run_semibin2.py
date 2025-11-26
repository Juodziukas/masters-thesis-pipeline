"""
run_semibin2.py
================

Run SemiBin2 (a deep learning‑based binning tool) on filtered contigs.
This script invokes the ``SemiBin`` command to perform single-sample
binning.  After SemiBin2 finishes, the script collects the resulting
bins and writes a manifest of contig membership.

Input (snakemake.input):
    fasta: path to filtered contigs FASTA
    depth: path to JGI depth file (optional, used for abundance information)

Output (snakemake.output):
    directory: directory where SemiBin2 writes bins
    manifest: tab‑separated file with two columns: bin filename and contig ID

Params (snakemake.params):
    outdir: output directory

The SemiBin2 CLI is still evolving; this script assumes the
``SemiBin`` executable is available and uses the ``single_easy_bin``
workflow.  Adjust the command if you use a different version of
SemiBin.
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

# Construct SemiBin2 command.  We use the single_easy_bin pipeline which
# accepts a fasta and depth file.  SemiBin writes bins into
# ``<outdir>/bins`` by default.
cmd = [
    "SemiBin",
    "single_easy_bin",
    "--input", str(fasta),
    "--depth", str(depth),
    "--output", str(outdir),
    "--threads", str(threads)
]

# Run SemiBin2
subprocess.run(cmd, check=True)

# Collect bins.  SemiBin creates a ``bins`` directory under the output
bins_dir = outdir / "bins"
bins = sorted(list(bins_dir.glob("*.fa")) + list(bins_dir.glob("*.fasta")))

# Build manifest
with open(manifest_path, "w") as mf:
    mf.write("bin\tcontigs\n")
    for bin_file in bins:
        with open(bin_file) as fh:
            for line in fh:
                if line.startswith(">"):
                    contig = line[1:].strip().split()[0]
                    mf.write(f"{bin_file.name}\t{contig}\n")