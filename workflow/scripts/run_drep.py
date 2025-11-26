"""
run_drep.py
===========

Dereplicate and refine bins using dRep.  dRep clusters genomes based
on average nucleotide identity (ANI) and chooses the best representative
from each cluster.  This script invokes ``dRep dereplicate`` on the
bins produced by DAS Tool and writes the dereplicated bins to an
output directory along with a manifest of contig membership.

Input (snakemake.input):
    bins_dir: directory containing bins to dereplicate

Output (snakemake.output):
    directory: directory where dRep writes dereplicated bins
    manifest: manifest file listing dereplicated bins and their contigs

Params (snakemake.params):
    outdir: output directory

Notes
-----
dRep will create several subdirectories in the output directory.  The
dereplicated genomes are typically found under ``<outdir>/dereplicated_genomes``.
This script collects any FASTA/FA files from that directory and builds
a manifest.
"""

import shutil
import subprocess
from pathlib import Path

# Snakemake inputs and parameters
bins_dir: Path = Path(str(snakemake.input.bins_dir))    # type: ignore
outdir: Path = Path(str(snakemake.params.outdir))        # type: ignore
threads: int = int(snakemake.threads)                    # type: ignore
manifest_path: Path = Path(str(snakemake.output[1]))     # type: ignore

outdir.mkdir(parents=True, exist_ok=True)

# Invoke dRep.  The dereplicated genomes will be placed under
# <outdir>/dereplicated_genomes
subprocess.run([
    "dRep",
    "dereplicate",
    str(outdir),
    str(bins_dir),
    "-p", str(threads)
], check=True)

# Find dereplicated bins
rep_dir = outdir / "dereplicated_genomes"
bins = sorted(list(rep_dir.glob("*.fa*")))

# Copy dereplicated bins into the top-level output directory for easier access
for bin_file in bins:
    shutil.copy2(bin_file, outdir / bin_file.name)

# Build manifest
with open(manifest_path, "w") as mf:
    mf.write("bin\tcontigs\n")
    for bin_file in bins:
        with open(bin_file) as fh:
            for line in fh:
                if line.startswith(">"):
                    contig = line[1:].strip().split()[0]
                    mf.write(f"{bin_file.name}\t{contig}\n")