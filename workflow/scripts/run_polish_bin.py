"""
run_polish_bin.py
=================

Perform per‑bin polishing on dereplicated bins.  Each bin is polished
with a single round of Racon followed by Medaka.  The rationale is
that after global polishing and dereplication, individual bins may
still contain residual errors or chimeras.  Re‑mapping the reads to
each bin and polishing improves consensus quality without altering the
assembly structure significantly.

Input (snakemake.input):
    manifest: path to the manifest file produced by the dRep stage
    reads: path to the trimmed reads FASTQ

Output (snakemake.output):
    directory: directory containing polished bins (FASTA files)
    manifest: manifest listing polished bin files and contig IDs

Params (snakemake.params):
    outdir: directory where polished bins will be written

The script assumes that minimap2, racon and medaka_consensus are
available in the active conda environment.  It runs one iteration of
Racon followed by a single Medaka consensus call.  If you wish to
adjust the number of Racon rounds or specify a Medaka model, modify
the commands below or extend the rule parameters accordingly.
"""

import csv
import shutil
import subprocess
from pathlib import Path

# Snakemake inputs and parameters
manifest_path: Path = Path(str(snakemake.input.manifest))    # type: ignore
reads: Path = Path(str(snakemake.input.reads))               # type: ignore
outdir: Path = Path(str(snakemake.params.outdir))            # type: ignore
threads: int = int(snakemake.threads)                        # type: ignore
polished_manifest_path: Path = Path(str(snakemake.output[1]))# type: ignore

outdir.mkdir(parents=True, exist_ok=True)

# Determine the directory containing the dereplicated bins.  The
# manifest lists bin filenames relative to this directory.
bins_dir = manifest_path.parent

# Read unique bin filenames from the manifest
bin_names = []
with open(manifest_path) as mf:
    reader = csv.DictReader(mf, delimiter="\t")
    for row in reader:
        name = row["bin"]
        if name not in bin_names:
            bin_names.append(name)

# For each bin, perform Racon followed by Medaka
for bin_name in bin_names:
    bin_path = bins_dir / bin_name
    if not bin_path.exists():
        # Skip missing bins
        continue
    # Temporary files for this bin
    base = bin_name.split('.')[0]
    paf = outdir / f"{base}.paf"
    racon_out = outdir / f"{base}_racon.fasta"
    medaka_dir = outdir / f"{base}_medaka"

    # 1. Map reads to the bin using minimap2
    subprocess.run([
        "minimap2",
        "-t", str(threads),
        "-x", "map-ont",
        str(bin_path),
        str(reads),
        "-o", str(paf)
    ], check=True)

    # 2. Racon polishing
    subprocess.run([
        "racon",
        "-t", str(threads),
        str(reads),
        str(paf),
        str(bin_path),
        str(racon_out)
    ], check=True)

    # 3. Medaka consensus; do not specify model to let medaka choose default
    subprocess.run([
        "medaka_consensus",
        "-i", str(reads),
        "-d", str(racon_out),
        "-o", str(medaka_dir),
        "-t", str(threads)
    ], check=True)

    # Copy final consensus to output directory with the original bin name
    final_consensus = medaka_dir / "consensus.fasta"
    dest = outdir / bin_name
    shutil.copy2(final_consensus, dest)

# Build manifest for polished bins
with open(polished_manifest_path, "w") as mf:
    mf.write("bin\tcontigs\n")
    for bin_file in sorted(outdir.glob("*.fa*")):
        with open(bin_file) as fh:
            for line in fh:
                if line.startswith(">"):
                    contig = line[1:].strip().split()[0]
                    mf.write(f"{bin_file.name}\t{contig}\n")