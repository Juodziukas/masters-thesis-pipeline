"""
run_racon.py
================

This script performs one or more rounds of Racon polishing on an assembly
before subsequent polishing by Medaka.  It reads the draft assembly and
trimmed Oxford Nanopore reads from Snakemake inputs and writes a
consensus FASTA after the specified number of iterations.  Each round
consists of mapping the reads to the current assembly with minimap2
and polishing with racon.

Input (snakemake.input):
    draft: path to the draft assembly FASTA
    reads: path to the trimmed FASTQ/FAST5 reads

Output (snakemake.output):
    polished: path to the final racon-polished FASTA

Params (snakemake.params):
    iterations: number of racon iterations (default: 3)
    outdir: directory in which intermediate files and final consensus
            will be written

This script assumes that minimap2 and racon are available in the
current conda environment.  It uses the "map-ont" preset for Oxford
Nanopore reads.  If you wish to modify the minimap2 preset, update
the `mm2_preset` variable below.
"""

import shutil
import subprocess
from pathlib import Path

# extract Snakemake parameters with type hints (ignored by runtime)
draft: Path = Path(str(snakemake.input.draft))           # type: ignore
reads: Path = Path(str(snakemake.input.reads))           # type: ignore
outdir: Path = Path(str(snakemake.params.outdir))        # type: ignore
iterations: int = int(snakemake.params.iterations)       # type: ignore
threads: int = int(snakemake.threads)                    # type: ignore
polished_out: Path = Path(str(snakemake.output.polished))# type: ignore

# Ensure output directory exists
outdir.mkdir(parents=True, exist_ok=True)

# Set minimap2 preset for ONT reads
mm2_preset = "map-ont"

# Start polishing
current = draft
for i in range(iterations):
    round_num = i + 1
    # Paths for this iteration
    paf = outdir / f"racon_round{round_num}.paf"
    racon_out = outdir / f"racon_round{round_num}.fasta"

    # Map reads to the current draft assembly
    # minimap2 emits PAF by default with the -x preset; we use -t for threads
    subprocess.run([
        "minimap2",
        "-t", str(threads),
        "-x", mm2_preset,
        str(current),
        str(reads),
        "-o", str(paf)
    ], check=True)

    # Perform racon polishing
    subprocess.run([
        "racon",
        "-t", str(threads),
        str(reads),
        str(paf),
        str(current),
        str(racon_out)
    ], check=True)

    # Use the newly polished assembly for the next iteration
    current = racon_out

# Copy the final polished assembly to the desired output
shutil.copy2(current, polished_out)