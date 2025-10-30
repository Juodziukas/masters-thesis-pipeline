import os
import subprocess
from pathlib import Path

fastq = str(snakemake.input.fastq)           # type: ignore # provided by Snakemake
outdir = str(snakemake.output.outdir)        # type: ignore
threads = int(snakemake.threads)             # type: ignore

Path(outdir).mkdir(parents=True, exist_ok=True)

cmd = [
    "NanoPlot",
    "--fastq", fastq,
    "--threads", str(threads),
    "-o", outdir,
    "--loglength",
    "--plots", "hex", "dot",
]

subprocess.run(cmd, check=True)
