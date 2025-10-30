import subprocess
from pathlib import Path

contigs = Path(str(snakemake.input.contigs))      # type: ignore
reads   = Path(str(snakemake.input.reads))        # type: ignore
depth   = Path(str(snakemake.output.depth))       # type: ignore
bam     = Path(str(snakemake.output.bam))         # type: ignore
bai     = Path(str(snakemake.output.bai))         # type: ignore
threads = int(snakemake.threads)                  # type: ignore
workdir = Path(str(snakemake.params.workdir))     # type: ignore

workdir.mkdir(parents=True, exist_ok=True)

# minimap2 → BAM (sorted)
mm_cmd = [
    "minimap2", "-t", str(threads), "-x", "map-ont",
    str(contigs), str(reads)
]

# Run minimap2 and capture SAM on stdout
mm_proc = subprocess.run(mm_cmd, check=True, capture_output=True)
# Sort to BAM
subprocess.run(
    ["samtools", "sort", "-@", str(threads), "-o", str(bam)],
    input=mm_proc.stdout,
    check=True
)
# Index BAM
subprocess.run(["samtools", "index", str(bam)], check=True)

# JGI depth (MetaBAT2 utility)
subprocess.run([
    "jgi_summarize_bam_contig_depths",
    "--outputDepth", str(depth),
    str(bam)
], check=True)
