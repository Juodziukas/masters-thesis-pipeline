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
    "minimap2", "-t", str(threads), "-ax", "map-ont",
    str(contigs), str(reads)
]

# Run minimap2 with Popen to stream output
mm_proc = subprocess.Popen(mm_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# Sort to BAM (reading from minimap2's stdout)
sort_proc = subprocess.run(
    ["samtools", "sort", "-@", str(threads), "-o", str(bam)],
    stdin=mm_proc.stdout,
    check=True
)

# Close the pipe and wait for minimap2 to finish
mm_proc.stdout.close()
mm_proc.wait()

# Check if minimap2 had any errors
if mm_proc.returncode != 0:
    raise subprocess.CalledProcessError(mm_proc.returncode, mm_cmd)

# Index BAM
subprocess.run(["samtools", "index", str(bam)], check=True)

# JGI depth (MetaBAT2 utility)
subprocess.run([
    "jgi_summarize_bam_contig_depths",
    "--outputDepth", str(depth),
    str(bam)
], check=True)
