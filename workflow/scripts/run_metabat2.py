import csv
import subprocess
from pathlib import Path

contigs  = Path(str(snakemake.input.contigs))     # type: ignore
depth    = Path(str(snakemake.input.depth))       # type: ignore
bins_dir = Path(str(snakemake.output.bins_dir))   # type: ignore
manifest = Path(str(snakemake.output.manifest))   # type: ignore
prefix   = str(snakemake.params.prefix)           # type: ignore
threads  = int(snakemake.threads)                 # type: ignore

bins_dir.mkdir(parents=True, exist_ok=True)

# Run metabat2
mb_cmd = [
    "metabat2",
    "-i", str(contigs),
    "-a", str(depth),
    "-o", prefix,
    "-t", str(threads)
]
subprocess.run(mb_cmd, check=True)

# Collect produced bins (bin.*.fa) and write manifest
bin_files = sorted(bins_dir.glob("bin.*.fa"))
with manifest.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["bin_id", "path"])
    for bf in bin_files:
        w.writerow([bf.stem, str(bf.resolve())])
