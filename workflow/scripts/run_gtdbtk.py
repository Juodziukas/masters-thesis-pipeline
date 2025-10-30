import os
import shutil
import subprocess
from pathlib import Path

outdir   = Path(str(snakemake.params.outdir))      # type: ignore
bins_dir = Path(str(snakemake.params.bins_dir))    # type: ignore
data_dir = Path(str(snakemake.params.data_dir))    # type: ignore
tsv_out  = Path(str(snakemake.output.tsv))         # type: ignore
threads  = int(snakemake.threads)                  # type: ignore

outdir.mkdir(parents=True, exist_ok=True)

cmd = [
    "gtdbtk", "classify_wf",
    "--genome_dir", str(bins_dir),
    "--out_dir", str(outdir),
    "--data_dir", str(data_dir),
    "--cpus", str(threads),
    "--skip_ani_screen"
]
subprocess.run(cmd, check=True)

src = outdir / "classification.tsv"
if not src.exists():
    cands = list(outdir.rglob("classification.tsv"))
    if not cands:
        raise FileNotFoundError(f"Could not find classification.tsv under {outdir}")
    src = cands[0]

tsv_out.parent.mkdir(parents=True, exist_ok=True)
shutil.copy2(src, tsv_out)
