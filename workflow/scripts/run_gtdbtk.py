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

# Check if any bin files exist
bin_files = list(bins_dir.glob("bin.*.fa"))
if not bin_files:
    raise FileNotFoundError(f"No bin files found in {bins_dir}. MetaBAT2 may not have produced any bins.")

# Set GTDBTK_DATA_PATH as environment variable
env = os.environ.copy()
env["GTDBTK_DATA_PATH"] = str(data_dir)

cmd = [
    "gtdbtk", "classify_wf",
    "--genome_dir", str(bins_dir),
    "--out_dir", str(outdir),
    "--extension", "fa",  # Explicitly specify .fa extension
    "--cpus", str(threads),
    "--skip_ani_screen"
]
subprocess.run(cmd, check=True, env=env)

# GTDB-Tk v2 no longer outputs classification.tsv
# Instead it writes gtdbtk.*.summary.tsv files.
summary_files = list(outdir.glob("gtdbtk.*.summary.tsv"))

if not summary_files:
    summary_files = list(outdir.rglob("gtdbtk.*.summary.tsv"))

if not summary_files:
    raise FileNotFoundError(f"Could not find any gtdbtk.*.summary.tsv under {outdir}")

src = summary_files[0]

tsv_out.parent.mkdir(parents=True, exist_ok=True)
shutil.copy2(src, tsv_out)
