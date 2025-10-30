import shutil
import subprocess
from pathlib import Path

out_tsv  = Path(str(snakemake.output.tsv))        # type: ignore
dbpath   = Path(str(snakemake.params.dbpath))     # type: ignore
outdir   = Path(str(snakemake.params.outdir))     # type: ignore
bins_dir = Path(str(snakemake.params.bins_dir))   # type: ignore
threads  = int(snakemake.threads)                 # type: ignore

outdir.mkdir(parents=True, exist_ok=True)

cmd = [
    "checkm2", "predict",
    "--threads", str(threads),
    "--input", str(bins_dir),
    "--output-directory", str(outdir),
    "--database", str(dbpath)
]
subprocess.run(cmd, check=True)

# Normalize to expected TSV path
cands = [
    outdir / "quality_report.tsv",
    outdir / "checkm2_quality.tsv",
    outdir / "quality.tsv"
]
src = next((p for p in cands if p.exists()), None)
if src is None:
    tsvs = list(outdir.glob("*.tsv"))
    if not tsvs:
        raise FileNotFoundError(f"No CheckM2 TSV found in {outdir}")
    src = tsvs[0]

out_tsv.parent.mkdir(parents=True, exist_ok=True)
shutil.copy2(src, out_tsv)
