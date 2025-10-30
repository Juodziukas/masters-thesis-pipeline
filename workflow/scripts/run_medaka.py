import shutil, subprocess
from pathlib import Path

draft = Path(str(snakemake.input.draft))          # type: ignore
reads = Path(str(snakemake.input.reads))          # type: ignore
polished = Path(str(snakemake.output.polished))   # type: ignore
threads = int(snakemake.threads)                  # type: ignore
model = str(snakemake.params.model)               # type: ignore
outdir = Path(str(snakemake.params.outdir))       # type: ignore

outdir.mkdir(parents=True, exist_ok=True)

cmd = ["medaka_consensus", "-i", str(reads), "-d", str(draft), "-o", str(outdir), "-t", str(threads)]
if model:
    cmd += ["-m", model]
subprocess.run(cmd, check=True)

src = outdir / "consensus.fasta"
if not src.exists():
    raise FileNotFoundError(f"Medaka did not produce {src}")
polished.parent.mkdir(parents=True, exist_ok=True)
shutil.copy2(src, polished)
