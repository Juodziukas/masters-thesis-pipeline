import shutil
import subprocess
from pathlib import Path

# Snakemake objects
reads_in = Path(str(snakemake.input.reads))        # type: ignore
contigs  = Path(str(snakemake.output.contigs))     # type: ignore
threads  = int(snakemake.threads)                  # type: ignore
params   = snakemake.params                        # type: ignore

outdir = Path(str(params.outdir))
outdir.mkdir(parents=True, exist_ok=True)

# optional params from config/rule
asm_cov     = getattr(params, "asm_cov", None)
genome_size = getattr(params, "genome_size", None)
use_meta    = True   # our rule always runs Flye in meta mode right now

def is_gzip(p: Path) -> bool:
    try:
        with open(p, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except OSError:
        return False

# fix mislabeled .gz
reads_arg = reads_in
if reads_in.suffix == ".gz" and not is_gzip(reads_in):
    fixed = outdir / "reads.fastq"
    shutil.copy2(reads_in, fixed)
    reads_arg = fixed

# base command
cmd = [
    "flye",
    "--nano-raw", str(reads_arg),
    "--out-dir", str(outdir),
    "--threads", str(threads),
]

if use_meta:
    cmd.append("--meta")

# add genome size if we have it
if genome_size:
    cmd += ["--genome-size", str(genome_size)]

# only add asm-coverage if we are NOT in meta mode
if (not use_meta) and asm_cov:
    cmd += ["--asm-coverage", str(asm_cov)]

print("[INFO] Running Flye with command:", " ".join(cmd))

subprocess.run(cmd, check=True)

# collect output
src = outdir / "assembly.fasta"
if not src.exists():
    alt = outdir / "contigs.fasta"
    if alt.exists():
        src = alt
    else:
        raise FileNotFoundError(f"Flye output not found in {outdir}")

contigs.parent.mkdir(parents=True, exist_ok=True)
shutil.copy2(src, contigs)
