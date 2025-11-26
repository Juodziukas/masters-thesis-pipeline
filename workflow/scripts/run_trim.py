import gzip
import subprocess
import shlex
from pathlib import Path

inp = Path(str(snakemake.input[0]))           # type: ignore
outp = Path(str(snakemake.output[0]))         # type: ignore
tool = str(snakemake.params.tool)             # type: ignore
cfg  = snakemake.config.get("trim", {})       # type: ignore
threads = int(snakemake.threads)

outp.parent.mkdir(parents=True, exist_ok=True)

def is_gzip(p: Path) -> bool:
    try:
        with open(p, "rb") as f:
            return f.read(2) == b"\x1f\x8b"
    except OSError:
        return False

def open_maybe_gzip(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return open(path, mode)

def write_valid_fastq(src_path: Path, dst_path: Path):
    """Small streaming validator: keeps 4-line FASTQ structure and gzips."""
    with open_maybe_gzip(src_path, "rt") as fin, gzip.open(dst_path, "wt") as fout:
        buf = []
        for line in fin:
            buf.append(line)
            if len(buf) == 4:
                if not buf[0].startswith("@"):
                    raise ValueError(f"Invalid FASTQ header in {src_path}: {buf[0]!r}")
                if not buf[2].startswith("+"):
                    raise ValueError(f"Invalid FASTQ '+' line in {src_path}: {buf[2]!r}")
                fout.writelines(buf)
                buf = []
        if buf:
            # don’t silently pass broken FASTQ to Flye
            raise ValueError(f"Incomplete FASTQ record at end of {src_path}")

if tool == "none":
    # just normalize and gzip
    write_valid_fastq(inp, outp)

elif tool == "porechop":
    # 1) ensure plain FASTQ input
    if is_gzip(inp):
        tmp_in = outp.parent / "tmp_trim_in.fastq"
        with gzip.open(inp, "rt") as fin, open(tmp_in, "w") as fout:
            for line in fin:
                fout.write(line)
    else:
        tmp_in = inp

    # 2) run porechop on files (your version wants this)
    tmp_out = outp.parent / "tmp_trim_out.fastq"
    extra = cfg.get("porechop_extra", "")
    cmd = f"porechop -i {tmp_in} -o {tmp_out} -t {threads} {extra}"
    proc = subprocess.run(shlex.split(cmd))
    if proc.returncode != 0:
        raise RuntimeError(f"porechop failed with exit code {proc.returncode}")

    # 3) validate + gzip for pipeline
    write_valid_fastq(tmp_out, outp)

elif tool == "filtlong":
    # 1) run filtlong to a temp plain fastq
    tmp_out = outp.parent / "tmp_trim_out.fastq"
    min_len = int(cfg.get("filtlong_min_length", 1000))
    keep_pct = int(cfg.get("filtlong_keep_percent", 90))
    cmd = [
        "filtlong",
        f"--min_length={min_len}",
        f"--keep_percent={keep_pct}",
        str(inp),
    ]
    with open(tmp_out, "w") as fout:
        proc = subprocess.run(cmd, stdout=fout)
    if proc.returncode != 0:
        raise RuntimeError(f"filtlong failed with exit code {proc.returncode}")

    # 2) validate + gzip
    write_valid_fastq(tmp_out, outp)

else:
    raise ValueError(f"Unknown trim tool: {tool}")
