import subprocess
from pathlib import Path

fasta = Path(str(snakemake.input.fasta))          # type: ignore
depth = Path(str(snakemake.input.depth))          # type: ignore
outdir = Path(str(snakemake.params.outdir))       # type: ignore
test_mode = bool(snakemake.params.test_mode)      # type: ignore
wildcards = snakemake.wildcards                   # type: ignore
threads = int(snakemake.threads)                  # type: ignore

outdir.mkdir(parents=True, exist_ok=True)

if test_mode:
    # make 1 fake bin and a tiny manifest
    with open(outdir / "bin.1.fa", "w") as fh:
        fh.write(f">bin1_{wildcards.sample}\nATGCATGCATGC\n")
    with open(outdir / "manifest.tsv", "w") as mf:
        mf.write("bin\tcontigs\n")
        mf.write(f"bin.1.fa\tcontig1_{wildcards.sample}\n")
else:
    # real metabat2
    subprocess.run([
        "metabat2",
        "-i", str(fasta),
        "-a", str(depth),
        "-o", str(outdir / "bin"),
        "-t", str(threads)
    ], check=True)
    
    # Create manifest
    bins = sorted(outdir.glob("bin.*.fa"))
    with open(outdir / "manifest.tsv", "w") as mf:
        mf.write("bin\tcontigs\n")
        for bin_file in bins:
            # Read contigs from bin file
            with open(bin_file) as f:
                for line in f:
                    if line.startswith(">"):
                        contig = line[1:].strip().split()[0]
                        mf.write(f"{bin_file.name}\t{contig}\n")
