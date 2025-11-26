"""
run_das_tool.py
================

Combine bins from multiple binning algorithms using DAS Tool.  This script
invokes the ``DAS_Tool`` program on the manifests produced by MetaBAT2,
MaxBin2 and SemiBin2.  DAS Tool selects the best bins from each
algorithm based on marker genes and produces a refined set of bins.

Input (snakemake.input):
    metabat_manifest: path to the MetaBAT2 manifest
    maxbin_manifest: path to the MaxBin2 manifest
    semibin_manifest: path to the SemiBin2 manifest

Output (snakemake.output):
    directory: directory where DAS Tool writes its bins
    manifest: tab‑separated file listing the refined bins and contigs

Params (snakemake.params):
    outdir: output directory

Notes
-----
DAS Tool expects bin FASTA files rather than manifest files.  However,
it accepts a comma‑separated list of manifest files via the ``-i``
option.  This script uses that mode and instructs DAS Tool to search
for the corresponding FASTA bins in the parent directories of the
manifests.  After DAS Tool finishes, bins are stored in a directory
with the suffix ``_DASTool_bins`` in the working directory; the script
moves them into the desired output directory and builds a manifest.
"""

import shutil
import subprocess
from pathlib import Path

# Snakemake inputs and parameters
metabat_manifest: Path = Path(str(snakemake.input.metabat_manifest))     # type: ignore
maxbin_manifest: Path = Path(str(snakemake.input.maxbin_manifest))       # type: ignore
semibin_manifest: Path = Path(str(snakemake.input.semibin_manifest))     # type: ignore
outdir: Path = Path(str(snakemake.params.outdir))                        # type: ignore
threads: int = int(snakemake.threads)                                    # type: ignore
manifest_path: Path = Path(str(snakemake.output[1]))                     # type: ignore

outdir.mkdir(parents=True, exist_ok=True)

# Build DAS Tool command.  The labels correspond to each manifest and
# will be reflected in the DAS Tool summary.  Use a temporary
# basename for output so we can move the bins afterwards.
basename = outdir / "das_tool"
cmd = [
    "DAS_Tool",
    "--input", ",".join([str(metabat_manifest), str(maxbin_manifest), str(semibin_manifest)]),
    "--labels", ",".join(["metabat2", "maxbin2", "semibin2"]),
    "--output", str(basename),
    "--threads", str(threads)
]

subprocess.run(cmd, check=True)

# DAS Tool writes bins into <basename>_DASTool_bins
das_bins_dir = Path(f"{basename}_DASTool_bins")

# Move bins into the target directory
for bin_file in das_bins_dir.glob("*.fa*"):
    shutil.move(str(bin_file), str(outdir / bin_file.name))

# Build manifest from the refined bins
with open(manifest_path, "w") as mf:
    mf.write("bin\tcontigs\n")
    for bin_file in sorted(outdir.glob("*.fa*")):
        with open(bin_file) as fh:
            for line in fh:
                if line.startswith(">"):
                    contig = line[1:].strip().split()[0]
                    mf.write(f"{bin_file.name}\t{contig}\n")