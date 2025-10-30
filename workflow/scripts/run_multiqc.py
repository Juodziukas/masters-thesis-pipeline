import os
import shutil
import subprocess
from pathlib import Path

outdir = Path(str(snakemake.params.outdir))      # type: ignore
scan_root = Path(str(snakemake.params.scan_root))# type: ignore
target_html = Path(str(snakemake.output.html))   # type: ignore

outdir.mkdir(parents=True, exist_ok=True)

cmd = [
    "multiqc",
    str(scan_root),
    "-o", str(outdir),
]
subprocess.run(cmd, check=True)


default_html = outdir / "multiqc_report.html"
if default_html.exists():
    if target_html != default_html:
        shutil.copy2(default_html, target_html)
else:
    htmls = list(outdir.glob("*.html"))
    if not htmls:
        raise FileNotFoundError(f"MultiQC did not produce an HTML report in {outdir}")
    shutil.copy2(htmls[0], target_html)
