
import os

    rule gtdbtk_classify_py:
        input:
            manifest = f"{OUTDIR}/binning/{{sample}}/polished_bins/manifest.tsv"
        output:
            tsv = f"{OUTDIR}/taxonomy/{{sample}}/gtdbtk.tsv"
        benchmark:
            f"benchmarks/binning_{{sample}}.tsv"
        threads: THREADS
        conda: "../envs/gtdbtk.yaml"
        params:
            outdir   = f"{OUTDIR}/taxonomy/{{sample}}",
            bins_dir = f"{OUTDIR}/binning/{{sample}}/polished_bins",
            data_dir = os.environ.get("GTDBTK_DATA_PATH", config["dbs"]["gtdbtk"])  # prefer env if set
        script:
            "../scripts/run_gtdbtk.py"
