rule checkm2_predict_py:
    input:
        f"{OUTDIR}/binning/{{sample}}/polished_bins/manifest.tsv"
    output:
        f"{OUTDIR}/mag_qc/{{sample}}/checkm2.tsv"
    threads: 4
    params:
        db = config["dbs"]["checkm2"],
        test_mode = config.get("test_mode", False)
    run:
        import os
        import shutil
        import subprocess
        from pathlib import Path

        outdir = Path(f"{OUTDIR}/mag_qc/{wildcards.sample}")
        outdir.mkdir(parents=True, exist_ok=True)

        if params.test_mode:
            # fake output for laptop tests
            with open(outdir / "checkm2.tsv", "w") as fh:
                fh.write("bin\tcompleteness\tcontamination\n")
                fh.write("bin.1.fa\t96.5\t1.8\n")
        else:
            tmp = outdir / "checkm2_run"
            if tmp.exists():
                shutil.rmtree(tmp)
            tmp.mkdir(parents=True, exist_ok=True)

            cmd = [
                "/mnt/workspace/miniconda3/envs/checkmv2/bin/checkm2", "predict",
                "--threads", str(threads),
                "--input", f"{OUTDIR}/binning/{wildcards.sample}/polished_bins",
                "--output-directory", str(tmp),
                "--database_path", params.db,
                "-x", "fa",
            ]
            env = os.environ.copy()
            env["PATH"] = "/mnt/workspace/miniconda3/envs/checkmv2/bin:" + env["PATH"]

            subprocess.run(cmd, check=True, env=env)

            # CheckM2 usually writes a TSV in tmp, so copy it to the final output
            tsvs = list(tmp.glob("*.tsv"))
            if not tsvs:
                raise FileNotFoundError(f"No TSV produced by CheckM2 in {tmp}")
            shutil.copy2(tsvs[0], outdir / "checkm2.tsv")
