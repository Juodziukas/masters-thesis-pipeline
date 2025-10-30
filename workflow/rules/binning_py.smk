rule coverage_depth_py:
    input:
        contigs = f"{OUTDIR}/polish/{{sample}}/medaka.fasta",
        reads   = f"{OUTDIR}/trim/{{sample}}.fastq.gz"
    output:
        depth = f"{OUTDIR}/binning/{{sample}}/jgi_depth.txt",
        bam   = temp(f"{OUTDIR}/binning/{{sample}}/reads.sorted.bam"),
        bai   = temp(f"{OUTDIR}/binning/{{sample}}/reads.sorted.bam.bai")
    threads: THREADS
    conda: "../envs/binning.yaml"
    params:
        workdir = f"{OUTDIR}/binning/{{sample}}"
    script:
        "../scripts/run_depth.py"

rule metabat2_py:
    input:
        fasta = f"{OUTDIR}/polish/{{sample}}/medaka.fasta",
        depth = f"{OUTDIR}/binning/{{sample}}/jgi_depth.txt"
    output:
        directory(f"{OUTDIR}/binning/{{sample}}/metabat2"),
        f"{OUTDIR}/binning/{{sample}}/metabat2/manifest.tsv"
    conda: "../envs/binning.yaml"
    threads: 4
    params:
        test_mode = config.get("test_mode", False)
    run:
        import os
        from pathlib import Path

        outdir = Path(f"{OUTDIR}/binning/{wildcards.sample}/metabat2")
        outdir.mkdir(parents=True, exist_ok=True)

        if params.test_mode:
            # make 1 fake bin and a tiny manifest
            with open(outdir / "bin.1.fa", "w") as fh:
                fh.write(f">bin1_{wildcards.sample}\nATGCATGCATGC\n")
            with open(outdir / "manifest.tsv", "w") as mf:
                mf.write("bin\tcontigs\n")
                mf.write(f"bin.1.fa\tcontig1_{wildcards.sample}\n")
        else:
            # real metabat2
            snakemake.script("../scripts/run_metabat2.py")
