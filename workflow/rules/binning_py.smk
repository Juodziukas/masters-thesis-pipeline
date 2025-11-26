rule coverage_depth_py:
    input:
        contigs = f"{OUTDIR}/polish/{{sample}}/medaka/consensus.fasta",
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
        fasta = f"{OUTDIR}/polish/{{sample}}/medaka/consensus.fasta",
        depth = f"{OUTDIR}/binning/{{sample}}/jgi_depth.txt"
    output:
        directory(f"{OUTDIR}/binning/{{sample}}/metabat2"),
        f"{OUTDIR}/binning/{{sample}}/metabat2/manifest.tsv"
    conda: "../envs/binning.yaml"
    threads: 4
    params:
        test_mode = config.get("test_mode", False),
        outdir = lambda wildcards: f"{OUTDIR}/binning/{wildcards.sample}/metabat2"
    script:
        "../scripts/run_metabat2.py"
