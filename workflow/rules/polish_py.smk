rule medaka_polish_py:
    input:
        draft = f"{OUTDIR}/assembly/{{sample}}/contigs.fasta",
        reads = f"{OUTDIR}/trim/{{sample}}.fastq.gz"
    output:
        polished = f"{OUTDIR}/polish/{{sample}}/medaka/consensus.fasta"
    benchmark:
        f"benchmarks/polish_medaka_{{sample}}.tsv"
    threads: THREADS
    conda: "../envs/polish.yaml"
    params:
        model = config.get("polish", {}).get("medaka_model", ""),
        outdir = f"{OUTDIR}/polish/{{sample}}/medaka"
    script:
        "../scripts/run_medaka.py"
