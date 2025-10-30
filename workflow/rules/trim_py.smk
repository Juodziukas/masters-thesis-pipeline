rule trim_reads_py:
    input:
        lambda wc: SAMPLE_FASTQ[wc.sample]
    output:
        f"{OUTDIR}/trim/{{sample}}.fastq.gz"
    threads: THREADS
    conda: "../envs/trim.yaml"
    params:
        tool = config.get("trim_tool", "none")
    script:
        "../scripts/run_trim.py"


