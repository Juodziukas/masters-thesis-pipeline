rule assemble_metaflye_py:
    input:
        reads = f"{OUTDIR}/trim/{{sample}}.fastq.gz"
    output:
        contigs = f"{OUTDIR}/assembly/{{sample}}/contigs.fasta"
    benchmark:
        f"benchmarks/assembly_{{sample}}.tsv"
    threads: 4
    resources:
        mem_mb = int(config["flye"].get("memory", 4)) * 1024,
        flye_mem = int(config["flye"].get("parallel_jobs", 1))
    conda: "../envs/assemble_flye.yaml"
    params:
        outdir      = f"{OUTDIR}/assembly/{{sample}}",
        asm_cov     = config["flye"].get("asm_coverage", None),
        genome_size = config["flye"].get("genome_size", None)
    script:
        "../scripts/run_flye.py"
