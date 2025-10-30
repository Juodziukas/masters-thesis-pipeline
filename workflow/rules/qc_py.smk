def fq_of(sample):
    return SAMPLE_FASTQ[sample]

rule nanoplot_py:
    input:
        fastq = lambda wc: fq_of(wc.sample)
    output:
        outdir = directory(f"{OUTDIR}/qc/{{sample}}/nanoplot")
    threads: min(2, THREADS)
    conda: "../envs/qc.yaml"
    script:
        "../scripts/run_nanoplot.py"

rule multiqc_py:
    input:
        nanoplots = expand(f"{OUTDIR}/qc/{{sample}}/nanoplot", sample=SAMPLES)
    output:
        html = f"{OUTDIR}/summary/multiqc/multiqc_report.html"
    threads: 2
    conda: "../envs/qc.yaml"
    params:
        outdir = f"{OUTDIR}/summary/multiqc",
        scan_root = OUTDIR
    script:
        "../scripts/run_multiqc.py"
