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
        nanoplots = expand(f"{OUTDIR}/qc/{{sample}}/nanoplot",                      sample=SAMPLES),
        trimmed   = expand(f"{OUTDIR}/trim/{{sample}}.fastq.gz",                    sample=SAMPLES),
        assemblies= expand(f"{OUTDIR}/assembly/{{sample}}/contigs.fasta",           sample=SAMPLES),
        medaka    = expand(f"{OUTDIR}/polish/{{sample}}/medaka/consensus.fasta",    sample=SAMPLES),
        # Use manifest of dereplicated, polished bins instead of MetaBAT2-only bins
        manifests = expand(f"{OUTDIR}/binning/{{sample}}/polished_bins/manifest.tsv",    sample=SAMPLES),
        checkm2   = expand(f"{OUTDIR}/mag_qc/{{sample}}/checkm2.tsv",               sample=SAMPLES),
        gtdb      = expand(f"{OUTDIR}/taxonomy/{{sample}}/gtdbtk.tsv",              sample=SAMPLES),
        summary   = f"{OUTDIR}/summary/metrics/summary.csv"
    output:
        html = f"{OUTDIR}/summary/multiqc/multiqc_report.html"
    threads: 2
    conda: "../envs/qc.yaml"
    params:
        outdir    = f"{OUTDIR}/summary/multiqc",
        scan_root = OUTDIR
    script:
        "../scripts/run_multiqc.py"
