rule summary_aggregate_py:
    input:
        contigs   = expand(f"{OUTDIR}/polish/{{sample}}/medaka/consensus.fasta",        sample=SAMPLES),
        # Use manifest of dereplicated, polished bins instead of MetaBAT2-only bins
        manifests = expand(f"{OUTDIR}/binning/{{sample}}/polished_bins/manifest.tsv", sample=SAMPLES),
        checkm2   = expand(f"{OUTDIR}/mag_qc/{{sample}}/checkm2.tsv",         sample=SAMPLES),
        gtdb      = expand(f"{OUTDIR}/taxonomy/{{sample}}/gtdbtk.tsv",        sample=SAMPLES)
    output:
        csv = f"{OUTDIR}/summary/metrics/summary.csv"
    threads: 2
    conda: "../envs/summary.yaml"
    params:
        samples = SAMPLES,
        outdir  = f"{OUTDIR}/summary/metrics"
    script:
        "../scripts/summarize_metrics.py"
