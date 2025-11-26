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
    # Run MetaBAT2 on filtered contigs to generate initial bins
    input:
        fasta = f"{OUTDIR}/binning/{{sample}}/filtered/contigs.fasta",
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

# Filter contigs before binning: remove short contigs (<2.5 kb), low‑coverage
# contigs (<5×), and extreme GC (<20 % or >80 %).  This greatly improves bin
# quality and reduces contamination.
rule filter_contigs_py:
        input:
            fasta = f"{OUTDIR}/polish/{{sample}}/medaka/consensus.fasta",
            depth = f"{OUTDIR}/binning/{{sample}}/jgi_depth.txt"
        output:
            filtered = f"{OUTDIR}/binning/{{sample}}/filtered/contigs.fasta"
        threads: 2
        conda: "../envs/filter.yaml"
        params:
            min_length   = int(config.get("filter", {}).get("min_length", 2500)),
            min_coverage = float(config.get("filter", {}).get("min_coverage", 5)),
            min_gc       = float(config.get("filter", {}).get("min_gc", 0.2)),
            max_gc       = float(config.get("filter", {}).get("max_gc", 0.8)),
            outdir       = f"{OUTDIR}/binning/{{sample}}/filtered"
        script:
            "../scripts/run_filter_contigs.py"

# Run MaxBin2 to obtain an alternative bin set.  Different binning
# algorithms complement each other and catch genomes missed by MetaBAT2.
rule maxbin2_py:
        input:
            fasta = f"{OUTDIR}/binning/{{sample}}/filtered/contigs.fasta",
            depth = f"{OUTDIR}/binning/{{sample}}/jgi_depth.txt"
        output:
            directory(f"{OUTDIR}/binning/{{sample}}/maxbin2"),
            f"{OUTDIR}/binning/{{sample}}/maxbin2/manifest.tsv"
        conda: "../envs/maxbin2.yaml"
        threads: 4
        params:
            outdir = lambda wildcards: f"{OUTDIR}/binning/{wildcards.sample}/maxbin2"
        script:
            "../scripts/run_maxbin2.py"

# Run SemiBin2 (deep learning based binning) for long‑read data.
rule semibin2_py:
        input:
            fasta = f"{OUTDIR}/binning/{{sample}}/filtered/contigs.fasta",
            depth = f"{OUTDIR}/binning/{{sample}}/jgi_depth.txt"
        output:
            directory(f"{OUTDIR}/binning/{{sample}}/semibin2"),
            f"{OUTDIR}/binning/{{sample}}/semibin2/manifest.tsv"
        conda: "../envs/semibin2.yaml"
        threads: 4
        params:
            outdir = lambda wildcards: f"{OUTDIR}/binning/{wildcards.sample}/semibin2"
        script:
            "../scripts/run_semibin2.py"

# Combine bins from MetaBAT2, MaxBin2 and SemiBin2 using DAS Tool.  This
# selects the best bins from each algorithm and produces a refined bin set.
rule das_tool_py:
        input:
            metabat_manifest = f"{OUTDIR}/binning/{{sample}}/metabat2/manifest.tsv",
            maxbin_manifest  = f"{OUTDIR}/binning/{{sample}}/maxbin2/manifest.tsv",
            semibin_manifest = f"{OUTDIR}/binning/{{sample}}/semibin2/manifest.tsv"
        output:
            directory(f"{OUTDIR}/binning/{{sample}}/das_tool"),
            f"{OUTDIR}/binning/{{sample}}/das_tool/manifest.tsv"
        conda: "../envs/das_tool.yaml"
        threads: 4
        params:
            outdir = lambda wildcards: f"{OUTDIR}/binning/{wildcards.sample}/das_tool"
        script:
            "../scripts/run_das_tool.py"

# Dereplicate/refine the DAS Tool bins with dRep to remove overlapping or
# redundant genomes.
rule drep_py:
        input:
            bins_dir = f"{OUTDIR}/binning/{{sample}}/das_tool"
        output:
            directory(f"{OUTDIR}/binning/{{sample}}/drep"),
            f"{OUTDIR}/binning/{{sample}}/drep/manifest.tsv"
        conda: "../envs/drep.yaml"
        threads: 4
        params:
            outdir = lambda wildcards: f"{OUTDIR}/binning/{wildcards.sample}/drep"
        script:
            "../scripts/run_drep.py"

# Perform per‑bin polishing (Racon+Medaka) on each dereplicated bin.
rule polish_bins_py:
        input:
            manifest = f"{OUTDIR}/binning/{{sample}}/drep/manifest.tsv",
            reads    = f"{OUTDIR}/trim/{{sample}}.fastq.gz"
        output:
            directory(f"{OUTDIR}/binning/{{sample}}/polished_bins"),
            f"{OUTDIR}/binning/{{sample}}/polished_bins/manifest.tsv"
        conda: "../envs/polish_bin.yaml"
        threads: THREADS
        params:
            outdir = lambda wildcards: f"{OUTDIR}/binning/{wildcards.sample}/polished_bins"
        script:
            "../scripts/run_polish_bin.py"
