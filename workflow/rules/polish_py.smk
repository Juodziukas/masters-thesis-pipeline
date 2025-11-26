# First perform Racon polishing on the assembly contigs.  This rule runs
# a configurable number of rounds of Racon (default: 3) to correct raw
# long‑read errors before Medaka.  The result is written under
# results/polish/<sample>/racon/consensus.fasta.
rule racon_polish_py:
        input:
            draft = f"{OUTDIR}/assembly/{{sample}}/contigs.fasta",
            reads = f"{OUTDIR}/trim/{{sample}}.fastq.gz"
        output:
            polished = f"{OUTDIR}/polish/{{sample}}/racon/consensus.fasta"
        benchmark:
            f"benchmarks/polish_racon_{{sample}}.tsv"
        threads: THREADS
        conda: "../envs/racon.yaml"
        params:
            iterations = int(config.get("polish", {}).get("racon_rounds", 3)),
            outdir = f"{OUTDIR}/polish/{{sample}}/racon"
        script:
            "../scripts/run_racon.py"

# After Racon has corrected systematic ONT errors, run Medaka to further
# polish the consensus.  The input to Medaka is the Racon output,
# producing the final polished consensus at
# results/polish/<sample>/medaka/consensus.fasta.
rule medaka_polish_py:
        input:
            draft = f"{OUTDIR}/polish/{{sample}}/racon/consensus.fasta",
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
