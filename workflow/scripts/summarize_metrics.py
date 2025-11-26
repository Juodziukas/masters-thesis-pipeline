from pathlib import Path
from collections import defaultdict
import csv
import math

try:
    import pandas as pd
except ImportError:
    pd = None

def fasta_lengths(fa_path):
    lens = []
    seq_len = 0
    with open(fa_path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                if seq_len > 0:
                    lens.append(seq_len)
                    seq_len = 0
            else:
                seq_len += len(line.strip())
    if seq_len > 0:
        lens.append(seq_len)
    return lens

def n50(lengths):
    if not lengths:
        return 0
    s = sorted(lengths, reverse=True)
    total = sum(s)
    cum = 0
    for L in s:
        cum += L
        if cum >= total/2:
            return L
    return 0

def sample_from_path(path_str):
    p = Path(path_str)
    parts = p.parts
    # find 'results' and take next element as stage, then next as sample
    if "results" in parts:
        i = parts.index("results")
        if i + 2 < len(parts):
            return parts[i+2]
    return p.parent.name

samples = list(snakemake.params.samples)   # type: ignore
outdir = Path(str(snakemake.params.outdir))# type: ignore
out_csv = Path(str(snakemake.output.csv))  # type: ignore

outdir.mkdir(parents=True, exist_ok=True)

# --- Assembly stats ---
asm_rows = []
for fa in snakemake.input.contigs:  # type: ignore
    s = sample_from_path(str(fa))
    lens = fasta_lengths(str(fa)) if Path(fa).exists() else []
    total = sum(lens)
    N50 = n50(lens)
    asm_rows.append({"sample": s, "assembly_bp": total, "n50": N50, "contigs": len(lens)})

# --- Binning: count bins per sample ---
bins_count = defaultdict(int)
for mf in snakemake.input.manifests:  # type: ignore
    s = sample_from_path(str(mf))
    c = 0
    if Path(mf).exists():
        with open(mf) as fh:
            rdr = csv.reader(fh, delimiter="\t")
            header = next(rdr, None)
            for _ in rdr:
                c += 1
    bins_count[s] = c

# --- CheckM2: summarize completeness/contamination, HQ/LQ ---
# MIMAG-ish cutoffs: HQ >=90% comp & <=5% contam; MQ >=50% & <=10% contam
cm2 = defaultdict(lambda: {"mean_comp": math.nan, "mean_contam": math.nan,
                           "hq_mag": 0, "mq_mag": 0, "total_bins": 0})

# Consolidate multiple files per sample (usually one)
cm2_by_sample = defaultdict(list)
for tsv in snakemake.input.checkm2:  # type: ignore
    s = sample_from_path(str(tsv))
    cm2_by_sample[s].append(tsv)

for s, paths in cm2_by_sample.items():
    comps, conts = [], []
    hq, mq, tot = 0, 0, 0
    for tsv in paths:
        if not Path(tsv).exists():
            continue
        with open(tsv) as fh:
            # try to detect header columns
            header = fh.readline().strip().split("\t")
            # normalize column names to lower
            cols = {name.lower(): i for i, name in enumerate(header)}
            # expected names (allow variants)
            comp_key = next((k for k in cols if "completeness" in k), None)
            con_key  = next((k for k in cols if "contamination" in k), None)
            # iterate rows
            for line in fh:
                if not line.strip(): continue
                parts = line.rstrip("\n").split("\t")
                try:
                    comp = float(parts[cols[comp_key]]) if comp_key else float("nan")
                    con  = float(parts[cols[con_key]])  if con_key  else float("nan")
                except Exception:
                    continue
                comps.append(comp); conts.append(con); tot += 1
                if comp >= 90.0 and con <= 5.0:
                    hq += 1
                elif comp >= 50.0 and con <= 10.0:
                    mq += 1
    mcomp = sum(comps)/len(comps) if comps else float("nan")
    mcont = sum(conts)/len(conts) if conts else float("nan")
    cm2[s]["mean_comp"] = mcomp
    cm2[s]["mean_contam"] = mcont
    cm2[s]["hq_mag"] = hq
    cm2[s]["mq_mag"] = mq
    cm2[s]["total_bins"] = tot

# --- GTDB: count classified genomes per sample ---
gtdb_counts = defaultdict(int)
for tsv in snakemake.input.gtdb:  # type: ignore
    s = sample_from_path(str(tsv))
    c = 0
    if Path(tsv).exists():
        with open(tsv) as fh:
            rdr = csv.reader(fh, delimiter="\t")
            header = next(rdr, None)
            for _ in rdr:
                c += 1
    gtdb_counts[s] = c

# --- Merge rows per sample ---
rows = []
for s in samples:
    a = next((r for r in asm_rows if r["sample"] == s), {"assembly_bp": 0, "n50": 0, "contigs": 0})
    row = {
        "sample": s,
        "assembly_bp": a["assembly_bp"],
        "n50": a["n50"],
        "contigs": a["contigs"],
        "bins_metabat2": bins_count.get(s, 0),
        "checkm2_mean_completeness": cm2[s]["mean_comp"],
        "checkm2_mean_contamination": cm2[s]["mean_contam"],
        "checkm2_HQ_MAGs": cm2[s]["hq_mag"],
        "checkm2_MQ_MAGs": cm2[s]["mq_mag"],
        "checkm2_total_bins": cm2[s]["total_bins"],
        "gtdb_classified": gtdb_counts.get(s, 0),
    }
    rows.append(row)

out_csv.parent.mkdir(parents=True, exist_ok=True)
with open(out_csv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["sample"])
    w.writeheader()
    for r in rows:
        w.writerow(r)
