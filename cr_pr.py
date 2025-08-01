#!/usr/bin/env python3

import argparse
import pandas as pd
from collections import defaultdict

# Define genotype categories
HET_GENOS = {"0/1", "0|1", "1/0", "1|0", "./1", ".|1", "1/.", "1|."}
HOM_GENOS = {"1/1", "1|1"}

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate Carrier and Prevalence Rates.")
    parser.add_argument("-gf", "--GTfile", required=True, help="Path to genotype matrix file")
    parser.add_argument("-qf", "--Queryfile", required=True, help="Path to query file")
    parser.add_argument("-o", "--Output", required=True, help="Output file path")
    return parser.parse_args()

def load_data(gt_file, query_file):
    gt_df = pd.read_csv(gt_file, sep="\t")
    query_df = pd.read_csv(query_file, sep="\t")
    return gt_df, query_df

def merge_and_align(gt_df, query_df):
    merged = pd.merge(query_df, gt_df, on="variant", how="inner")
    return merged

def is_het(geno): return geno in HET_GENOS
def is_hom(geno): return geno in HOM_GENOS

def count_raw_variant(row, samples):
    het_count = sum(is_het(row[s]) for s in samples)
    hom_count = sum(is_hom(row[s]) for s in samples)
    return pd.Series([het_count, hom_count])

def build_gene_dict(df, samples):
    gene_variants = defaultdict(list)
    for _, row in df.iterrows():
        gene_variants[row['gene']].append(row)
    return gene_variants

def unique_counts_per_gene_moi(df, samples):
    grouped = df.groupby(['gene', 'moi'])
    het_per_gene_moi = {}
    hom_per_gene_moi = {}
    for (gene, moi), group in grouped:
        het = {s: 0 for s in samples}
        hom = {s: 0 for s in samples}
        for _, row in group.iterrows():
            for s in samples:
                g = row[s]
                if is_het(g):
                    het[s] += 1
                elif is_hom(g):
                    hom[s] += 1

        uniq_het = sum(1 for s in samples if het[s] > 0 and hom[s] == 0)
        uniq_hom = sum(1 for s in samples if hom[s] > 0)
        for idx in group.index:
            het_per_gene_moi[idx] = uniq_het
            hom_per_gene_moi[idx] = uniq_hom
    return het_per_gene_moi, hom_per_gene_moi

def ar_carrier_diseased(df, samples):
    grouped = df[df['moi'] == 'AR'].groupby(['condition', 'moi'])
    uniq_het_rec = {}
    uniq_hom_rec = {}

    for (cond, moi), group in grouped:
        het = {s: 0 for s in samples}
        hom = {s: 0 for s in samples}
        for _, row in group.iterrows():
            for s in samples:
                g = row[s]
                if is_het(g):
                    het[s] += 1
                elif is_hom(g):
                    hom[s] += 1
        final_het = sum(1 for s in samples if het[s] > 0 and hom[s] == 0)
        final_hom = sum(1 for s in samples if hom[s] > 0)
        for idx in group.index:
            uniq_het_rec[idx] = final_het
            uniq_hom_rec[idx] = final_hom
    return uniq_het_rec, uniq_hom_rec

def ad_diseased(df, samples):
    grouped = df[df['moi'] == 'AD'].groupby(['condition', 'moi'])
    het_hom_dom = {}
    for (cond, moi), group in grouped:
        seen = {s: 0 for s in samples}
        for _, row in group.iterrows():
            for s in samples:
                g = row[s]
                if is_hom(g) or is_het(g):
                    seen[s] += 1
        count = sum(1 for s in samples if seen[s] > 0)
        for idx in group.index:
            het_hom_dom[idx] = count
    return het_hom_dom

def gene_level_counts(df, samples):
    het_gene = {}
    hom_gene = {}
    grouped = df.groupby("gene")
    for gene, group in grouped:
        het_total = 0
        hom_total = 0
        for _, row in group.iterrows():
            for s in samples:
                g = row[s]
                if is_het(g): het_total += 1
                if is_hom(g): hom_total += 1
        for idx in group.index:
            het_gene[idx] = het_total
            hom_gene[idx] = hom_total
    return het_gene, hom_gene

def main():
    args = parse_args()
    gt_df, query_df = load_data(args.GTfile, args.Queryfile)
    merged_df = merge_and_align(gt_df, query_df)
    samples = [col for col in merged_df.columns if col not in ['variant', 'gene', 'moi', 'condition']]

    merged_df[['het_raw_per_variant', 'hom_raw_per_variant']] = merged_df.apply(lambda row: count_raw_variant(row, samples), axis=1)
    
    het_gene, hom_gene = gene_level_counts(merged_df, samples)
    merged_df['het_raw_per_gene'] = merged_df.index.map(het_gene)
    merged_df['hom_raw_per_gene'] = merged_df.index.map(hom_gene)

    uniq_het_moi, uniq_hom_moi = unique_counts_per_gene_moi(merged_df, samples)
    merged_df['uniq_het_for_moi_per_gene'] = merged_df.index.map(uniq_het_moi)
    merged_df['uniq_hom_for_moi_per_gene'] = merged_df.index.map(uniq_hom_moi)

    uniq_het_rec, uniq_hom_rec = ar_carrier_diseased(merged_df, samples)
    merged_df['uniq_het_for_rec_pmoi_CARRIER'] = merged_df.index.map(lambda x: uniq_het_rec.get(x, "."))
    merged_df['uniq_hom_for_rec_pmoi_DISEASED'] = merged_df.index.map(lambda x: uniq_hom_rec.get(x, "."))

    het_hom_dom = ad_diseased(merged_df, samples)
    merged_df['uniq_hethom_for_dom_pmoi_DISEASED'] = merged_df.index.map(lambda x: het_hom_dom.get(x, "."))

    merged_df.to_csv(args.Output, sep="\t", index=False)

if __name__ == "__main__":
    main()
