#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import argparse as ap
from sklearn.preprocessing import scale
from sklearn.naive_bayes import MultinomialNB
from sklearn.model_selection import StratifiedKFold
from scipy.stats import mannwhitneyu

def seq_to_kmerdict(df, k):
    df_sub_copy = df[["Proband", "sequence", "GC_100bp"]].copy()
    seqs = list(df_sub_copy["sequence"])

    kmer_dict = {} # Key: kmer, value: array of counts in each sequence
    for i in range(len(seqs)):
        seq = seqs[i]
        for j in range(len(seq)):
            kmer = seq[j:j + k]
            if len(kmer) != k:
                break
            if kmer in kmer_dict:
                kmer_dict[kmer][i] += 1
            else:
                kmer_dict[kmer] = np.zeros(len(seqs))
                kmer_dict[kmer][i] += 1

    all_kmers = []
    for kmer in kmer_dict:
        df_sub_copy[kmer] = kmer_dict[kmer]
        all_kmers.append(kmer)
    return df_sub_copy, all_kmers


def find_k_closest(gene, variants, k, gene_matrix):
    # Find the k variants closest to target gene
    cnt = 0
    top_genes = gene_matrix.loc[gene].sort_values(ascending=False)
    top_genes_list = [item for item in top_genes.index]
    current_gene_index = 0
    variants_index = []
    while cnt < k:
        current_gene = top_genes_list[current_gene_index]
        variants_with = variants[variants["GeneNAMEOuterTSS"] == current_gene]
        if cnt + len(variants_with) <= k:
            variants_index += list(variants_with.index)
        else:
            variants_index += list(np.random.choice(variants_with.index, size=k - cnt, replace=False))
        cnt += len(variants_with)
        current_gene_index += 1
    variants_index = [item for item in variants_index]
    return variants.loc[variants_index], top_genes[current_gene]


def model_2fold(X, y, gcs, seed=None):
    cv = StratifiedKFold(n_splits=2, shuffle=True, random_state=seed)
    clf = MultinomialNB(fit_prior=False)
    y_preds = np.zeros(len(y))
    fold_ps = []
    fold_ups = []
    fold_gcs = []
    fold_ugcs = []
    for train_index, test_index in cv.split(X, y):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        # clf is the naive bayes classifier
        clf.fit(X_train, y_train)

        y_pred_fold = clf.predict_proba(X_test)[:, 1]
        y_preds[test_index] = y_pred_fold
        proband_scores = list(y_pred_fold[y_test == True])
        sibling_scores = list(y_pred_fold[y_test == False])
        u, p = mannwhitneyu(proband_scores, sibling_scores, alternative="greater")
        es = u / (len(proband_scores) * len(sibling_scores))
        fold_ps.append(p)
        fold_ups.append(es)

        gc_test = gcs[test_index]
        proband_scores = list(gc_test[y_test == True])
        sibling_scores = list(gc_test[y_test == False])
        u, p = mannwhitneyu(proband_scores, sibling_scores, alternative="greater")
        es = u / (len(proband_scores) * len(sibling_scores))
        fold_gcs.append(p)
        fold_ugcs.append(es)
    return fold_ps, fold_ups, fold_gcs, fold_ugcs


parser = ap.ArgumentParser()
parser.add_argument('-v', '--variants', help="Variants file", nargs='?')
parser.add_argument('-g', '--genes', help="Gene matrix file", nargs='?')
parser.add_argument('-o', '--output', help="Output path", nargs='?')
parser.add_argument('-n', '--neighbors', help="Size of each neighborhood, default 1000", nargs='?', type=int, default=1000, const=1000)
parser.add_argument('-k', '--kmer', help="Length of k-mers, default 6", nargs='?', type=int, default=6, const=6)


args = parser.parse_args()
output = args.output
k = args.kmer
n_neighbors = args.neighbors
gene_file = args.genes
variants = args.variants

sys.stdout = open(output, "w+")

# np.random.seed(644)

gene_matrix = pd.read_csv(gene_file, sep='\t', index_col=0)
df = pd.read_csv(variants, index_col=0).reset_index(drop=True)
df_kmers, features = seq_to_kmerdict(df, k)

all_genes = gene_matrix.index

sys.stdout.write("Gene\tGC_u\tGC_p\tKmer_u_fold\tKmer_p_fold\tGC_u_fold\tGC_p_fold\n")
for gene in all_genes:
    try:
        neighbors, corr = find_k_closest(gene, df, n_neighbors, gene_matrix)
    except KeyError:
        continue
    train_data = df_kmers.iloc[list(neighbors.index)].copy()
    output_line = [gene]
    X = train_data[features].values
    y = train_data["Proband"].values
    gcs = train_data["GC_100bp"].values
    
    # GC MWU
    proband_gcs = train_data[train_data["Proband"] == True]["GC_100bp"]
    sibling_gcs = train_data[train_data["Proband"] == False]["GC_100bp"]
    u, p = mannwhitneyu(proband_gcs, sibling_gcs, alternative="greater")
    es = u / (len(proband_gcs) * len(sibling_gcs))
    output_line += [es, p]

    # 2-fold GC vs. K-mers
    fold_ps, fold_ups, fold_gcs, fold_ugcs = model_2fold(X, y, gcs, seed=np.random.randint(0, 9999999))
    output_line += [fold_ups[0], fold_ps[0]]
    output_line += [fold_ugcs[0], fold_gcs[0]]

    sys.stdout.write('\t'.join([str(item) for item in output_line]) + '\n')
