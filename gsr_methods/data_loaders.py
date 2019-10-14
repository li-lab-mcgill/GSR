import os
import pandas as pd
import numpy as np
from pandas_plink import read_plink


def get_gene_to_chromo_dict(gene_chromo_dir):
    gene_to_chromo = {}

    with open(gene_chromo_dir) as gene_chromo_file:
        for l in gene_chromo_file:
            l = l.rstrip().split(',')
            if l[1] not in ['X', 'Y', 'MT']:
                gene_to_chromo[l[0].split('.')[0]] = int(l[1])

    return gene_to_chromo


def get_gene_map(gene_map_dir):
    entrez_to_ensemble = {}
    ensemble_to_entrez = {}

    with open(gene_map_dir) as gene_map_file:
        gene_map_file.readline()
        for l in gene_map_file:
            l = l.split()
            entrez_to_ensemble[l[1]] = l[0]
            ensemble_to_entrez[l[0]] = l[1]

    return entrez_to_ensemble, ensemble_to_entrez


def get_pathway_dict(pathway_dir, entrez_to_ensemble_df, gene_ensemble_id_list, reference_pathways=None):
    pathway_dict = {}

    with open(pathway_dir) as pathway_file:
        for l in pathway_file:
            ensemble_gene_ids = []
            l = l.split()
            for entrez_id in l[2:]:
                # select only only those entrez_ids that have the corresponding ensemble id
                if entrez_id in entrez_to_ensemble_df.index:
                    ensemble_id = entrez_to_ensemble_df.loc[entrez_id, "ensemble_id"]
                    # select only only those ensemble ids that have the corresponding weights
                    if ensemble_id in gene_ensemble_id_list:
                        ensemble_gene_ids.append(ensemble_id)

            if reference_pathways is None:
                if len(ensemble_gene_ids) >= 5:
                    pathway_dict[l[0]] = ensemble_gene_ids  # the gene id is in string format
            else:
                if l[0] in reference_pathways.keys():
                    pathway_dict[l[0]] = ensemble_gene_ids

    return pathway_dict


def get_chromo_snp_dict(thousand_G_dir):
    chromo_snp_dict = {}

    for i in range(1, 23):
        chromo_dir = os.path.join(thousand_G_dir, "1000G.EUR.{}".format(i))
        (bim, fam, bed) = read_plink(chromo_dir, verbose=False)
        chromo_snp = np.array(bim['snp'])
        X = bed.compute().T  # columns as SNP and row as number of individuals
        X_df = pd.DataFrame(data=X, columns=chromo_snp)
        chromo_snp_dict[i] = X_df

    return chromo_snp_dict