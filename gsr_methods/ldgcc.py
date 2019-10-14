import csv
import pandas as pd
import numpy as np
import fnmatch
import os
import statsmodels.api as sm


# this method help to get Wg and the corresponding SNP
def get_Wg_and_SNPs(gene_dir):
    with open(gene_dir, 'r') as csvfile:
        lines = [l for l in csvfile][:2]
        text_sample = "".join(lines)
        if len(lines) < 2:
            has_header = False
        else:
            has_header = csv.Sniffer().has_header(text_sample)

    if has_header:
        print(gene_dir)
        W_df = pd.read_csv(gene_dir, index_col=0)
    else:
        W_df = pd.read_csv(gene_dir, index_col=0, header=None, names=["w_hat"])

    non_zero_W_g = W_df["w_hat"].values
    selected_SNPs = W_df["w_hat"].index

    return non_zero_W_g, selected_SNPs


# this is the method to get Ag_hat for a given gene g
def get_Ag_hat(gene_ensemble_id, gene_dir, chromo_snp_dict, gene_to_chromo):
    non_zero_W_g, selected_SNPs = get_Wg_and_SNPs(gene_dir)

    X = chromo_snp_dict[gene_to_chromo[gene_ensemble_id]]   # get the chromosome where the gene located
    X = X[selected_SNPs].values   # obtain the corresponding SNP values
    standardized_X = (X - X.mean(axis=0)) / X.std(axis=0)  # standardize X

    Ag_hat = np.squeeze(np.matmul(standardized_X, non_zero_W_g))
    standardized_Ag_hat = (Ag_hat - Ag_hat.mean(axis=0)) / Ag_hat.std(axis=0)

    return  standardized_Ag_hat, selected_SNPs


def get_A_hat(gene_to_chromo, entrez_to_ensemble_df, chromo_snp_dict, tissue_gene_dir, sum_stats_df):
    # calculating A_hat
    global gene_ensemble_id
    print("[GSR]: start calculating A_hat.")

    N_1KG = chromo_snp_dict[1].shape[0]
    all_tissue_gene_files = fnmatch.filter(os.listdir(tissue_gene_dir), '*')
    A_hat = []
    gene_ensemble_id_list = []  # list of gene that has the corresponding weight and a corresponding entrez id
    all_SNPs_set = set()

    for g, gene_file in enumerate(all_tissue_gene_files):
        for gid in gene_file.split('.'):
            if "ENSG" in gid:
                gene_ensemble_id = gid
                break

        if gene_ensemble_id in entrez_to_ensemble_df["ensemble_id"].values:
            gene_dir = os.path.join(tissue_gene_dir, gene_file)
            Ag_hat, selected_SNPs = get_Ag_hat(gene_ensemble_id, gene_dir, chromo_snp_dict, gene_to_chromo)

            if not np.isnan(Ag_hat).any():
                gene_ensemble_id_list.append(gene_ensemble_id)
                A_hat.append(np.squeeze(Ag_hat))  # the g_th colum in A_hat is basically Ag_hat
                all_SNPs_set.update(selected_SNPs)  # saving all possible SNPs for indexing z_hat later

    A_hat = np.array(A_hat).T
    gene_ensemble_id_list = np.array(gene_ensemble_id_list)
    A_hat_df = pd.DataFrame(data=A_hat, columns=gene_ensemble_id_list)  # gene_ensemble_id_list has the same other as all_tissue_gene_files

    default_value = np.zeros(len(all_SNPs_set))
    all_SNPs_set = list(all_SNPs_set)
    all_SNPs_df = pd.DataFrame(data=default_value, columns=["Z"], index=all_SNPs_set)
    sum_stats_SNPs = sum_stats_df.index
    available_SNPs = np.intersect1d(sum_stats_SNPs, all_SNPs_df.index)
    all_SNPs_df.loc[available_SNPs, "Z"] = sum_stats_df.loc[available_SNPs, "Z"]

    print("[GSR]: finish calculating A_hat.")

    return A_hat_df, all_SNPs_df, gene_ensemble_id_list


def get_chi_squares(gene_ensemble_id_list, tissue_gene_dir, cs_cal, A_hat_df, all_SNPs_df, sum_stats_df):
    # calculating all the chi-square
    global Chi_g_square
    print("[GSR]: start calculating chi-square.")

    A_hat = A_hat_df.values
    N_1KG = A_hat_df.shape[0]
    Ngs = np.sum(A_hat ** 2, axis=0)
    N_total_population = sum_stats_df['N'].values[0]

    Chi_squares = np.zeros(len(gene_ensemble_id_list))
    gene_without_z_scores = []

    for g, gene_ensemble_id in enumerate(gene_ensemble_id_list):
        gene_file = fnmatch.filter(os.listdir(tissue_gene_dir), '*{}*'.format(gene_ensemble_id))[0]
        gene_dir = os.path.join(tissue_gene_dir, gene_file)

        non_zero_W_g, selected_SNPs = get_Wg_and_SNPs(gene_dir)  # get Wg and the SNP for gene g
        Zg_hat = all_SNPs_df.loc[selected_SNPs, 'Z'].values

        if cs_cal == "sum":
            Chi_g_square = (np.matmul(non_zero_W_g.T, Zg_hat) ** 2) * (N_1KG ** 2) / (
                    N_total_population * Ngs[g] ** 2)
        elif cs_cal == "max":
            Chi_g_square = (np.multiply(non_zero_W_g.squeeze(), Zg_hat).max() ** 2) * (N_1KG ** 2) / (
                    N_total_population * (Ngs[g] ** 2))

        if Chi_g_square == 0:
            gene_without_z_scores.append(gene_ensemble_id)

        Chi_squares[g] = Chi_g_square

    chi_square_df = pd.DataFrame(data=Chi_squares, columns=["Chi-Square"],
                                 index=gene_ensemble_id_list)  # gene_ensemble_id_list has the same other as all_tissue_gene_files
    chi_square_df = chi_square_df.drop(gene_without_z_scores)  # removing gene that has non zero square

    print("[GSR]: finish calculating chi-square.")

    return chi_square_df


def perform_ldgcc(pathway_dict, A_hat_df, gene_ensemble_id_list, chi_square_df, sum_stats_df):
    print("[GSR]: start doing ldgcc.")

    if A_hat_df.shape[1] != A_hat_df[gene_ensemble_id_list].shape[1]:
        A_hat_df = A_hat_df[gene_ensemble_id_list]

    A_hat = A_hat_df.values
    N_1KG = A_hat_df.shape[0]
    Ngs = np.sum(A_hat ** 2, axis=0)
    N_total_population = sum_stats_df['N'].values[0]

    pathway_list = list(pathway_dict.keys())
    ldgcc = np.zeros((A_hat.shape[1], len(pathway_dict) + 1))  # +1 for dummy pathway
    scaling_factor = N_total_population * (N_1KG ** 2) / (Ngs ** 2)
    empty_pathway = []

    # looping for each pathway c
    for c, pathway in enumerate(pathway_list):
        ensemble_gene_list = pathway_dict[pathway]

        selected_pathway_genes = np.intersect1d(gene_ensemble_id_list, ensemble_gene_list)

        if len(selected_pathway_genes) == 0:
            empty_pathway.append(pathway)
        else:
            pathway_correlation_matrix = np.matmul(A_hat.T, A_hat_df[selected_pathway_genes].values)  # A_hat_df[selected_pathway_genes] = A_hat_c

            l_c = scaling_factor * np.sum(pathway_correlation_matrix ** 2, axis=1) / (N_1KG ** 2)  # l(g) size = G x 1
            ldgcc[:, c] += l_c

    dummy_pathway = np.matmul(A_hat.T, A_hat)  # the dummy pathway has all the genes
    dummy_pathway_l_c = scaling_factor * np.sum(dummy_pathway ** 2, axis=1) / (N_1KG ** 2)  # calculating the l(g, c) where c is the dummy pathway
    ldgcc[:, -1] += dummy_pathway_l_c

    columns_list = pathway_list + ["Dummy Pathway"]
    ldgcc_df = pd.DataFrame(data=ldgcc, columns=columns_list, index=gene_ensemble_id_list)
    ldgcc_df = ldgcc_df.drop(columns=empty_pathway)


    # Start doing regression

    ldgcc_and_chi_square_df = pd.concat([ldgcc_df, chi_square_df], join="inner", axis=1)

    pathway_list = ldgcc_and_chi_square_df.columns[:-2]
    y = ldgcc_and_chi_square_df.values[:, -1]
    results_params = []
    results_tvalues = []
    results_pvalues = []

    for pathway in pathway_list:
        X = ldgcc_and_chi_square_df[[pathway, 'Dummy Pathway']].values
        X = sm.add_constant(X)  # add a constant column to absorb the bias
        pathway_results = sm.OLS(y, X).fit()
        results_params.append(pathway_results.params[1])
        results_tvalues.append(pathway_results.tvalues[1])
        results_pvalues.append(pathway_results.pvalues[1])

    regression_results_params_df = pd.DataFrame(data={"tau_c": results_params,
													  "t-value":results_tvalues,
													  "p-value":results_pvalues},
												index=pathway_list)

    print("[GSR]: finish doing ldgcc.")

    return regression_results_params_df

