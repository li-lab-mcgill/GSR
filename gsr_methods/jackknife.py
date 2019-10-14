import pandas as pd
import numpy as np
from gsr_methods import data_loaders as dl
from gsr_methods.ldgcc import perform_ldgcc


def perform_jackknife(regression_result_df, blocks_list, pathway_dir, entrez_to_ensemble_df, gene_ensemble_id_list, pathway_dict, A_hat_df,
					  chi_square_df, sum_stats_df):

	n_blocks = len(blocks_list)
	jackknife_df = regression_result_df["tau_c"]

	for b in range(n_blocks):
		blk = blocks_list[b]
		filtered_gene_ensemble_id_list = gene_ensemble_id_list[~np.isin(gene_ensemble_id_list, blk)]
		filtered_pathway_dict = dl.get_pathway_dict(pathway_dir, entrez_to_ensemble_df,filtered_gene_ensemble_id_list,
													reference_pathways=pathway_dict)
		filtered_genes = np.array(list(filtered_pathway_dict.values()))
		filtered_genes = np.concatenate(filtered_genes)
		filtered_gene_ensemble_id_list = gene_ensemble_id_list[np.isin(gene_ensemble_id_list,
																	   filtered_genes)]
		filtered_regression_results_df = perform_ldgcc(filtered_pathway_dict, A_hat_df,
													   filtered_gene_ensemble_id_list,
													   chi_square_df, sum_stats_df)
		filtered_regression_results_df.rename(columns={"tau_c": f"tau_c_{b}"}, inplace=True)

		jackknife_df = pd.concat([jackknife_df, filtered_regression_results_df[f"tau_c_{b}"]],
								 axis=1, sort=False)


	est = np.array([jackknife_df["tau_c"].values])
	blk_est = jackknife_df[jackknife_df.columns[1:]].values.T
	pseudovalues = n_blocks * est - (n_blocks - 1) * blk_est

	jknife_cov = np.atleast_2d(np.cov(pseudovalues.T, ddof=1) / n_blocks)
	jknife_var = np.atleast_2d(np.diag(jknife_cov))
	jknife_se = np.atleast_2d(np.sqrt(jknife_var))
	jknife_est = np.atleast_2d(np.mean(pseudovalues, axis=0))

	z_scores = np.squeeze(jknife_est) / np.squeeze(jknife_se)
	jackknife_result_df = pd.DataFrame(data={"Z-scores": z_scores}, index=jackknife_df.index)

	return jackknife_result_df
