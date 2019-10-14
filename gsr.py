import pandas as pd
import numpy as np
import os
import argparse
from gsr_methods import data_loaders as dl
from gsr_methods.ldgcc import perform_ldgcc, get_A_hat, get_chi_squares


# defining the arguments
parser = argparse.ArgumentParser(description="LDGSR Manual", formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 argument_default=argparse.SUPPRESS)
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument("--sumstats", help="directory for summary statistic")
optional.add_argument("--pathway", help="directory for gene-pathway information",
                      default=os.path.join("data", "msigdb", "msigBIOCARTA_KEGG_REACTOME.gmt"))
optional.add_argument("--gene-map", help="directory for the entrez to ensemble ID",
                      default=os.path.join("data", "msigdb", "entrez2ensembl_without_hla.txt"))
optional.add_argument("--gene_chrom", help="directory for gene chromosomal location information",
                      default=os.path.join("data", "gene_annotation", "gene_to_chromosome.csv"))
optional.add_argument("--reference", help="directory for reference genome folder",
                      default=os.path.join("data", "LDREF"))
optional.add_argument("--tissue", help="directory for gene expression files",
                      default=os.path.join("data", "tissues", "Whole_Blood"))
optional.add_argument("--cs-cal", help="the way to calculate chi-square (sum or max)",
                      default="sum")
optional.add_argument("--output", help="output directory",
                      default=os.path.join("./"))


if __name__ == '__main__':
    try:
        args = parser.parse_args()

        if not os.path.exists(args.sumstats):
            raise Exception("[GSR]: Provided summary statistic does not exists.")
        if not os.path.exists(args.pathway):
            raise Exception("[GSR]: Provided pathway file does not exists.")
        if not os.path.exists(args.gene_map):
            raise Exception("[GSR]: Provided gene code file does not exists.")
        if not os.path.exists(args.gene_chrom):
            raise Exception("[GSR]: Provided gene chromosome files do not exists.")
        if not os.path.exists(args.reference):
            raise Exception("[GSR]: Provided reference genome does not exists.")
        if not os.path.exists(args.tissue):
            raise Exception("[GSR]: Provided gene weight files do not exists.")
        if not os.path.exists(args.output):
            raise Exception("Provided output directory does not exists.")
        if args.cs_cal != "sum" and args.cs_cal != "max":
            raise Exception("[GSR]: Chi-Square is currently only calculated by sumation or maximization.")

        sum_stats_dir = args.sumstats
        pathway_dir = args.pathway
        gene_map_dir = args.gene_map
        gene_chromo_dir = args.gene_chrom
        thousand_G_dir = args.reference
        tissue_gene_dir = args.tissue
        output_dir = args.output
        cs_cal = args.cs_cal

    except Exception as e:
        print(e, end="\n\n")
        parser.parse_args(['-h'])

    else:
        try:
            print("[GSR]: start loadling necessary data.")
            # loading the summary statistic
            sum_stats_df = pd.read_csv(sum_stats_dir, delim_whitespace=True, index_col=0)
            sum_stats_df = sum_stats_df.loc[np.isfinite(sum_stats_df["Z"])]

            # loading gene chromosomal location information into a dictionary
            gene_chromo_dir = os.path.join("data", "gene_annotation", "gene_to_chromosome.csv")
            gene_to_chromo = dl.get_gene_to_chromo_dict(gene_chromo_dir)

            # loading genes indexing information into a dictionary
            entrez_to_ensemble, _ = dl.get_gene_map(gene_map_dir)
            entrez_to_ensemble_df = pd.DataFrame.from_dict(entrez_to_ensemble, orient='index', columns=["ensemble_id"])

            # getting all the SNP information by reading the plink files
            chromo_snp_dict = dl.get_chromo_snp_dict(thousand_G_dir)
            print("[GSR]: finish loadling necessary data.")

            # calculating A_hat
            A_hat_df, all_SNPs_df, gene_ensemble_id_list = get_A_hat(gene_to_chromo, entrez_to_ensemble_df,
                                                                     chromo_snp_dict, tissue_gene_dir, sum_stats_df)

            # calculating chi_square
            chi_square_df = get_chi_squares(gene_ensemble_id_list, tissue_gene_dir, cs_cal, A_hat_df, all_SNPs_df,
                                            sum_stats_df)

            # loading pathway-genes information into a dictionary
            pathway_dict = dl.get_pathway_dict(pathway_dir, entrez_to_ensemble_df, gene_ensemble_id_list)

            regression_results_df = perform_ldgcc(pathway_dict, A_hat_df, gene_ensemble_id_list,
                                                  chi_square_df, sum_stats_df)

            sum_stats = sum_stats_dir.split('/')[-1].replace('.', '_')
            file_to_save = "{}_regression_results.csv".format(sum_stats_dir.split('/')[-1])
            regression_results_df.sort_values(by=["t-value"], ascending=False).to_csv(os.path.join(output_dir, file_to_save))

        except Exception as e:
            print("[GSR]: Unexpected error occur: ", end="")
            print(e)

        finally:
            print("[GSR]: program finishes running!")
