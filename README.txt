Gene Score Regression (GSR) User Guid

Programming language requirement
    Python 3.7 + 
  
Python package
    numpy 1.16 + 
    statsmodels 0.10 +
    Pandas 0.25 +
    pandas_plink 2.0 +

GSR commands
  The general commands for using gsr is 
    Python gsr.py --sumstats SUMSTATS --pathway PATHWAY --gene-map GENE_MAP --gene_chrom GENE_CHROM --reference REFERENCE --tissue TISSUE --cs-cal CS_CAL --output OUTPUT 

  Requirement argument
    --sumstats SUMSTATS
        "SUMSTATS" is the directory for the input summary statistic. For example "--sumstats data/sumstats_formatted/PASS_LDL.sumstats". The .sumstats files need to be a tab seperated file, and it has to include the column "SNP", "N", "CHISQ", and "Z". Please refer to one of the .sumstats file in data/sumstats_formatted for reference.

  Optional arguments
    --pathway PATHWAY
        "PATHWAY" is the file directory for gene pathway or gene set. By default, gsr used the provided file in "data/msigdb/msigBIOCARTA_KEGG_REACTOME.gmt". The gene set file is tab seperated file, and it follows "NAME_OF_THE_SET\tINFORMATION\tENTREZ_GENE_ID1\tENTREZ_GENE_ID2...". Please refer to "msigBIOCARTA_KEGG_REACTOME.gmt" for reference.

    --gene-map GENE_MAP
        "GENE_MAP" is the file directory for a list of entrez IDs and their correpond ENSEMBLE ID. By default, gsr used the provided file in "data/msigdb/entrez2ensembl_without_hla.txt". The file is space seperated where the first row is "ensembl_gene_id entrezgene". Please refer to "entrez2ensembl_without_hla.txt" for reference.

    --gene_chrom GENE_CHROM
        "GENE_CHROM" is the file directory for the information for which chromosome the gene is found. By default, gsr used the provided file in "data/gene_annotation/gene_to_chromosome.csv". Please refer to "gene_to_chromosome.csv" for reference.

    --reference REFERENCE
        "REFERENCE" is the folder directory containing a reference population whole genotype except the sex chromosome genotype. By default, gsr used the provided folder in "data/LDREF" which is the 1KG project genotype. Every chromosome information files in the folder needs to follow a plink format, which means it needs to have ".bed, .bim, fam" for 1 chomosome genotype file. Please to "LDREF/" for reference.

    --tissue TISSUE
        "TISSUE" is the folder directory containing tissue TWAS gene information. By default, gsr used the provided folder in "data/tissues/Whole_Blood". The nomination of the files in the folder needs to follow the format "TISSUE_NAME.ENSEMBLE_ID.VER.wgt.csv", and each files contain two columns "SNP" and "estimated from TWAS". For example the file "Whole_Blood.ENSG00000273478.1.wgt.csv" contain the SNP and their TWAS weight for ENSG00000273478 gene expression in whole blood. Please to "Whole_Blood/" for reference.

    --cs-cal CS_CAL
        "CS_CAL" is the method for calculating the gene chi square. It can either be "sum" or "max". By default is set to be "sum", e.g. "--cs-cal sum".

    --output OUTPUT
        "OUTPUT" is the out directory for the gene regression result. By default is the currently directory where gsr.py is located, e.g. "--output ./"

    The arguments can follows any order. User can always use "Python gsr.py -h" for a quick overview of the software's arguments. 


Example of commands for quick start
    "python gsr.py --sumstats data/sumstats_formatted/PASS_LDL.sumstats". A file named 'PASS_LDL_regression_results.csv' will be generated in the current directory.

    "python gsr.py --sumstats data/sumstats_formatted/PASS_Schizophrenia.sumstats --pathway data/msigdb/c5.bp.v6.2.entrez.gmt --cs-cal max --output ~/gsr_result/". In this command, GSR will use the `c5.bp.v6.2.entrez.gmt` file for gene set information, and a fille named "PASS_Schizophrenia_regression_results.csv" will be generated in the `~/gsr_result/` directory.







