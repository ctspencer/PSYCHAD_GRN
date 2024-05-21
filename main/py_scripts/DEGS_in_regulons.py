import pandas as pd
import numpy as np

def calculate_fdr_log10fdr(data, p_value_col='p.value'):
    """
    Calculate False Discovery Rate (FDR) and log10(FDR) from a pandas DataFrame.

    Args:
        data (pandas.DataFrame): DataFrame containing the data.
        p_value_col (str, optional): Name of the column containing p-values. Default is 'p.value'.

    Returns:
        pandas.DataFrame: Input DataFrame with two new columns 'FDR' and 'log10(FDR)' added.
    """
    # Sort the p-values in ascending order
    data = data.sort_values(by=p_value_col)

    # Rank the sorted p-values
    ranks = np.arange(len(data)) + 1

    # Calculate FDR using the Benjamini-Hochberg procedure
    fdr = data[p_value_col] * len(data) / ranks

    # Calculate log10(FDR)
    log10_fdr = -np.log10(fdr)

    # Add the FDR and log10(FDR) columns to the DataFrame
    data['FDR'] = fdr
    data['-log10(FDR)'] = log10_fdr

    return data

def filter_and_aggregate_genes(df, pc_gene_list, pvalue_threshold, coefficients, t_condition, annolevel='subclass', celltype='Micro'):
    """
    Filters and aggregates gene data based on specific criteria and statistical thresholds.
    
    Parameters:
    - df (DataFrame): The dataframe containing gene data with columns like 'p.value', 'coef', 'statistic', etc.
    - pc_gene_list (list): List of pre-defined genes of interest.
    - pvalue_threshold (float): Threshold for the p-value; genes with a p-value below this threshold are considered.
    - coefficients (list): A list of coefficients used to filter genes in the dataframe.
    - t_condition (str): Condition for filtering based on the statistic value. It can be '<=0', '>=0', or 'all'.
    - annolevel (str, optional): Annotation level used to further filter genes; defaults to 'subclass'.
    - celltype (str, optional): Specific cell type to further refine the filtering; defaults to 'Micro'.

    Returns:
    DataFrame: A dataframe containing aggregated results for common genes across the different filters, 
    including mean statistics and p-values, sorted by ascending p-value.
    """
    filtered_dfs = []
    for coef in coefficients:
        if t_condition == '<=0':
            filtered_df = df[(df['p.value'] < pvalue_threshold) & (df['coef'] == coef) & (df['statistic'] <= 0) & (df['ID'].isin(pc_gene_list)) & (df['AnnoLevel'] == annolevel) & (df['method'] == 'FE') & (df['assay'] == celltype)]
        elif t_condition == '>=0':
            filtered_df = df[(df['p.value'] < pvalue_threshold) & (df['coef'] == coef) & (df['statistic'] >= 0) & (df['ID'].isin(pc_gene_list) & (df['AnnoLevel'] == annolevel) & (df['method'] == 'FE')) & (df['assay'] == celltype)]
        elif t_condition == 'all':
            filtered_df = df[(df['p.value'] < pvalue_threshold) & (df['coef'] == coef) & (df['ID'].isin(pc_gene_list) & (df['AnnoLevel'] == annolevel) & (df['method'] == 'FE')) & (df['assay'] == celltype)]
        else:
            raise ValueError("Invalid t condition.")
        filtered_dfs.append(set(filtered_df['ID']))

    common_genes = set.intersection(*filtered_dfs)
    result_df = df[df['ID'].isin(common_genes)]
    result_df = result_df.groupby('ID', as_index=False).agg({'statistic':'mean', 'p.value': 'mean'}).sort_values(by='p.value', ascending=True)
    return result_df






def significant_targets_in_celltypes(GRN, TF_list, celltype_DEGs, quantile=0.0):
    # Initialize dictionary to store results
    TF_targets_dict = {}

    # Iterate through each cell type and its corresponding differentially expressed genes (DEGs)
    for celltype, DEGs in celltype_DEGs.items():
        # Prepare a sub-dictionary for the current cell type if not already present
        if celltype not in TF_targets_dict:
            TF_targets_dict[celltype] = {}
        
        # Iterate through each transcription factor (TF)
        for TF in TF_list:
            # Filter the data for the current TF
            tf_data = GRN[GRN['TF'] == TF]

            # Sort the data by importance score in descending order
            sorted_tf_data = tf_data.sort_values(by='importance', ascending=False)

            # Calculate the threshold percentile (default is 50th percentile)
            percentile_threshold = sorted_tf_data['importance'].quantile(quantile)

            # Select regulon genes above the percentile threshold
            regulon_genes = sorted_tf_data[sorted_tf_data['importance'] > percentile_threshold]['target'].tolist()

            # Determine significant targets by finding the intersection of DEGs and regulon genes
            SIGNIFICANT_TARGETS = list(set(DEGs) & set(regulon_genes))

            # Store results in the dictionary
            TF_targets_dict[celltype][TF] = {'SIGNIFICANT_TARGETS': SIGNIFICANT_TARGETS}
    
    return TF_targets_dict

import pandas as pd

def generate_gene_count_dataframe(tf_targets_dict):
    """
    Processes a dictionary containing transcription factors (TFs) and their significant targets to produce a DataFrame
    summarizing the count of significant targets for each TF in each cell type.

    Parameters:
    tf_targets_dict (dict): A dictionary where keys are cell types and values are dictionaries of TFs and their significant targets.

    Returns:
    DataFrame: A pandas DataFrame where rows correspond to cell types and columns to TFs, populated with the count of significant targets.
    """
    # Initialize an empty dictionary to store gene counts
    data = {}

    # Process each cell type and TF to calculate the count of significant targets
    for celltype, tfs in tf_targets_dict.items():
        data[celltype] = {}
        for tf, targets in tfs.items():
            target_count = len(targets['SIGNIFICANT_TARGETS'])
            data[celltype][tf] = target_count
            print(tf, 'SIGNIFICANT_TARGETS:', target_count)

    # Convert the dictionary to a DataFrame and fill missing values with 0
    gene_count_df = pd.DataFrame.from_dict(data, orient='index').fillna(0)
    gene_count_df = gene_count_df.astype(int)  # Ensure data type is integer

    return gene_count_df

# Example usage:
# Assuming you have a dictionary `TF_targets_dict` structured appropriately:
# df = generate_gene_count_dataframe(TF_targets_dict)
# print(df)



def calculate_gene_overlap(tf_targets_dict):
    """
    Calculates and returns a DataFrame representing the overlap of significant genes across different cell types based
    on transcription factor (TF) target data.

    Parameters:
    tf_targets_dict (dict): A dictionary where keys are cell types and values are dictionaries of TFs and their significant targets.

    Returns:
    DataFrame: A pandas DataFrame where rows and columns correspond to cell types, populated with the counts of overlapping significant genes.
    """
    # Initialize an empty dictionary to store all significant genes per cell type
    all_genes_per_celltype = {}

    # Compile all significant genes for each cell type into a set
    for celltype, tfs in tf_targets_dict.items():
        genes_set = set()
        for tf, info in tfs.items():
            genes_set.update(info['SIGNIFICANT_TARGETS'])
        all_genes_per_celltype[celltype] = genes_set

    # Initialize a DataFrame to store the overlap counts
    overlap_df = pd.DataFrame(index=all_genes_per_celltype.keys(), columns=all_genes_per_celltype.keys(), data=0)

    # Calculate the overlap for each pair of cell types
    for celltype1, genes1 in all_genes_per_celltype.items():
        for celltype2, genes2 in all_genes_per_celltype.items():
            overlap_count = len(genes1.intersection(genes2))
            overlap_df.loc[celltype1, celltype2] = overlap_count

    return overlap_df

# Example usage:
# Assuming you have a dictionary `TF_targets_dict` structured appropriately:
# overlap_df = calculate_gene_overlap(TF_targets_dict)
# print(overlap_df)


# Define a function to select regulon genes above the 50th percentile of importance scores
def significant_targets(GRN, TF_list):
    TF_targets_dict = {}
    for TF in TF_list:
        # Filter the data for the current TF
        tf_data = GRN[GRN['TF'] == TF]

        # Sort the data by importance score in descending order
        sorted_tf_data = tf_data.sort_values(by='importance', ascending=False)

        # Calculate the 50th percentile of importance scores
        percentile_50 = sorted_tf_data['importance'].quantile(quantile)

        # Select regulon genes above the 50th percentile
        regulon_genes = sorted_tf_data[sorted_tf_data['importance'] > percentile_50][['target']].target.tolist()

        UPREGULATED_TARGETS = list(set(UPREGULATED_GENES.ID.tolist()) & set(regulon_genes))
        DOWNREGULATED_TARGETS = list(set(DOWNREGULATED_GENES.ID.tolist()) & set(regulon_genes))
        SIGNIFICANT_TARGETS = list(set(SIGNIFICANT_GENES.ID.tolist()) & set(regulon_genes))
        print(TF, 'UPREGULATED_TARGETS:', len(UPREGULATED_TARGETS))
        print(TF, 'DOWNREGULATED_TARGETS:', len(DOWNREGULATED_TARGETS))
        print(TF, 'SIGNIFICANT_TARGETS:', len(SIGNIFICANT_TARGETS))

        TF_targets_dict[TF] = {'UPREGULATED_TARGETS': UPREGULATED_TARGETS, 'DOWNREGULATED_TARGETS': DOWNREGULATED_TARGETS, 'SIGNIFICANT_TARGETS': SIGNIFICANT_TARGETS}
    return TF_targets_dict

#TF_list = ['MITF', 'KLF12', 'GLIS3']
# TF_list= ['MITF']
#TF_targets_dict = significant_targets(SCENIC_GRN, TF_list)