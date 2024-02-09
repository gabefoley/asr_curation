import pandas as pd
import numpy as np
import ast
from Bio import SeqIO
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import random
from scipy.stats import fisher_exact

dbscan_params_eps = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 2, 5, 10, 20]
dbscan_params_min_samples = [1,2,3,4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 100, 200, 500, 1000, 5000]




def generate_dbscan_embeddings(df, embedding_col,
                               dbscan_params_eps = dbscan_params_eps,
                               dbscan_params_min_samples = dbscan_params_min_samples
                               ):
    df = df.fillna('None')

    if embedding_col not in df:
        print (f"Requested embedding column {embedding_col} is not in dataframe")
        return

    # df['embeddings_encoded'] = df['Prot_T5 Embed Encoded'].apply(lambda s: ast.literal_eval(s))
    #
    #
    # embeddings = np.vstack(df['embeddings_encoded'].to_numpy())

    # df['embeddings_encoded'] = df['Prot_T5 Embed Encoded'].apply(lambda s: ast.literal_eval(s))


    embeddings = np.vstack(df[embedding_col].to_numpy())


    for eps_param in dbscan_params_eps:
        for min_samples_param in dbscan_params_min_samples:
            dbscan = DBSCAN(eps=eps_param, min_samples=min_samples_param)
            df[f'dbscan_eps={eps_param}_minsamples={min_samples_param}'] = dbscan.fit_predict(embeddings)

    return df




def compare_categorical_distributions(df1, df2, column_name, alpha=0.05):

    df1_counts = df1[column_name].value_counts().sort_index()
    df2_counts = df2[column_name].value_counts().sort_index()

    contingency_table = pd.concat([df1_counts, df2_counts], axis=1, keys=['df1', 'df2']).fillna(0).astype(int)

    if contingency_table.shape != (2, 2):
        return

    odds_ratio, p_value = fisher_exact(contingency_table)

    if p_value < alpha:
        print(contingency_table)
        print("There is a statistically significant difference between the distributions.")
        return column_name

    else:
        return

def generate_parameter_subsets(df, cols_to_check):
    '''
    Generate a dictionary with each parameter combination, and the list of annotation columns implied as to be excluded on'''
    parameter_subsets = defaultdict(list)

    for eps_param in dbscan_params_eps:
        for min_samples_param in dbscan_params_min_samples:
            noise_indices = set(df[df[f'dbscan_eps={eps_param}_minsamples={min_samples_param}'] == -1].index)
            noise_df = df.loc[list(noise_indices)]

            # Merge df1 and df2 with an indicator column
            merged = pd.merge(df, noise_df, on=['info'] + cols_to_check, how='left', indicator=True,
                              suffixes=('', '_df2'))

            # Filter rows where the indicator column is 'left_only'
            subset_wo_noise = merged[merged['_merge'] == 'left_only']

            # Drop the indicator column and reset the index of the result DataFrame
            subset_wo_noise = subset_wo_noise.drop(columns=['_merge']).reset_index(drop=True)
            subset_wo_noise = subset_wo_noise[[col for col in subset_wo_noise.columns if not col.endswith('_df2')]]

            diff_columns = []

            if subset_wo_noise.shape[0] != 0:

                for column in cols_to_check:
                    alpha = 0.05  # Significance level

                    diff_column = compare_categorical_distributions(df, noise_df, column, alpha)

                    if diff_column:
                        diff_columns.append(diff_column)
                parameter_subsets[f'dbscan_eps={eps_param}_minsamples={min_samples_param}'] = diff_columns

    return parameter_subsets

def generate_db_scan_coverage_plot(df, parameter_subsets, outpath):
    # Create a scatter plot
    plt.figure(figsize=(20, 8))

    # Create a dictionary to store the color mapping based on feature sets
    colour_mapping = {}
    feature_count = defaultdict(int)

    # Create a grid of points
    for i, eps_param in enumerate(dbscan_params_eps):
        for j, min_samples_param in enumerate(dbscan_params_min_samples):
            combination_name = f'dbscan_eps={eps_param}_minsamples={min_samples_param}'
            features = set(parameter_subsets.get(combination_name, []))  # Get the set of features for the combination

            # Convert the set to a frozenset to make it hashable
            features_frozen = frozenset(features)

            # Determine the color based on the feature set
            if features_frozen in colour_mapping:
                colour = colour_mapping[features_frozen]
                feature_count[features_frozen] += 1
            else:
                colour = plt.cm.tab20(len(colour_mapping))  # Use a different color for each unique feature set
                colour_mapping[features_frozen] = colour
                feature_count[features_frozen] += 1

            # Calculate the position on the grid
            x = i + 1  # Adjusted by 1 to start from 1
            y = j + 1  # Adjusted by 1 to start from 1

            # Plot a point at the (x, y) coordinates with the determined color
            plt.scatter(x, y, c=[colour], s=500, label=None)  # No label for points

            num_outliers = df[combination_name].tolist().count(-1)
            plt.text(x, y, str(num_outliers), fontsize=8, color='white', ha='center', va='center_baseline')

    # Sort the dictionary by values in ascending order
    sorted_features = dict(sorted(feature_count.items(), key=lambda item: -item[1]))
    legend_labels = []

    # Iterate through the sorted dictionary
    for features in sorted_features.keys():
        colour = colour_mapping[features]
        # Create legend entries for unique feature sets

        legend_labels.append(
            plt.Line2D([0], [0], marker='o', color='w', label=[x for x in list(features)], markersize=10,
                       markerfacecolor=colour))

    # Label the axes with values while keeping the distances/placement the same
    plt.xticks(range(1, len(dbscan_params_eps) + 1), dbscan_params_eps)
    plt.yticks(range(1, len(dbscan_params_min_samples) + 1), dbscan_params_min_samples)

    # Customize the plot
    plt.xlabel('DBSCAN Parameter - eps')
    plt.ylabel('DBSCAN Parameter - min_samples')
    plt.title('Combinations of DBSCAN Parameters')
    # plt.grid(True)
    # plt.legend(handles=legend_labels, title='Unique Feature Sets', bbox_to_anchor=(-10, -1, 5, 4))

    legend = plt.legend(handles=legend_labels, title='Unique Feature Sets', loc='upper left', bbox_to_anchor=(1, 1),
                        ncol=1)
    plt.subplots_adjust(right=0.7)  # Adjust the right margin to make space for the legend

    plt.savefig(f'{outpath}'.png)


def generate_dbscan_coverage(df, embedding_col, outpath, skip_cols=[]):
    embedding_cols = ['Prot_T5 Embed Encoded', 'Prot_T5 Embed Mean', 'Prot_T5 Embed CLS']
    df = generate_dbscan_embeddings(df, embedding_col)

    # We always want to skip these columns as they will be uninformative
    skip_cols += ['info', 'truncated_info', 'extracted_id',
                 'extracted_name', 'sequence_version', 'protein_existence',
                 'reviewed', 'UniProt_DB', 'cdd_evalue', 'cdd_incomplete', 'embeddings',
                 'embeddings2', 'IQR_marked', 'curated', 'Routine_exclusion', 'cdd_pssm_id',
                 'cdd_from', 'cdd_to', 'cdd_bitscore', 'lineage_thermo', 'organism_name', 'organism_id', 'length',
                 'mass', 'sequence_version', 'version', 'Length_2', 'annotation_score',
                 'lineage_subkingdom', 'volker_sadht', 'lineage_thermo', 'thermo_bacteria', 'lineage_thermo_no_metha',
                 'model_name', 'Prot_T5 Embed Encoded', 'Prot_T5 Embed Mean', 'Prot_T5 Embed CLS', 'embeddings_encoded',
                  'BRENDA', 'dbscan']

    cols_to_check = [x for x in df.keys() if x not in skip_cols and not any(x.startswith(term) for term in skip_cols)]



    parameter_subsets = generate_parameter_subsets(df, cols_to_check)

    # Write out the parameter subsets
    with open(f'{outpath}'.txt, 'w') as file:
        for item in parameter_subsets:
            file.write("%s\n" % item)

    generate_db_scan_coverage_plot(df, parameter_subsets, outpath)

