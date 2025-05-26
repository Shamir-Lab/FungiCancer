import pandas as pd
import numpy as np
import random
from scipy.stats import linregress
from scipy import stats
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection
import math
from sklearn.linear_model import LogisticRegression
from statsmodels.formula.api import ols
import statsmodels.api as sm
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt
#from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
from scipy.stats import pearsonr
from openpyxl import load_workbook
from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl.styles import Font
def extract_kingdom_genus(col_name):
    kingdom = next((item.split('__')[1] for item in col_name.split('.') if item.startswith('k__')), None)
    genus = next((item.split('__')[1] for item in col_name.split('.') if item.startswith('g__')), None)
    return kingdom, genus
	
	

def sum_counts(df, count_columns, dataset_name):
    """
    Sum the counts for each column.

    Parameters:
    - df: The DataFrame containing the data for this dataset.
    - count_columns: List of column names that contain the count data.
    - dataset_name: String to label the resulting dataset.

    Returns:
    A Series with aggregated counts.
    """
    # Sum the counts for each column
    summed_counts = df[count_columns].sum()
    summed_counts.name = dataset_name
    return summed_counts

def sum_counts_by_cancer(df, count_columns, dataset_name, cancer_col='investigation'):
    """
    Sum the counts for each column grouped by cancer type.

    Parameters:
    - df: The DataFrame containing the data for this dataset.
    - count_columns: List of column names that contain the count data.
    - dataset_name: String to label the resulting dataset.
    - cancer_col: Name of the column that indicates the cancer type.

    Returns:
    A DataFrame with aggregated counts per cancer type.
    """
    summed_counts = df.groupby(cancer_col)[count_columns].sum()
    # Add a level of index for the dataset name
    summed_counts = summed_counts.set_index([pd.Series([dataset_name] * summed_counts.shape[0], name='Dataset')], append=True)
    return summed_counts



def calculate_relative_abundance(df, species_cols):
    # Calculate relative abundances
    rel_abundance = df[species_cols].div(df[species_cols].sum(axis=1), axis=0)
    # Keep all other columns
    for col in df.columns:
        if col not in species_cols:
            rel_abundance[col] = df[col]
    return rel_abundance




def compute_ratios_by_cancer_filtered(df, dataset_name, fungus_columns, side=True):
    """
    Compute the pairwise ratios of counts for each fungus pair and each cancer.
    """
    ratios = []
    for cancer_type in df.index.get_level_values('investigation').unique():
        # Select the data for this cancer type
        df_cancer = df.loc[(slice(None), cancer_type), :].droplevel('investigation').loc[dataset_name,]

        for i in range(len(fungus_columns)):
            for j in range(i + 1, len(fungus_columns)):
                if (side):
                    species1 = fungus_columns[i]
                    species2 = fungus_columns[j]
                else:
                    species1 = fungus_columns[j]
                    species2 = fungus_columns[i]
                if df_cancer[species1].sum() != 0 and df_cancer[species2].sum() != 0:
                    ratio = df_cancer[species1].sum() / df_cancer[species2].sum()
                    if not np.isnan(ratio) and not np.isinf(ratio):
                        ratios.append((cancer_type, f'{species1}/{species2}', ratio))
    return pd.DataFrame(ratios, columns=['Cancer', 'Pair', 'Ratio'])


def compute_ratios_by_cancer(df, dataset_name, fungus_columns):
    """
    Compute the pairwise ratios of counts for each fungus pair and each cancer.
    """
    ratios = []
    for cancer_type in df.index.get_level_values('investigation').unique():
        # Select the data for this cancer type
        df_cancer = df.loc[(slice(None), cancer_type), :].droplevel('investigation').loc[dataset_name,]

        for i in range(len(fungus_columns)):
            #for j in range(i + 1, len(fungus_columns)):
            for j in range(len(fungus_columns)):
                if i != j: #new row
                    species1 = fungus_columns[i]
                    species2 = fungus_columns[j]
                    if df_cancer[species1].sum() != 0 and df_cancer[species2].sum() != 0:
                        ratio = df_cancer[species1].sum() / df_cancer[species2].sum()
                        if not np.isnan(ratio) and not np.isinf(ratio):
                            ratios.append((cancer_type, f'{species1}/{species2}', ratio))
    return pd.DataFrame(ratios, columns=['Cancer', 'Pair', 'Ratio'])

def count_zeros_by_cancer(df, species_col, cancer_col='investigation'):
    """Return the number of zero counts for the given species column, grouped by cancer type."""
    return (df[df[species_col] == 0].groupby(cancer_col).size())

def compute_correlations_and_zero_counts(df1, df2, species_list, cancer_col='investigation'):
    correlations = pd.DataFrame(index=species_list)
    for species in species_list:
        correlations.at[species, 'corr_raw_count_all_cancer'] = df1[species].corr(df2[species])

        # Loop through cancer types
        cancer_types = df1[cancer_col].unique()
        for cancer_type in cancer_types:
            sub_df1 = df1[df1[cancer_col] == cancer_type]
            sub_df2 = df2[df2[cancer_col] == cancer_type]

            correlations.at[species, f'corr_raw_count_{cancer_type}'] = sub_df1[species].corr(sub_df2[species])

            # Zero count calculations
            zero_counts_df1 = count_zeros_by_cancer(sub_df1, species, cancer_col)
            zero_counts_df2 = count_zeros_by_cancer(sub_df2, species, cancer_col)

            correlations.at[species, f'zero_count_{cancer_type}_salzberg'] = zero_counts_df1.get(cancer_type, 0)
            correlations.at[species, f'zero_count_{cancer_type}_poore'] = zero_counts_df2.get(cancer_type, 0)

    return correlations

def generate_correlation_table_v2(df1, df2, species, cancer_col='investigation'):
    correlations = pd.DataFrame(index=species)

    # Overall correlations and zero counts
    for species_name in species:
        correlations.at[species_name, 'corr_raw_count_all_cancer'] = df1[species_name].corr(df2[species_name])
        correlations.at[species_name, 'zero_count_salzberg'] = (df1[species_name] == 0).sum()
        correlations.at[species_name, 'zero_count_poore'] = (df2[species_name] == 0).sum()

    # Per-cancer correlations and zero counts
    cancer_types = df1[cancer_col].unique()
    for cancer_type in cancer_types:
        df1_cancer = df1[df1[cancer_col] == cancer_type]
        df2_cancer = df2[df2[cancer_col] == cancer_type]
        for species_name in species:
            correlations.at[species_name, f'corr_raw_count_{cancer_type}'] = df1_cancer[species_name].corr(df2_cancer[species_name])
            correlations.at[species_name, f'zero_count_{cancer_type}_salzberg'] = (df1_cancer[species_name] == 0).sum()
            correlations.at[species_name, f'zero_count_{cancer_type}_poore'] = (df2_cancer[species_name] == 0).sum()
            correlations.at[species_name, f'zero_count_both_{cancer_type}'] = ((df1_cancer[species_name] == 0) & (df2_cancer[species_name] == 0)).sum()

            # If you want to include relative abundances
            df1_ra = calculate_relative_abundance(df1_cancer, [species_name])
            df2_ra = calculate_relative_abundance(df2_cancer, [species_name])
            correlations.at[species_name, f'corr_RA_{cancer_type}'] = df1_ra[species_name].corr(df2_ra[species_name])

    return correlations


def generate_correlation_table_sorted(df1, df2, species, cancer_col='investigation'):
    correlations = pd.DataFrame(index=species)

    # Overall correlations and zero counts
    for species_name in species:
        correlations.at[species_name, 'corr_raw_count_all_cancer'] = df1[species_name].corr(df2[species_name])
        correlations.at[species_name, 'zero_count_salzberg'] = (df1[species_name] == 0).sum()
        correlations.at[species_name, 'zero_count_poore'] = (df2[species_name] == 0).sum()
        correlations.at[species_name, 'non_zero_samples_both_datasets'] = ((df1[species_name] != 0) & (df2[species_name] != 0)).sum()

    # Per-cancer correlations and zero counts
    cancer_types = df1[cancer_col].unique()
    for cancer_type in cancer_types:
        df1_cancer = df1[df1[cancer_col] == cancer_type]
        df2_cancer = df2[df2[cancer_col] == cancer_type]
        for species_name in species:
            correlations.at[species_name, f'corr_raw_count_{cancer_type}'] = df1_cancer[species_name].corr(df2_cancer[species_name])
            correlations.at[species_name, f'zero_count_{cancer_type}_salzberg'] = (df1_cancer[species_name] == 0).sum()
            correlations.at[species_name, f'zero_count_{cancer_type}_poore'] = (df2_cancer[species_name] == 0).sum()
            correlations.at[species_name, f'zero_count_both_{cancer_type}'] = ((df1_cancer[species_name] == 0) & (df2_cancer[species_name] == 0)).sum()
            correlations.at[species_name, f'non_zero_samples_{cancer_type}_both_datasets'] = ((df1_cancer[species_name] != 0) & (df2_cancer[species_name] != 0)).sum()
            # If you want to include relative abundances
            df1_ra = calculate_relative_abundance(df1_cancer, [species_name])
            df2_ra = calculate_relative_abundance(df2_cancer, [species_name])
            correlations.at[species_name, f'corr_RA_{cancer_type}'] = df1_ra[species_name].corr(df2_ra[species_name])

    mean_abundances = df1[species].mean()
    correlations_sorted = correlations.loc[mean_abundances.sort_values(ascending=False).index]
    return correlations_sorted

def pivot_fungi_counts(df_res, p_val_threshold):
    significant_fungi = df_res[df_res['p_value_fdr'] < p_val_threshold]

    grouped = significant_fungi.groupby(['fungus_id', 'cancer_investigation']).size().reset_index(name='counts')

    # Pivot the data to get counts of unique analyses in which each fungi was significant
    pivot_table = grouped.pivot_table(index='cancer_investigation', columns='counts', aggfunc='size', fill_value=0)
    return pivot_table

def map_genome_to_species_fungusid(df, taxonomy_df):
    # Create a mapping from genomeID to Species name, removing the "s__" prefix
    mapping = taxonomy_df.set_index('genomeID')['Species'].str.split("__").str[-1].to_dict()
    
    # Use the mapping to replace genomeID with Species name
    df['fungus_id'] = df['fungus_id'].map(mapping)
    
    return df

def extract_significant_taxa(df_res, dataframes, threshold, p_val_threshold):
    # Filter significant fungi
    significant_fungi = df_res[df_res['p_value_fdr'] < p_val_threshold]
    
    # Group by 'fungus_id' and 'cancer_investigation', then count
    grouped = significant_fungi.groupby(['fungus_id', 'cancer_investigation']).size().reset_index(name='counts')
    
    # Filter those that have counts greater than or equal to the threshold
    filtered_group = grouped[grouped['counts'] >= threshold]
    result_df = map_genome_to_species_fungusid(filtered_group, dataframes["poore_2022_fungi_species"])
    
    return result_df[['cancer_investigation', 'fungus_id']]

def compute_ratios(df, dataset_a, dataset_b):
    """
    Compute the pairwise ratios of counts for each fungus pair.
    """

    ratios = []
    for i, species1 in enumerate(df.index):
        for j, species2 in enumerate(df.index):
            if i != j:
                ratio_poore = df.loc[species1, dataset_a] / df.loc[species2, dataset_a]
                ratio_salz = df.loc[species1, dataset_b] / df.loc[species2, dataset_b]

                if not (np.isnan(ratio_poore) or np.isinf(ratio_poore)) and not (np.isnan(ratio_salz) or np.isinf(ratio_salz)):
                    ratios.append((f'{species1}/{species2}', ratio_poore, ratio_salz))

    ratios_df = pd.DataFrame(ratios, columns=['Pair', 'Ratio_'+dataset_a, 'Ratio_'+dataset_b])

    ratios_df_no_zeros = ratios_df[(ratios_df['Ratio_'+dataset_a] != 0) & (ratios_df['Ratio_'+dataset_b] != 0)]

    print(f'list of species that were removed because their total count was zero: {list(set(ratios_df.Pair).difference(set(ratios_df_no_zeros.Pair)))}')

    return ratios_df_no_zeros

def compute_ratios_filtered(df, dataset_a, dataset_b, side_type=True):
    """
    Compute the pairwise ratios of counts for each fungus pair.
    """
    ratios = []
    for i in range(len(df.index)):
      for j in range(i + 1, len(df.index)):
        if (side_type):
          species1 = df.index[i]
          species2 = df.index[j]
        else:
          species1 = df.index[j]
          species2 = df.index[i]

        ratio_poore = df.loc[species1, dataset_a] / df.loc[species2, dataset_a]
        ratio_salz = df.loc[species1, dataset_b] / df.loc[species2, dataset_b]

        if not (np.isnan(ratio_poore) or np.isinf(ratio_poore)) and not (np.isnan(ratio_salz) or np.isinf(ratio_salz)):
            ratios.append((f'{species1}/{species2}', ratio_poore, ratio_salz))

    ratios_df = pd.DataFrame(ratios, columns=['Pair', 'Ratio_'+dataset_a, 'Ratio_'+dataset_b])

    ratios_df_no_zeros = ratios_df[(ratios_df['Ratio_'+dataset_a] != 0) & (ratios_df['Ratio_'+dataset_b] != 0)]

    print(f'list of species that were removed because their total count was zero: {list(set(ratios_df.Pair).difference(set(ratios_df_no_zeros.Pair)))}')

    return ratios_df_no_zeros

def calc_indepedent_correlations(resulting_df_allcancers,ratios_df_1, dataset_1_a, dataset_1_b,ratios_df_2,dataset_2_a, dataset_2_b, n_iterations):
  corr_list_1 = []
  corr_list_2 = []

  for iter in range(n_iterations):

    n_samples = 154
    #create a permutation of the species types
    rand_permutation = random.sample(sorted(resulting_df_allcancers.index),n_samples)

    #create a list of ratios
    new_list = []
    for i in range(int(n_samples/2)):
      element = ''.join([rand_permutation[i],'/', rand_permutation[i+1]])
      new_list = new_list +[element]

    # take only samoples ratios
    merged_ratios_1 = ratios_df_1[ratios_df_1.Pair.isin(new_list)]
    merged_ratios_2 = ratios_df_2[ratios_df_2.Pair.isin(new_list)]

    # calc pearson R
    slope, intercept, r_value_1,  p_value, std_err = linregress(merged_ratios_1['Ratio_'+dataset_1_a], merged_ratios_1['Ratio_'+dataset_1_b])
    slope, intercept, r_value_2,  p_value, std_err = linregress(merged_ratios_2['Ratio_'+dataset_2_a], merged_ratios_2['Ratio_'+dataset_2_b])

    corr_list_1 = corr_list_1+[r_value_1]
    corr_list_2 = corr_list_2+[r_value_2]

  return(corr_list_1,corr_list_2)


def generate_correlation_table_sorted_new(df1, df2, species, cancer_col='investigation'):
    correlations = pd.DataFrame(index=species)

    # Overall correlations and zero counts
    for species_name in species:
        correlations.at[species_name, 'corr_raw_count_all_cancer'] = df1[species_name].corr(df2[species_name])
        correlations.at[species_name, 'zero_count_salzberg'] = (df1[species_name] == 0).sum()
        correlations.at[species_name, 'zero_count_poore'] = (df2[species_name] == 0).sum()
        non_zero_vec = (df1[species_name] != 0) & (df2[species_name] != 0)
        correlations.at[species_name, 'non_zero_samples_both_datasets'] = (non_zero_vec).sum()
        correlations.at[species_name, 'corr_non_zero_raw_count_all_cancer'] = df1[non_zero_vec][species_name].corr(df2[non_zero_vec][species_name])

    # Per-cancer correlations and zero counts
    cancer_types = df1[cancer_col].unique()
    for cancer_type in cancer_types:
        df1_cancer = df1[df1[cancer_col] == cancer_type]
        df2_cancer = df2[df2[cancer_col] == cancer_type]
        for species_name in species:
            correlations.at[species_name, f'corr_raw_count_{cancer_type}'] = df1_cancer[species_name].corr(df2_cancer[species_name])
            correlations.at[species_name, f'zero_count_salzberg_{cancer_type}'] = (df1_cancer[species_name] == 0).sum()
            correlations.at[species_name, f'zero_count_poore_{cancer_type}'] = (df2_cancer[species_name] == 0).sum()
            non_zero_vec_cacner  = (df1_cancer[species_name] != 0) & (df2_cancer[species_name] != 0)
            correlations.at[species_name, f'non_zero_samples_both_datasets_{cancer_type}'] = (non_zero_vec_cacner).sum()
            correlations.at[species_name, f'corr_non_zero_raw_count_all_cancer_{cancer_type}'] = df1_cancer[non_zero_vec_cacner][species_name].corr(df2_cancer[non_zero_vec_cacner][species_name])

    mean_abundances = df1[species].mean()
    correlations_sorted = correlations.loc[mean_abundances.sort_values(ascending=False).index]
    return correlations_sorted

def generate_full_species_corr_table(salzberg_merged_samesamples_sametaxa,poore_2022_wisoverlap_merged_sampesamples_sametaxa,poore_2020_merged_sampesamples_sametaxa,all_together_same_taxa):
  # revise col names Poore -> G23, NH22, P20 and remove 'TCGA-' from headers

  c2_sorted = generate_correlation_table_sorted_new(salzberg_merged_samesamples_sametaxa,poore_2022_wisoverlap_merged_sampesamples_sametaxa, all_together_same_taxa)
  c2_sorted.columns = c2_sorted.columns.str.replace('poore','NH22')
  c2_sorted.columns = c2_sorted.columns.str.replace('salzberg','G23')
  c2_sorted.columns = c2_sorted.columns.str.replace('raw_count_all_cancer','NH22_G23')
  c2_sorted.columns = c2_sorted.columns.str.replace('both_datasets','NH22_G23')
  c2_sorted.columns = c2_sorted.columns.str.replace('_TCGA-','_')
  c2_sorted.columns = c2_sorted.columns.str.replace('corr_raw_count_','corr_NH22_G23_')
  c2_sorted = c2_sorted.sort_index()

  # extract cols of P20 and merge and remove 'TCGA-' from headers
  c3_sorted = generate_correlation_table_sorted_new(salzberg_merged_samesamples_sametaxa,poore_2020_merged_sampesamples_sametaxa, all_together_same_taxa)
  c3_sorted.columns = c3_sorted.columns.str.replace('poore','P20')
  c3_sorted.columns = c3_sorted.columns.str.replace('raw_count_all_cancer','P20_G23')
  c3_sorted.columns = c3_sorted.columns.str.replace('both_datasets','P20_G23')
  c3_sorted.columns = c3_sorted.columns.str.replace('_TCGA-','_')
  c3_sorted.columns = c3_sorted.columns.str.replace('corr_raw_count_','corr_P20_G23_')

  c3_sorted = c3_sorted.sort_index()

  c3_sorted_partial = c3_sorted[['corr_P20_G23','zero_count_P20','non_zero_samples_P20_G23','corr_non_zero_P20_G23','corr_P20_G23_BRCA','zero_count_P20_BRCA','non_zero_samples_P20_G23_BRCA','corr_non_zero_P20_G23_BRCA',
                                'corr_P20_G23_HNSC','zero_count_P20_HNSC','non_zero_samples_P20_G23_HNSC','corr_non_zero_P20_G23_HNSC','corr_P20_G23_BLCA','zero_count_P20_BLCA','non_zero_samples_P20_G23_BLCA','corr_non_zero_P20_G23_BLCA']]

  c_merged = c3_sorted_partial.join(c2_sorted)
  all_cancer_cols = ['corr_P20_G23','corr_NH22_G23', 'zero_count_P20', 'zero_count_NH22', 'zero_count_G23','non_zero_samples_P20_G23','corr_non_zero_P20_G23','non_zero_samples_NH22_G23','corr_non_zero_NH22_G23']
  per_cancer_cols = list(itertools.chain(*[[f'corr_P20_G23_{cancer_type}', f'corr_NH22_G23_{cancer_type}',f'zero_count_P20_{cancer_type}',f'zero_count_NH22_{cancer_type}', f'zero_count_G23_{cancer_type}', f'non_zero_samples_P20_G23_{cancer_type}',f'corr_non_zero_P20_G23_{cancer_type}',f'non_zero_samples_NH22_G23_{cancer_type}',f'corr_non_zero_NH22_G23_{cancer_type}'] for cancer_type in ['BRCA','BLCA','HNSC']]))
  c_merged = c_merged[all_cancer_cols+per_cancer_cols]
  return(c_merged)

#Function for voom-snm
def extract_ppv_per_threshold(df_nromalized_with_cancer,cancer_name,only_rna_seq,pos_class,genome_id,threshold):
  x = df_nromalized_with_cancer[(df_nromalized_with_cancer.investigation == cancer_name)][[genome_id.values[0],'experimental_strategy','sample_type']]
  x = x[x.sample_type.isin(['Primary Tumor','Solid Tissue Normal'])]

  x[genome_id.values[0]] = round(x[genome_id.values[0]],4)

  if(only_rna_seq):
    x = x[x.experimental_strategy.isin(experimental_strategy_name)]
  print('Total Samples=', len(x))
  conteg_df = pd.crosstab(x[genome_id.values[0]]>=threshold, x.sample_type==pos_class)

  conteg_df = conteg_df.add_prefix('tot ').reset_index().rename_axis(None, axis=1)
  print('PPV=',round((conteg_df['tot True'][1] / conteg_df[conteg_df[genome_id.values[0]] == True][['tot False','tot True']].sum(axis=1)[1])*100,2),'%')
  print('Sensitivity=', round((conteg_df['tot True'][1]/(conteg_df['tot True'][1]+conteg_df['tot True'][0]))*100,2),'%')
  print('Error Rate (FP+FN)=', round((conteg_df['tot False'][1]+conteg_df['tot True'][0])/len(x)*100,2),'%')

  return(conteg_df)

def calculate_p_values_per_factor_total_counts(df_raw,meta_data,res_path,comparisons_dict,raw_fungi_for_analysis):
  #calculate total redad counts per sample
  df_raw['total_counts'] = df_raw.iloc[0:df_raw.shape[0],1:225].sum(axis=1)

  #Take only samples of primary tumor and RNA-Sew, and in our cohort
  meta_data = meta_data[meta_data.sample_type == 'Primary Tumor']
  meta_data = meta_data[meta_data.experimental_strategy == 'RNA-Seq']
  meta_data = meta_data[meta_data.sampleid.isin(raw_fungi_for_analysis.Sample)]

  #crate empty dataframe
  df_total_counts_per_factor = pd.DataFrame(columns=['cancer_investigation']+list(comparisons_dict.keys()))
  df_total_counts_per_factor['cancer_investigation'] = meta_data['investigation'].unique()

  #for factor_type in ['race', 'gender', 'over_65', 'obese']:
  for ind in range(len(comparisons_dict.keys())):
    variable = list(comparisons_dict.keys())[ind]
    factor = comparisons_dict[variable][4]

    op1 = comparisons_dict[variable][0]
    op2 = comparisons_dict[variable][1]
    confounders = comparisons_dict[variable][2]
    need_dummies = comparisons_dict[variable][3]

    df_raw_with_race = df_raw.merge(meta_data[['sampleid','investigation','sample_type','experimental_strategy',factor]], on='sampleid')

    #drop rows with no value
    df_raw_with_race = df_raw_with_race[~df_raw_with_race[factor].isna()]

    if (factor =='race'):
      df_raw_with_race = df_raw_with_race[df_raw_with_race[factor].isin(['WHITE','BLACK OR AFRICAN AMERICAN','ASIAN'])]

    for cancer_type in meta_data['investigation'].unique():
      df_raw_with_race_cancer = df_raw_with_race[df_raw_with_race.investigation == cancer_type]

      #statistical test
      if (len(df_raw_with_race_cancer[factor].value_counts()) > 1):
        vec_a = df_raw_with_race_cancer[df_raw_with_race_cancer[factor] ==op1]['total_counts']
        vec_b = df_raw_with_race_cancer[df_raw_with_race_cancer[factor] == op2]['total_counts']
        #U1, p_val = mannwhitneyu(vec_a,vec_b, method="auto")
        U1, p_val = stats.ttest_ind(vec_a,vec_b)

      else:
        p_val = np.NaN

      df_total_counts_per_factor.loc[(df_total_counts_per_factor.cancer_investigation == cancer_type), variable] = round(p_val, 4)

  res = (df_total_counts_per_factor[['cancer_investigation']+list(comparisons_dict.keys())])
  res.to_csv(f'{res_path}/sup_table_3_p_values_total_counts.csv')
  return (res)

def calculate_p_values_per_factor_total_counts_anova(df_raw,meta_data,res_path,raw_fungi_for_analysis):
  #calculate total redad counts per sample
  df_raw['total_counts'] = df_raw.iloc[0:df_raw.shape[0],1:225].sum(axis=1)

  #Take only samples of primary tumor and RNA-Sew, and in our cohort
  meta_data = meta_data[meta_data.sample_type == 'Primary Tumor']
  meta_data = meta_data[meta_data.experimental_strategy == 'RNA-Seq']
  meta_data = meta_data[meta_data.sampleid.isin(raw_fungi_for_analysis.Sample)]

  #crate empty dataframe
  df_total_counts_per_factor = pd.DataFrame(columns=['cancer_investigation']+['race', 'gender', 'over_65', 'obese'])
  df_total_counts_per_factor['cancer_investigation'] = meta_data['investigation'].unique()

  for factor_type in ['race', 'gender', 'over_65', 'obese']:
    df_raw_with_race = df_raw.merge(meta_data[['sampleid','investigation','sample_type','experimental_strategy',factor_type]], on='sampleid')

    #drop rows with no value
    df_raw_with_race = df_raw_with_race[~df_raw_with_race[factor_type].isna()]

    if (factor_type =='race'):
      df_raw_with_race = df_raw_with_race[df_raw_with_race[factor_type].isin(['WHITE','BLACK OR AFRICAN AMERICAN','ASIAN'])]

    for cancer_type in meta_data['investigation'].unique():
      df_raw_with_race_cancer = df_raw_with_race[df_raw_with_race.investigation == cancer_type]

      #statistical test
      if (len(df_raw_with_race_cancer[factor_type].value_counts()) > 1):
        if (factor_type=='race'):
          white = df_raw_with_race_cancer[df_raw_with_race_cancer.race == 'WHITE']['total_counts']
          black = df_raw_with_race_cancer[df_raw_with_race_cancer.race == 'BLACK OR AFRICAN AMERICAN']['total_counts']
          asian = df_raw_with_race_cancer[df_raw_with_race_cancer.race == 'ASIAN']['total_counts']
          p_val = stats.f_oneway(white,asian,black)[1]
        else:
          vec_a = df_raw_with_race_cancer[df_raw_with_race_cancer[factor_type] == df_raw_with_race_cancer[factor_type].value_counts().index[0]]['total_counts']
          vec_b = df_raw_with_race_cancer[df_raw_with_race_cancer[factor_type] == df_raw_with_race_cancer[factor_type].value_counts().index[1]]['total_counts']
          #U1, p_val = mannwhitneyu(vec_a,vec_b, method="auto")
          U1, p_val = stats.ttest_ind(vec_a,vec_b)

      else:
        p_val = np.NaN

      df_total_counts_per_factor.loc[(df_total_counts_per_factor.cancer_investigation == cancer_type), factor_type] = round(p_val, 4)

  res = df_total_counts_per_factor[['cancer_investigation','race', 'gender', 'over_65', 'obese']]
  res.to_csv(f'{res_path}/sup_table_3_p_values_total_counts_anova.csv')
  return(res)



def add_bmi_and_age_groups(merged):
    bins = [0, 18.5, 25, 30, 35, 40, float('inf')]
    labels = ['<18.5', '18.5-25', '25-30', '30-35', '35-40', '40+']
    merged['BMI_Group'] = pd.cut(merged['BMI'], bins=bins, labels=labels, right=False)
    merged['BMI_Group'].value_counts()

    bins2 = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    labels2 = ['0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90', '90-100']
    merged['age_at_diagnosis_group'] = pd.cut(merged['age_at_diagnosis'], bins=bins2, labels=labels2, right=False)

    return merged

def format_with_percentage(count_df):
    total_per_cancer = count_df.sum(axis=1)
    percentage_df = (count_df.div(total_per_cancer, axis=0) * 100).round(2)
    formatted_df = count_df.astype(str) + " (" + percentage_df.astype(str) + "%)"
    return formatted_df

def get_top_n_categories_for_each_cancer(df, column, n=5):
    all_cancers = df['investigation'].unique()
    top_categories_df = pd.DataFrame(index=[f"{i+1}st" if i == 0 else f"{i+1}nd" if i == 1 else f"{i+1}rd" if i == 2 else f"{i+1}th" for i in range(n)], columns=all_cancers)

    for cancer in all_cancers:
        top_counts = df[df['investigation'] == cancer][column].value_counts().head(n)
        total_counts = df[df['investigation'] == cancer][column].count()
        formatted_categories = top_counts.index + " (" + top_counts.astype(str) + ", " + (100 * top_counts / total_counts).round(2).astype(str) + "%)"
        top_categories_df[cancer][:len(formatted_categories)] = formatted_categories.values

    return top_categories_df

replacements = {
    "gender": "No gender",
    "age_at_diagnosis_group": "No age group",
    "race": "No race",
    "stage_numbered": "No stage",
    "BMI_Group": "No BMI Group",
    "data_submitting_center_label": "No Data Submitting Center",
    "tissue_source_site_label": "No Tissue Source Site Label"
}

def add_category_if_categorical(df, column, category):
    if pd.api.types.is_categorical_dtype(df[column]):
        if category not in df[column].cat.categories:
            df[column] = df[column].cat.add_categories([category])
    return df[column]

def create_conf_and_vars_xlsx(merged,res_path):
    """
    receives a df with the metadata and fungi counts
    saves a table with the counts for each confounder and variable of interest to the res_path
    """
    merged['investigation'] = add_category_if_categorical(merged, 'investigation', 'No Investigation')

    for column, replacement in replacements.items():
        merged[column] = add_category_if_categorical(merged, column, replacement)

    merged.fillna(replacements, inplace=True)

    vars_and_conf = ["gender", "age_at_diagnosis_group", "race", "stage_numbered", "BMI_Group", "data_submitting_center_label"]

    pivoted_counts = {}
    for variable in vars_and_conf:
        counts = merged.groupby(["investigation", variable]).size().unstack(fill_value=0)
        formatted_counts = format_with_percentage(counts)
        pivoted_counts[variable] = formatted_counts.T

    for column in ['histological_diagnosis_label', 'tissue_source_site_label']:
        top_n_df = get_top_n_categories_for_each_cancer(merged, column, n=5)
        ordered_columns = pivoted_counts[vars_and_conf[0]].columns
        top_n_df = top_n_df.reindex(columns=ordered_columns)
        pivoted_counts[column] = top_n_df

    with pd.ExcelWriter("all_conf_and_vars_counts.xlsx", mode="w", engine="openpyxl") as writer:
        start_row = 0
        for variable, counts in pivoted_counts.items():
            counts.to_excel(writer, sheet_name="Combined", startrow=start_row, startcol=0)
            start_row += counts.shape[0] + 3

        workbook = writer.book
        sheet = workbook['Combined']

        sheet['A50'] = 'histological_diagnosis_label'
        sheet['A50'].font = Font(bold=True)

        sheet['A58'] = 'tissue_source_site_label'
        sheet['A58'].font = Font(bold=True)

        for col in sheet.columns:
            max_length = 0
            column = col[0].column_letter
            for cell in col:
                try:
                    if len(str(cell.value)) > max_length:
                        max_length = len(cell.value)
                except:
                    pass
            adjusted_width = (max_length + 2)
            sheet.column_dimensions[column].width = adjusted_width

        for row in sheet.iter_rows():
            max_height = 0
            for cell in row:
                if cell.value:
                    max_height = max(max_height, len(str(cell.value)) // adjusted_width)
            sheet.row_dimensions[row[0].row].height = max_height * 15

    workbook.save(f"{res_path}/sup_table_1.xlsx")

