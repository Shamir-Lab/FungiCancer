import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import linregress,mannwhitneyu
import math
from sklearn.linear_model import LogisticRegression
from statsmodels.formula.api import ols
import statsmodels.api as sm
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles
from scipy.stats import pearsonr

def single_barplot_for_PSM(results, all_cancers_color_dict, title):
    try:
        counts_PSM = results[results['p_value_fdr'] < 0.05].groupby('cancer_investigation')['fungus_id'].count().reset_index(name='counts_PSM')

        if counts_PSM.empty:
            raise ValueError("min() arg is an empty sequence")

        plt.figure(figsize=(12, 6))
        ax = sns.barplot(x='cancer_investigation', y='counts_PSM', data=counts_PSM, palette=all_cancers_color_dict)

        for i, bar in enumerate(ax.patches):
            bar.set_facecolor(all_cancers_color_dict[counts_PSM.loc[i, 'cancer_investigation']])
            bar.set_alpha(1.0)  # Set full opacity
            bar.set_edgecolor("black")  # Set outline color
            bar.set_linewidth(1)  # Set outline width

        plt.xticks(rotation=90, fontsize=9)
        plt.xlabel('Cancer Type')
        plt.ylabel('Count of Significant p-values')
        plt.title(title)

        plt.show()
    except ValueError as e:
        if "min() arg is an empty sequence" in str(e):
            print("No data to plot")
        else:
            raise  # raise the error if it's not the one we expect


def paired_three_barplot_for_PSM(results_norm, results_not_norm, results_no_logcpm, all_cancers_color_dict, title):

    # Extracting counts for each investigation and dataset
    counts_PSM_norm = results_norm[results_norm['p_value_fdr'] < 0.05].groupby('cancer_investigation')['fungus_id'].count().reset_index(name='counts_PSM_norm')
    counts_PSM_not_norm = results_not_norm[results_not_norm['p_value_fdr'] < 0.05].groupby('cancer_investigation')['fungus_id'].count().reset_index(name='counts_PSM_not_norm')
    counts_PSM_no_logcpm = results_no_logcpm[results_no_logcpm['p_value_fdr'] < 0.05].groupby('cancer_investigation')['fungus_id'].count().reset_index(name='counts_PSM_no_logcpm')

    # Merging the counts
    merged_counts = counts_PSM_norm.merge(counts_PSM_not_norm, on='cancer_investigation', how='outer').merge(counts_PSM_no_logcpm, on='cancer_investigation', how='outer').fillna(0)

    # Setting up the plot
    bar_width = 0.25
    r1 = np.arange(len(merged_counts))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]

    plt.figure(figsize=(12, 6))

    # Plotting bars
    for idx, row in merged_counts.iterrows():
        plt.bar(r1[idx], row['counts_PSM_norm'], color=all_cancers_color_dict[row['cancer_investigation']], width=bar_width, edgecolor='grey', label='logCPM+Norm' if idx == 0 else "")
        plt.bar(r2[idx], row['counts_PSM_not_norm'], color=all_cancers_color_dict[row['cancer_investigation']], width=bar_width, edgecolor='grey', label='logCPM' if idx == 0 else "", hatch='/', alpha=0.7)
        plt.bar(r3[idx], row['counts_PSM_no_logcpm'], color=all_cancers_color_dict[row['cancer_investigation']], width=bar_width, edgecolor='grey', label='Raw' if idx == 0 else "", alpha=0.5)

    # Adding labels and legend
    plt.xlabel('Cancer Type', fontweight='bold')
    plt.xticks([r + bar_width for r in range(len(merged_counts))], merged_counts['cancer_investigation'].tolist(), rotation=90)
    plt.ylabel('Count of Significant p-values')
    plt.legend()
    plt.title(title)

    plt.show()


def paired_four_barplot_for_PSM_4analyses(results1, results2, results3, results4, all_cancers_color_dict, title):

    counts1 = results1[results1['p_value_fdr'] < 0.05].groupby('cancer_investigation')['fungus_id'].count().reset_index(name='counts1')
    counts2 = results2[results2['p_value_fdr'] < 0.05].groupby('cancer_investigation')['fungus_id'].count().reset_index(name='counts2')
    counts3 = results3[results3['p_value_fdr'] < 0.05].groupby('cancer_investigation')['fungus_id'].count().reset_index(name='counts3')
    counts4 = results4[results4['p_value_fdr'] < 0.05].groupby('cancer_investigation')['fungus_id'].count().reset_index(name='counts4')

    merged_counts = counts1.merge(counts2, on='cancer_investigation', how='outer').merge(counts3, on='cancer_investigation', how='outer').merge(counts4, on='cancer_investigation', how='outer').fillna(0)

    bar_width = 0.2
    r1 = np.arange(len(merged_counts))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]
    r4 = [x + bar_width for x in r3]

    plt.figure(figsize=(15, 6))

    for idx, row in merged_counts.iterrows():
        plt.bar(r1[idx], row['counts1'], color=all_cancers_color_dict[row['cancer_investigation']], width=bar_width, edgecolor='grey', label='Raw' if idx == 0 else "", alpha=1)
        plt.bar(r2[idx], row['counts2'], color=all_cancers_color_dict[row['cancer_investigation']], width=bar_width, edgecolor='grey', label='logCPM+Norm' if idx == 0 else "", hatch='//', alpha=0.8)
        plt.bar(r3[idx], row['counts3'], color=all_cancers_color_dict[row['cancer_investigation']], width=bar_width, edgecolor='grey', label='relBW' if idx == 0 else "", hatch='..', alpha=0.7)
        plt.bar(r4[idx], row['counts4'], color=all_cancers_color_dict[row['cancer_investigation']], width=bar_width, edgecolor='grey', label='mmuphinBW' if idx == 0 else "", hatch='\\', alpha=0.6)

    plt.xlabel('Cancer Type', fontweight='bold')
    plt.xticks([r + 1.5*bar_width for r in range(len(merged_counts))], merged_counts['cancer_investigation'].tolist(), rotation=90)
    plt.ylabel('Count of Significant p-values')
    plt.legend()
    plt.title(title)

    plt.tight_layout()
    plt.show()

def paired_four_barplot_for_PSM_4analyses_titles(results1,t1,results2,t2,results3,t3,results4,t4, all_cancers_color_dict, title):

    counts1 = results1[results1['p_value_fdr'] < 0.05].groupby('cancer_investigation')['fungus_id'].count().reset_index(name='counts1')
    counts2 = results2[results2['p_value_fdr'] < 0.05].groupby('cancer_investigation')['fungus_id'].count().reset_index(name='counts2')
    counts3 = results3[results3['p_value_fdr'] < 0.05].groupby('cancer_investigation')['fungus_id'].count().reset_index(name='counts3')
    counts4 = results4[results4['p_value_fdr'] < 0.05].groupby('cancer_investigation')['fungus_id'].count().reset_index(name='counts4')

    merged_counts = counts1.merge(counts2, on='cancer_investigation', how='outer').merge(counts3, on='cancer_investigation', how='outer').merge(counts4, on='cancer_investigation', how='outer').fillna(0)

    bar_width = 0.2
    r1 = np.arange(len(merged_counts))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]
    r4 = [x + bar_width for x in r3]

    plt.figure(figsize=(15, 6))

    for idx, row in merged_counts.iterrows():
        plt.bar(r1[idx], row['counts1'], color=all_cancers_color_dict[row['cancer_investigation']], width=bar_width, edgecolor='grey', label=t1 if idx == 0 else "", alpha=1)
        plt.bar(r2[idx], row['counts2'], color=all_cancers_color_dict[row['cancer_investigation']], width=bar_width, edgecolor='grey', label=t2 if idx == 0 else "", hatch='//', alpha=0.8)
        plt.bar(r3[idx], row['counts3'], color=all_cancers_color_dict[row['cancer_investigation']], width=bar_width, edgecolor='grey', label=t3 if idx == 0 else "", hatch='..', alpha=0.7)
        plt.bar(r4[idx], row['counts4'], color=all_cancers_color_dict[row['cancer_investigation']], width=bar_width, edgecolor='grey', label=t4 if idx == 0 else "", hatch='\\', alpha=0.6)

    plt.xlabel('Cancer Type', fontweight='bold')
    plt.xticks([r + 1.5*bar_width for r in range(len(merged_counts))], merged_counts['cancer_investigation'].tolist(), rotation=90)
    plt.ylabel('Count of Significant p-values')
    plt.legend()
    plt.title(title)

    plt.tight_layout()
    plt.show()
    
    
def apply_analysis_and_barplot(merged_df, columns_to_test, title = '', var1='BLACK OR AFRICAN AMERICAN', var2='WHITE'): #TODO - the proper importing here
    """
    columns_to_test are the fungal \ bacterial columns - usually I take the original counts_df before merge, columns 1:the end (sample name cols excluded) - previously list_of_bacteria
    assumes only race diagnosis for our matter, with same counfounders + need for dummies
    pass a pair of vars out of black, white, asian, default is black and white
    returns a result df with p value corrected fdr

    """
    merged_df_for_analysis = merged_df[merged_df['race'].isin([var1,var2])].copy()
    res = IPTW_weighting(merged_df_for_analysis, confounders_for_race, "race",need_dummies_race, list(poore_2020.columns[1:]), op1 = var1, op2 =var2,  weight = "IPTW")
    single_barplot_for_PSM(res,color_mapping, title )
    return res

def bar_plot_by_cancer(pt_merged_filtered,res_path):
    plt.rcParams.update({'font.size': 18})

    # Remove substring of 'TCGA-' from cancer types
    pt_merged_filtered.investigation = pt_merged_filtered.investigation.str.replace("TCGA-","")
    # Group data by 'investigation' and count
    cancer_counts = pt_merged_filtered.groupby('investigation').size().sort_values(ascending=False)

    # Plot
    ax = cancer_counts.plot(kind='bar', figsize=(15, 10), color='royalblue')

    for i, count in enumerate(cancer_counts):
        ax.text(i, count + 10, str(int(count)), ha='center', va='center')  # The "+10" is to position the label a little above the bar

    plt.xlabel('Cancer Type')
    plt.ylabel('Number of Samples')
    #plt.title('Number of Samples by Cancer Type')
    plt.tight_layout()
    plt.savefig(f'{res_path}/Figure_1B.png')
    plt.show()

def mean_total_count_per_sample(poore_2020_merged_sampesamples_sametaxa,poore_2022_wisoverlap_merged_sampesamples_sametaxa,salzberg_merged_samesamples_sametaxa,res_path):
  mean_data = pd.melt(pd.concat([
      poore_2020_merged_sampesamples_sametaxa[['investigation', 'mean_per_patient']].assign(dataset='P20'),
      poore_2022_wisoverlap_merged_sampesamples_sametaxa[['investigation', 'mean_per_patient']].assign(dataset='NH22'),
      salzberg_merged_samesamples_sametaxa[['investigation', 'mean_per_patient']].assign(dataset='G23')
  ]), id_vars=['investigation', 'dataset'], value_vars=['mean_per_patient'], value_name='value')

  mean_data.investigation = mean_data.investigation.str.replace("TCGA-","")

  plt.rcParams.update({'font.size': 18})

  plt.figure(figsize=(10, 8))
  sns.boxplot(x='investigation', y='value', hue='dataset', data=mean_data)
  #plt.title('Mean per Patient - Same samples and same taxa')
  plt.ylabel('Log(Mean Total Counts per Sample)')
  plt.xlabel('Cancer Type')
  plt.legend(title='Cohort', fontsize="10", loc ="upper right")
  legend = plt.legend(title='Cohort', fontsize="14", loc ="upper right", labelspacing=1.1)


   # Set the font size of the legend items
  for text in legend.get_texts():
        text.set_fontsize(18)  # You can set the font size to your desired value

  #plt.xticks(rotation=90)  # Rotate the x-axis labels
  plt.yscale('log')  # Optional: Apply logarithmic scale
  plt.savefig(f'{res_path}/Figure_2B.png')

  plt.show()

def total_count_per_cancer_type(poore_2020_merged_sampesamples_sametaxa,poore_2022_wisoverlap_merged_sampesamples_sametaxa,salzberg_merged_samesamples_sametaxa,res_path):
  total_data = pd.melt(pd.concat([
      poore_2020_merged_sampesamples_sametaxa[['investigation', 'total_counts']].assign(dataset='P20'),
      poore_2022_wisoverlap_merged_sampesamples_sametaxa[['investigation', 'total_counts']].assign(dataset='NH22'),
      salzberg_merged_samesamples_sametaxa[['investigation', 'total_counts']].assign(dataset='G23')
  ]), id_vars=['investigation', 'dataset'], value_vars=['total_counts'], value_name='value')
  total_data.investigation = total_data.investigation.str.replace("TCGA-","")

  plt.figure(figsize=(10, 6))
  sns.barplot(x='investigation', y='value', hue='dataset', data=total_data, ci = None)
  #plt.title('Total Counts by Investigation and Dataset - same samples and taxa')
  plt.ylabel('Log(Total Counts)')
  plt.xlabel('Cancer Type')
  plt.legend(title='Cohort', fontsize=16)
  plt.yscale('log')  # Optional: Apply logarithmic scale
  plt.savefig(f'{res_path}/Figure_2A.png')
  plt.show()

def scatter_totalcounts(resulting_df_allcancers, dataset_a,dataset_b):
    plt.figure(figsize=(8, 6))
    plt.rcParams.update({'font.size': 14})

    # Filter out zero values
    mask = (resulting_df_allcancers[dataset_a] != 0) & (resulting_df_allcancers[dataset_b] != 0)
    filtered_data = resulting_df_allcancers[mask]

    print(f'list of species that were removed because their total count was zero: {list(set(resulting_df_allcancers.index).difference(set(filtered_data.index)))}')

    x_data = filtered_data[dataset_a]
    y_data = filtered_data[dataset_b]

    #log_x_data = np.log10(x_data)
    #log_y_data = np.log10(y_data)
    #  regression
    #slope, intercept, r_value, p_value, std_err = linregress(log_x_data, log_y_data)

    #new regression
    slope, intercept, r_value, p_value, std_err = linregress(x_data, y_data)

    # Plotting
    plt.scatter(x_data, y_data,color='purple')
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel(dataset_b+" [log(Total Counts)]")
    plt.xlabel(dataset_a+" [log(Total Counts)]")
    #plt.title("Salzberg ~ Poore 2022 total count per Species")
    plt.xlim(min(pd.concat([x_data,y_data])),max(pd.concat([x_data,y_data])))
    plt.ylim(min(pd.concat([x_data,y_data])),max(pd.concat([x_data,y_data])))
    # square plot
    #plt.axis('square')
    #x_vals = np.array([x_data.min(), x_data.max()])
    #y_vals = np.exp(intercept + slope * np.log(x_vals))
    #plt.plot(x_vals, y_vals, '--', color='red')

    plt.annotate(f'\u03C1  = {round(r_value,2)} \nP-value = {"{:.2E}".format(p_value)}', xy=(0.1, 0.8), xycoords='axes fraction')
    plt.show()

def scatter_totalcounts_per_cancer(resulting_df_bycancer_flat, colors, dataset_a, dataset_b):
  plt.figure(figsize=(8, 6))
  cancer_types = resulting_df_bycancer_flat['investigation'].unique()

  for cancer_type in cancer_types:
      group = resulting_df_bycancer_flat[resulting_df_bycancer_flat['investigation'] == cancer_type]
      total_count_a = group[group.Dataset == dataset_a].reset_index().drop(['Dataset', 'investigation'], axis=1).transpose()[0]
      total_count_b = group[group.Dataset == dataset_b].reset_index().drop(['Dataset', 'investigation'], axis=1).transpose()[0]

      slope, intercept, r_value, p_value, std_err = linregress(total_count_a, total_count_b)

      #x_vals = np.linspace(group['Ratio_'+dataset_a].min(), group['Ratio_'+dataset_a].max(), 100)
      #y_vals = intercept + slope * x_vals
      plt.plot([0],[0],color=colors[cancer_type], label=f'{cancer_type.replace("TCGA-","")} (\u03C1={r_value:.3f}, P-value={"{:.2E}".format(p_value)}')
      #plt.annotate(f'\u03C1   = {round(r_value,2)} \nP-value<{p_value}', xy=(0.1, 0.8), xycoords='axes fraction')

      plt.scatter(total_count_a, total_count_b, color=colors[cancer_type], alpha=0.6)
      plt.yscale("log")
      plt.xscale("log")
  plt.xlabel(dataset_a+' [log(Total Counts)]')
  plt.ylabel(dataset_b+' [log(Total Counts)]')
  plt.legend(loc="upper left",fontsize="10")
  plt.show()


def plot_ratios(merged_ratios, dataset_a, dataset_b):
    epsilon = 0
    plt.rcParams.update({'font.size': 18})

    #merged_ratios['Log_'+dataset_a] = np.log10(merged_ratios['Ratio_'+dataset_a] + epsilon)
    #merged_ratios['Log_'+dataset_b] = np.log10(merged_ratios['Ratio_'+dataset_b] + epsilon)

    plt.figure(figsize=(12, 8))
    plt.scatter(merged_ratios['Ratio_'+dataset_a], merged_ratios['Ratio_'+dataset_b], alpha=0.6, color='purple')
    plt.yscale("log")
    plt.xscale("log")

    # Linear Regression
    slope, intercept, r_value,  p_value, std_err = linregress(merged_ratios['Ratio_'+dataset_a], merged_ratios['Ratio_'+dataset_b])
    print(r_value)
    if(p_value == 0):
      p_value = 7.14E-293
    #x_vals = np.linspace(merged_ratios['Ratio_'+dataset_a].min(), merged_ratios['Ratio_'+dataset_a].max(), 100)
    #y_vals = intercept + slope * x_vals
    #plt.plot(x_vals, y_vals, '--',color='red')
    if(p_value == 7.14E-293):
      plt.annotate(f'\u03C1   = {round(r_value,2)} \nP-value<{p_value}', xy=(0.1, 0.8), xycoords='axes fraction')
    else:
      plt.annotate(f'\u03C1  = {round(r_value,2)} \nP-value = {"{:.2E}".format(p_value)}', xy=(0.1, 0.8), xycoords='axes fraction')

    #plt.title('Log Scale Scatter Plot of Pairwise Ratios between Species Counts')
    plt.xlabel(dataset_a+ ' log(Pairwise ratio of species)')
    plt.ylabel(dataset_b+ ' log(Pairwise ratio of species)')
    #plt.legend()
    plt.show()

def pearson_ind_corellations(list_corr_P20_G23, list_corr_NH22_G23):
  U1, p_value = mannwhitneyu(list_corr_P20_G23, list_corr_NH22_G23)
  sns.distplot(list_corr_P20_G23, hist=False, label = "P20 vs. G23", color='red')
  sns.distplot(list_corr_NH22_G23, hist=False, label="NH22 vs. G23", color='purple')
  plt.ylabel('Density')
  plt.xlabel('Pearson Correlation')
  plt.annotate(f'P-value = {"{:.2E}".format(p_value)}', xy=(0.25, 0.02), xycoords='axes fraction')

  plt.legend(loc='upper right')
  plt.show()

def plot_ratios_bycancer_v2(merged_ratios, colors, dataset_a, dataset_b):

    plt.figure(figsize=(12, 8))
    cancer_types = merged_ratios['Cancer'].unique()

    for cancer_type in cancer_types:
        group = merged_ratios[merged_ratios['Cancer'] == cancer_type]

        slope, intercept, r_value, p_value, std_err = linregress(group['Ratio_'+dataset_a], group['Ratio_'+dataset_b])
        #x_vals = np.linspace(group['Ratio_'+dataset_a].min(), group['Ratio_'+dataset_a].max(), 100)
        #y_vals = intercept + slope * x_vals
        plt.plot([0],[0],color=colors[cancer_type], label=f'{cancer_type.replace("TCGA-","")} (\u03C1={r_value:.3f}, P-value={"{:.2E}".format(p_value)}')
        #plt.annotate(f'\u03C1   = {round(r_value,2)} \nP-value<{p_value}', xy=(0.1, 0.8), xycoords='axes fraction')

        plt.scatter(group['Ratio_'+dataset_a], group['Ratio_'+dataset_b], color=colors[cancer_type], alpha=0.6)
        plt.yscale("log")
        plt.xscale("log")

    #plt.title('Log Scale Scatter Plot of Pairwise Ratios between Species Counts ')
    plt.xlabel(dataset_a+' log(Pairwise ratio of species)')
    plt.ylabel(dataset_b+' log(Pairwise ratio of species)')
    plt.legend(loc="lower right",fontsize="12")
    plt.show()

def scatter_totalcounts_per_sample_per_cancer(total_count_data, dataset_a,dataset_b,colors):
  plt.figure(figsize=(8, 6))
  plt.rcParams.update({'font.size': 14})

  #x_data = total_count_data[total_count_data.dataset ==dataset_a]['value']
  #y_data = total_count_data[total_count_data.dataset ==dataset_b]['value']

  cancer_types = total_count_data['investigation'].unique()

  for cancer_type in cancer_types:
      group = total_count_data[total_count_data['investigation'] == cancer_type]
      total_count_a = group[group.dataset == dataset_a]['value']
      total_count_b = group[group.dataset == dataset_b]['value']

      slope, intercept, r_value, p_value, std_err = linregress(total_count_a, total_count_b)

      #x_vals = np.linspace(group['Ratio_'+dataset_a].min(), group['Ratio_'+dataset_a].max(), 100)
      #y_vals = intercept + slope * x_vals
      plt.plot([0],[0],color=colors[cancer_type], label=f'{cancer_type.replace("TCGA-","")} (\u03C1={r_value:.3f}, P-value={"{:.2E}".format(p_value)}')
      #plt.annotate(f'\u03C1   = {round(r_value,2)} \nP-value<{p_value}', xy=(0.1, 0.8), xycoords='axes fraction')

      plt.scatter(total_count_a, total_count_b, color=colors[cancer_type], alpha=0.6)
      plt.yscale("log")
      plt.xscale("log")
  plt.xlabel(dataset_a+' [log(Total Counts)]')
  plt.ylabel(dataset_b+' [log(Total Counts)]')
  plt.legend(loc="upper left",fontsize="10")
  plt.show()

def scatter_totalcounts_per_sample(total_count_data, dataset_a,dataset_b):
    plt.figure(figsize=(8, 6))
    plt.rcParams.update({'font.size': 14})

    x_data = total_count_data[total_count_data.dataset ==dataset_a]['value']
    y_data = total_count_data[total_count_data.dataset ==dataset_b]['value']

    #new regression
    slope, intercept, r_value, p_value, std_err = linregress(x_data, y_data)

    # Plotting
    plt.scatter(x_data, y_data,color='purple')
    plt.yscale("log")
    plt.xscale("log")
    plt.ylabel(dataset_b+" [log(Total Counts)]")
    plt.xlabel(dataset_a+" [log(Total Counts)]")
    #plt.title("Salzberg ~ Poore 2022 total count per Species")
    plt.xlim(min(pd.concat([x_data,y_data])),max(pd.concat([x_data,y_data])))
    plt.ylim(min(pd.concat([x_data,y_data])),max(pd.concat([x_data,y_data])))
    # square plot
    #plt.axis('square')
    #x_vals = np.array([x_data.min(), x_data.max()])
    #y_vals = np.exp(intercept + slope * np.log(x_vals))
    #plt.plot(x_vals, y_vals, '--', color='red')

    plt.annotate(f'\u03C1  = {round(r_value,2)} \nP-value = {"{:.2E}".format(p_value)}', xy=(0.1, 0.8), xycoords='axes fraction')
    plt.show()

def scatter_per_species_per_cancer(df1, df2, species,colors, cancer_type,ind, cancer_col='investigation', title = 'Counts',ds1_name = "NH22", ds2_name = "G23"):
    """
    Scatter plot the relative abundance of a given species in two datasets and
    draw regression lines for each cancer type, using a log scale.

    Parameters:
    - df1: DataFrame containing data from the first dataset.
    - df2: DataFrame containing data from the second dataset.
    - species: The species column name for which the scatter plot is to be generated.
    - cancer_col: Column name that indicates cancer type.
    """

    if (cancer_type != ''):
      df1 = df1[df1[cancer_col] == cancer_type]
      df2 = df2[df2[cancer_col] == cancer_type]

    plt.figure(figsize=(6, 4))

    # Get unique cancer types
    #unique_cancers = df1[cancer_col].unique()

    # Define colors for each cancer type
    #for cancer in unique_cancers:

    x_data = df1[species]
    y_data = df2[species]

    # Filter non-zero values for log scale
    mask = (x_data > 0) & (y_data > 0)
    #x_data_filtered = np.log10(x_data[mask])
    #y_data_filtered = np.log10(y_data[mask])

    if (cancer_type != ''):
      plt.scatter(x_data[mask], y_data[mask], color=colors[cancer_type], alpha=0.6)
    else:
      plt.scatter(x_data[mask], y_data[mask], color='purple', alpha=0.6)

    plt.yscale("log")
    plt.xscale("log")

    # Regression
    if len(x_data[mask]) > 0:  # Check if there's any data after filtering
        slope, intercept, r_value, p_value, std_err = linregress(x_data[mask], y_data[mask])
        x_vals = np.linspace(x_data[mask].min(), x_data[mask].max(), 100)
        y_vals = intercept + slope * x_vals
        #plt.plot([0], [0], color=colors[cancer_type], alpha=0.6, label=f"{cancer_type} (R^2 = {r_value:.3f})")
        if (cancer_type != ''):
          plt.plot([0],[0],color=colors[cancer_type], label=f'{cancer_type.replace("TCGA-","")} (\u03C1={r_value:.3f}, P-value={"{:.2E}".format(p_value)})')
        else:
          plt.plot([0],[0],color='purple', label=f'(\u03C1={r_value:.3f}, P-value={"{:.2E}".format(p_value)})')

    plt.xlabel(f'{ds1_name} Log(total count per sample)',fontsize=16)
    plt.ylabel(f'{ds2_name} Log(total count per sample)',fontsize=16)
    plt.title(f'{species}(count={len(y_data[mask])},rank={ind+1})',fontsize=18)
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)
    plt.legend(fontsize=14,loc="upper left")
    plt.show()

def top_5_fungi_per_cancer(poore_2022_wisoverlap_merged_sampesamples_sametaxa,salzberg_merged_samesamples_sametaxa,cancer_type, all_together_same_taxa,colors,df1_name, df2_name):

  cancer_types = salzberg_merged_samesamples_sametaxa['investigation'].unique()

  if (cancer_type != ''):
    total_counts_salzberg = salzberg_merged_samesamples_sametaxa[salzberg_merged_samesamples_sametaxa.investigation == cancer_type][all_together_same_taxa].fillna(0).astype(bool).sum(axis=0)
  else:
    total_counts_salzberg = salzberg_merged_samesamples_sametaxa[all_together_same_taxa].fillna(0).astype(bool).sum(axis=0)

  # the top 5 fungus types with highest total count
  top_5_fungi = total_counts_salzberg.nlargest(5).index.tolist()

  for ind, fungi in enumerate(top_5_fungi):
    scatter_per_species_per_cancer(poore_2022_wisoverlap_merged_sampesamples_sametaxa, salzberg_merged_samesamples_sametaxa, fungi, colors,cancer_type,ind,ds1_name=df1_name, ds2_name=df2_name)


#Function for voom-snm
def bar_plot_per_cancer_type(df_nromalized_with_cancer,cancer_name,genome_id,only_rna_seq,experimental_strategy_name, font_size=None):
  if (only_rna_seq):
    cancer_counts = round(df_nromalized_with_cancer[(df_nromalized_with_cancer.experimental_strategy.isin(experimental_strategy_name)) &
                                                (df_nromalized_with_cancer.investigation == cancer_name) &
                                                (df_nromalized_with_cancer.sample_type =='Solid Tissue Normal')][genome_id],4).value_counts().sort_index()

    other_cancers_counts = round(df_nromalized_with_cancer[(df_nromalized_with_cancer.experimental_strategy.isin(experimental_strategy_name)) &
                                                          (df_nromalized_with_cancer.investigation == cancer_name) &
                                                          (df_nromalized_with_cancer.sample_type=='Primary Tumor')][genome_id],4).value_counts().sort_index()
  else:
    cancer_counts = round(df_nromalized_with_cancer[(df_nromalized_with_cancer.investigation == cancer_name) &
                                                    (df_nromalized_with_cancer.sample_type =='Solid Tissue Normal')][genome_id],4).value_counts().sort_index()

    other_cancers_counts = round(df_nromalized_with_cancer[(df_nromalized_with_cancer.investigation == cancer_name) &
                                                          (df_nromalized_with_cancer.sample_type=='Primary Tumor')][genome_id],4).value_counts().sort_index()

  df_val_counts_cancer = pd.DataFrame(cancer_counts)
  df_value_counts_reset_cancer = df_val_counts_cancer.reset_index()
  df_value_counts_reset_cancer.columns = ['unique_values', 'counts']

  df_val_counts = pd.DataFrame(other_cancers_counts)
  df_value_counts_reset = df_val_counts.reset_index()
  df_value_counts_reset.columns = ['unique_values', 'counts']

  print("No. Solid Tissue Normal Samples: ", df_value_counts_reset_cancer.sum())
  print("No. Primary Tumor Samples: ", df_value_counts_reset.sum())

  compare = df_value_counts_reset.merge(df_value_counts_reset_cancer, on='unique_values', how='outer').sort_values(by='unique_values')

  #df_nromalized_with_cancer[genome_id].value_counts()
  #df_nromalized_with_cancer[(df_nromalized_with_cancer[genome_id] ==0)][genome_id].value_counts()
  x_axis = np.arange(len(compare.unique_values))
  plt.figure(figsize=(24,12), dpi = 1000)

  plt.xlabel('Voom-SNM normalized values')
  plt.ylabel('Number of Samples')
  plt.bar(x_axis -0.2, compare.counts_y, width=0.4, label = cancer_name+' '+'Solid Tissue Normal')
  plt.bar(x_axis +0.2, compare.counts_x, width=0.4, label = cancer_name+' '+'Primary Tumor')

  plt.xticks(x_axis, compare.unique_values, rotation='vertical', fontsize=font_size)
  plt.yticks(fontsize=font_size)

  plt.legend()

  # Display

  plt.show()
