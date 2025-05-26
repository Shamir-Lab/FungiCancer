import pandas as pd
import numpy as np

import plotting
STAGES_DICT ={"Stage I" :1, "Stage IIA" :2, "Stage IIB" :2, "Stage III" : 3, "Stage II":2, "Stage IV":4, "Stage IIIA":3, "Stage IB":1, "Stage IA": 1, "Stage IIIB":3, "Stage IIIC":3, "Stage IVA":4, "Stage IIC":2, "Stage IS":1, "Stage IVC":4, "!! or II NOS": np.nan, "Stage X":np.nan, "Stage IVB":4, "Stage 0":np.nan, "Stage Tis":np.nan}

def extract_kingdom_genus(col_name):
    kingdom = next((item.split('__')[1] for item in col_name.split('.') if item.startswith('k__')), None)
    genus = next((item.split('__')[1] for item in col_name.split('.') if item.startswith('g__')), None)
    return kingdom, genus

def extract_taxonomy_info(df):
    """Extract Kingdom and Genus information from poore_2020 DataFrame."""
    columns = list(df.columns)[1:]
    kingdoms, genera = zip(*[extract_kingdom_genus(col) for col in columns])

    return pd.DataFrame({
        'Kingdom': kingdoms,
        'Genus': genera
    })

def read_dataframes(root_path):
    """Load all dataframes and return them."""
    files = {
        'metadata_poore_2022': 'files_for_fungi_analysis/metadata_genus_WIS_overlapping_fungi_bacteria_14494samples.tsv',
        'bacteria_metadata': 'files_for_fungi_analysis/Metadata-TCGA-All-18116-Samples.csv',
        'salz_hnsc': 'files_for_fungi_analysis/TableS9_HNSC_all.xlsx',
        'salz_brca': 'files_for_fungi_analysis/TableS10_BRCA_WGS.xlsx',
        'salz_blca': 'files_for_fungi_analysis/TableS8_BLCA.all.xlsx',
        'poore_2022_wisoverlap': 'files_for_fungi_analysis/count_data_genus_raw_WIS_overlapping_fungi_bacteria_14494samples.tsv',
        'genus_poore_2022': 'files_for_fungi_analysis/taxonomy_table_WIS_overlapping_fungi_bacteria.tsv',
        'poore_2020': 'files_for_fungi_analysis/Kraken-TCGA-Raw-Data-17625-Samples.csv',
        'poore_2022_raw_fungi_counts': 'files_for_fungi_analysis/count_data_fungi_decontaminated_raw.tsv',
        'poore_2022_fungi_species': 'files_for_fungi_analysis/taxonomy_table_rep200.tsv',
        'poore_fungi_voom_snm':'files_for_fungi_analysis/count_data_fungi_decontaminated_voom_snm_corrected.tsv',
        'for_bmi':'files_for_fungi_analysis/for_bmi.xlsx',
        'poore_fungi_voom_snm':'files_for_fungi_analysis/count_data_fungi_decontaminated_voom_snm_corrected.tsv'
    }
        
    dataframes = {}
    for name, file in files.items():
        if file.endswith('.tsv'):
            dataframes[name] = pd.read_csv(f'{root_path}/{file}', sep='\t')
        elif file.endswith('.csv'):
            dataframes[name] = pd.read_csv(f'{root_path}/{file}')
        elif file.endswith('.xlsx'):
            dataframes[name] = pd.read_excel(f'{root_path}/{file}')
                # Transformations on poore_2020 after reading the dataframe
    poore_2020 = dataframes['poore_2020']

    # Extract kingdom and genus
    poore_2020_cols_old_before_transform = list(poore_2020.columns)[1:]
    kingdoms, genera = zip(*[extract_kingdom_genus(col) for col in poore_2020_cols_old_before_transform])

    poore_2020_taxonomy_df = pd.DataFrame({
        'Kingdom': kingdoms,
        'Genus': genera
    })

    # Rename columns to simplify taxa names
    samples_col_poore_2020 = poore_2020.columns[0]
    genus_to_column_mapping = {taxon.split('.')[-1].split('__')[-1]: taxon for taxon in poore_2020_cols_old_before_transform}
    column_to_genus_mapping = {v: k for k, v in genus_to_column_mapping.items()}
    new_columns = {col: column_to_genus_mapping.get(col, col) for col in poore_2020.columns if col != samples_col_poore_2020}
    poore_2020 = poore_2020.rename(columns=new_columns)
    poore_2020 = poore_2020.rename(columns={"Unnamed: 0": "Sample", "3alikevirus": "virus3alike"})

    # Store back the transformed dataframe
    dataframes['poore_2020'] = poore_2020
    dataframes['poore_2020_taxonomy'] = poore_2020_taxonomy_df

    # Reading and transforming for_bmi
    for_bmi = dataframes['for_bmi']
    for_bmi["height"] = pd.to_numeric(for_bmi["height"], errors="coerce")
    for_bmi["weight"] = pd.to_numeric(for_bmi["weight"], errors="coerce")
    for_bmi["number_pack_years_smoked"] = pd.to_numeric(for_bmi["number_pack_years_smoked"], errors="coerce")
    dataframes['for_bmi'] = for_bmi



    return dataframes


def replace_not_available_with_nan(dataframes):
    """Replaces 'Not available' with np.nan in all provided dataframes."""
    for key in dataframes:
        dataframes[key] = dataframes[key].replace('Not available', np.nan, inplace=True)
    return dataframes

def simplify_salz_df(df):
    """Simplify column names for Salzberg dataframes."""
    df.rename(columns=lambda x: x.split('_')[-1] if x.startswith('g_') else x, inplace=True)
    df.columns = [col.replace(' ', '_') for col in df.columns]
    return df

def compute_over_65(df):
    """Compute the 'over_65' column based on 'age_at_diagnosis'."""
    over_65 = (df['age_at_diagnosis'] > 70).astype(int)
    over_65[df['age_at_diagnosis'].isna()] = np.nan
    df['over_65'] = over_65
    return df

def assign_stage_numbers(df):
    """Assign numbers to the stages."""
    df['stage_numbered'] = df['pathologic_stage_label'].map(STAGES_DICT)
    return df

def add_bmi_and_obesity_cols(df, dataframes):
    """
    Adds BMI and obesity columns to the dataframe based on height and weight.

    Parameters:
    - df (DataFrame): The input dataframe which should have columns 'height' and 'weight'.

    Returns:
    - DataFrame: A modified dataframe with added 'BMI' and 'obese' columns.
    """
    
    df = pd.merge(df, dataframes["for_bmi"][['bcr_patient_uuid', 'height', 'weight']], how='left', left_on='cgc_case_uuid', right_on='bcr_patient_uuid')
    
    # Compute the BMI
    df.loc[(~df['height'].isna()) & (~df['weight'].isna()), 'BMI'] = df['weight'] / (df['height']/100)**2

    df["height"] = pd.to_numeric(df["height"], errors="coerce")
    df["weight"] = pd.to_numeric(df["weight"], errors="coerce")

    # Compute the BMI
    df.loc[(~df['height'].isna()) & (~df['weight'].isna()), 'BMI'] = df['weight'] / (df['height']/100)**2
    
    # Compute the obesity column
    df.loc[(~df['BMI'].isna()), 'obese'] = (df["BMI"] > 30).astype(int)


    return df

def replace_not_available_with_nan(dataframes):
    """Replaces 'Not available' with np.nan in all provided dataframes."""
    for key in dataframes:
        dataframes[key].replace("Not available", np.nan,inplace=True)
    return dataframes

def get_overlapping_taxa(dataframes):
    """Return overlapping taxa between two dataframes."""
    overlap_taxa_salz_poore2020 = list(set(dataframes['bacteria_counts_salzberg'].columns[1:]).intersection(dataframes['poore_2020'].columns[1:]))
    overlap_taxa_poore2020_poore2022 = list(set(dataframes['poore_2020'].columns[1:]).intersection(dataframes['poore_2022_wisoverlap'].columns[1:]))
    all_together_same_taxa = list(set(overlap_taxa_salz_poore2020).intersection(overlap_taxa_poore2020_poore2022))
    return all_together_same_taxa

def standardize_sample_column_name(dfs, col_name='Sample', alternative_col_name='sampleid'):
    """Standardize the sample column name across all dataframes."""
    for key, df in dfs.items():
        if alternative_col_name in df.columns:
            dfs[key] = df.rename(columns={alternative_col_name: col_name})
    return dfs

def filter_dataframes_on_overlap_using_filename(dfs, filenames_list, taxa_list, meta_cols_dict, filename_col='filename'):
    """Filter dataframes based on overlapping filenames and taxa."""
    filtered_dfs = {}
    for key, df in dfs.items():
        meta_cols = meta_cols_dict.get(key, [])  # Get the metadata columns for the specific dataframe
        columns = [filename_col] + taxa_list + meta_cols
        filtered_dfs[key] = df[df[filename_col].isin(filenames_list)][columns]
    return filtered_dfs

def get_overlapping_samples(dataframes,sample_col = 'filename' ):
    salz_and_poore_2020_same_filenames = set(dataframes['poore_2020_merged'][sample_col]).intersection(dataframes['salzberg_merged'][sample_col])

    final_overlapping_filenames = salz_and_poore_2020_same_filenames.intersection(dataframes['poore_2022_wisoverlap_merged'][sample_col])
    return final_overlapping_filenames

def raw_fungi_for_analysis(dataframes,specimen_type):
    """
    returns the for analysis and for BC files
    """
    raw_fungi_table_merged = dataframes['poore_2022_raw_fungi_counts_merged']
    pt_merged = raw_fungi_table_merged[raw_fungi_table_merged["sample_type"] =="Primary Tumor"]
    pt_number = len(pt_merged)
    print(f"Only Primary Tumor Samples: {pt_number}")
    pt_merged_rnaseq = pt_merged[pt_merged["experimental_strategy"] == specimen_type]
    pt_merged_rnaseq = pt_merged_rnaseq[pt_merged_rnaseq['data_submitting_center_label'].notna()]
    pt_rnaseq_number = len(pt_merged_rnaseq)
    print(f"Only RNA-Seq Samples: {pt_rnaseq_number}")
    columns_to_analyze = ["gender", "age_at_diagnosis", "race", "histological_diagnosis_label", "stage_numbered"]
    pt_merged_rnaseq_nonanconfs = pt_merged_rnaseq.dropna(subset=columns_to_analyze)
    pt_rnaseq_nonanconfs_number = len(pt_merged_rnaseq_nonanconfs)
    print(f"Samples with all clinical variables: {pt_rnaseq_nonanconfs_number}")
    return pt_merged_rnaseq_nonanconfs, pt_merged_rnaseq

def read_flow(root_path):
    dataframes = read_dataframes(root_path)
    dataframes = replace_not_available_with_nan(dataframes)
    dataframes = standardize_sample_column_name(dataframes)
    # Renaming for better readability
    dataframes['salz_blca'].rename(columns={"Unnamed: 0": "Sample"}, inplace=True)
    dataframes['bacteria_metadata'].rename(columns={"Unnamed: 0": "Sample"}, inplace=True)
    # Simplify Salzberg dataframes
    salz_dfs = ['salz_hnsc', 'salz_brca', 'salz_blca']
    for df_name in salz_dfs:
        dataframes[df_name] = simplify_salz_df(dataframes[df_name])

    # Merge Salzberg dataframes
    dataframes['bacteria_counts_salzberg'] = pd.concat([dataframes[name] for name in salz_dfs])

    # Merge the big tables
    merged_dfs = {
        'poore_2020_merged': ('poore_2020', 'bacteria_metadata'),
        'poore_2022_wisoverlap_merged': ('poore_2022_wisoverlap', 'metadata_poore_2022'),
        'salzberg_merged': ('bacteria_counts_salzberg', 'bacteria_metadata'),
        'poore_2022_raw_fungi_counts_merged': ('poore_2022_raw_fungi_counts', 'metadata_poore_2022'),
        'poore_fungi_voom_snm_merged':('poore_fungi_voom_snm','metadata_poore_2022')

    }

    for merged_name, (df1_name, df2_name) in merged_dfs.items(): # Removed the on_col here
        dataframes[merged_name] = dataframes[df1_name].merge(dataframes[df2_name], on='Sample', how='left')


    # Apply the over_65 and stages assignment
    for merged_df in merged_dfs.keys():
        dataframes[merged_df] = compute_over_65(dataframes[merged_df])
        dataframes[merged_df] = assign_stage_numbers(dataframes[merged_df])

    dataframes['poore_2020_taxonomy'] = extract_taxonomy_info(dataframes['poore_2020'])

    #Same sampe + taxa:
    samples_overlap = get_overlapping_samples(dataframes,sample_col = 'filename' )
    all_together_same_taxa = get_overlapping_taxa(dataframes)

    # Metadata columns
    metadata_cols_old = dataframes['salzberg_merged'].columns[-43:].tolist()
    metadata_cols_new = dataframes["poore_2022_wisoverlap_merged"].columns[-43:].to_list() #metdata column names in poore 2022
    # Filter dataframes based on overlap
    dfs_to_filter = {
        'salzberg_merged': dataframes['salzberg_merged'],
        'poore_2020_merged': dataframes['poore_2020_merged'],
        'poore_2022_wisoverlap_merged': dataframes['poore_2022_wisoverlap_merged']
    }

    metadata_cols_dict = {
        'salzberg_merged': metadata_cols_old,
        'poore_2020_merged': metadata_cols_old,
        'poore_2022_wisoverlap_merged': metadata_cols_new
    }

    sampesamples_sametaxa_dfs = filter_dataframes_on_overlap_using_filename(dfs_to_filter, samples_overlap, all_together_same_taxa, metadata_cols_dict)
    return dataframes, sampesamples_sametaxa_dfs


def load_batch_corrected_files(dataframes, root_path,specimen_type, cohort_definition):
    dfs = {}
    for method_transform in ['relative_abundance', 'clr_c', 'clr_offset', 'log_cpm_quantile']:
        #for method_bc in ['bmc', 'combat', 'mmuphin']:
        for method_bc in ['bmc', 'combat', 'mmuphin','plsda']:
            key_name = f'{method_transform}_{method_bc}'
            print(f'generate batch corrected and normalized file using {key_name}')

            if key_name not in ['relative_abundance_plsda', 'log_cpm_quantile_plsda']:
                dfs[key_name] = pd.read_csv(f'{root_path}/files_for_fungi_analysis/upated_batch_corrected_files/adjusted_{specimen_type}_{cohort_definition}_{method_transform}_{method_bc}.csv')
                #dfs[key_name] = pd.read_csv(f'{root_path}/files_for_fungi_analysis/new_batch_corrected_files/adjuste_fungi_onlypt_{method_transform}_{method_bc}.csv')
                dfs[key_name].rename(columns={"Unnamed: 0": "Sample"}, inplace=True)
                dfs[key_name] = dfs[key_name].merge(dataframes["metadata_poore_2022"], on="Sample", how="left")

                dfs[key_name] = compute_over_65(dfs[key_name])
                dfs[key_name] = assign_stage_numbers(dfs[key_name])
    return dfs


