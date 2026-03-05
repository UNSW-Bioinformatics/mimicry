from pathlib import Path
from scipy import stats

import pandas as pd
import numpy as np

def get_threshold(pvalues):
    return stats.false_discovery_control(pvalues)

def replace_semicolon(organism):
    return organism.replace(";", "\n")

app_dir = Path(__file__).parent

### Oral ###
df_hcc_oral = pd.read_csv(app_dir / "patient_match_all_oral.csv", index_col=0)
df_hcc_oral.head()

df_hcc_oral['log2FC_OCI'] = np.log2(df_hcc_oral['ratio OCI to control'])
df_hcc_oral['log2FC_OLN'] = np.log2(df_hcc_oral['ratio OLN to control'])
df_hcc_oral['log2FC_OLX'] = np.log2(df_hcc_oral['ratio OLX to control'])

df_hcc_oral['negative_log_pval_oci'] = np.log10(get_threshold(df_hcc_oral['pvalue OCI to control'])) * (-1)
df_hcc_oral['negative_log_pval_oln'] = np.log10(get_threshold(df_hcc_oral['pvalue OLN to control'])) * (-1)
df_hcc_oral['negative_log_pval_olx'] = np.log10(get_threshold(df_hcc_oral['pvalue OLX to control'])) * (-1)

df_hcc_oral['organism'] = df_hcc_oral['species'].apply(replace_semicolon)

### Stool ###
df_hcc_stool = pd.read_csv(app_dir / "patient_match_all_stool.csv", index_col=0)
df_hcc_stool.head()

df_hcc_stool['log2FC_OCI'] = np.log2(df_hcc_stool['ratio OCI to control'])
df_hcc_stool['log2FC_OLN'] = np.log2(df_hcc_stool['ratio OLN to control'])
df_hcc_stool['log2FC_OLX'] = np.log2(df_hcc_stool['ratio OLX to control'])

df_hcc_stool['negative_log_pval_oci'] = np.log10(get_threshold(df_hcc_stool['pvalue OCI to control'])) * (-1)
df_hcc_stool['negative_log_pval_oln'] = np.log10(get_threshold(df_hcc_stool['pvalue OLN to control'])) * (-1)
df_hcc_stool['negative_log_pval_olx'] = np.log10(get_threshold(df_hcc_stool['pvalue OLX to control'])) * (-1)

df_hcc_stool['organism'] = df_hcc_stool['species'].apply(replace_semicolon)