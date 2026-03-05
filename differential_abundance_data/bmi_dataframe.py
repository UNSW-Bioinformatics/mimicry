from pathlib import Path
from scipy import stats

import pandas as pd
import numpy as np
import re

def get_threshold(pvalues):
    return stats.false_discovery_control(pvalues)

def remove_bracket_details(organism):
    return re.sub("[\(\[].*?[\)\]]", "", organism)

def rename_header(df_to_rename):
    df_to_rename.rename(
        columns={
            "mean OCO": "mean healthy weight",
            "mean OLX": "mean obesity",
            "pvalue OLX to control": "pvalue Obesity to healthy weight",
            "ratio OLX to control": "ratio of Obesity to healthy weight",
            "negative_log_pval_olx": "negative_log_pval_obesity",
            "log2FC_OLX": "log2FC_Obesity",
            "mean OLN": "mean overweight",
            "pvalue OLN to control": "pvalue Overweight to healthy weight",
            "ratio OLN to control": "ratio of Overweight to healthy weight",
            "negative_log_pval_oln": "negative_log_pval_overweight",
            "log2FC_OLN": "log2FC_Overweight"
        },
        inplace=True
    )
    return df_to_rename

app_dir = Path(__file__).parent

### Liver ###
df_bmi_liver = pd.read_csv(app_dir / "bmi_liver_match.csv", index_col=0)
df_bmi_liver.head()

df_bmi_liver['log2FC_OCI'] = np.log2(df_bmi_liver['ratio OCI to control'])
df_bmi_liver['log2FC_OLN'] = np.log2(df_bmi_liver['ratio OLN to control'])
df_bmi_liver['log2FC_OLX'] = np.log2(df_bmi_liver['ratio OLX to control'])

df_bmi_liver['negative_log_pval_oci'] = np.log10(get_threshold(df_bmi_liver['pvalue OCI to control'])) * (-1)
df_bmi_liver['negative_log_pval_oln'] = np.log10(get_threshold(df_bmi_liver['pvalue OLN to control'])) * (-1)
df_bmi_liver['negative_log_pval_olx'] = np.log10(get_threshold(df_bmi_liver['pvalue OLX to control'])) * (-1)

df_bmi_liver['organism'] = df_bmi_liver['species'].apply(remove_bracket_details)

df_bmi_liver = rename_header(df_bmi_liver)

### Breast ###
df_bmi_breast = pd.read_csv(app_dir / "bmi_breast_match.csv", index_col=0)
df_bmi_breast.head()

df_bmi_breast['log2FC_OCI'] = np.log2(df_bmi_breast['ratio OCI to control'])
df_bmi_breast['log2FC_OLN'] = np.log2(df_bmi_breast['ratio OLN to control'])
df_bmi_breast['log2FC_OLX'] = np.log2(df_bmi_breast['ratio OLX to control'])

df_bmi_breast['negative_log_pval_oci'] = np.log10(get_threshold(df_bmi_breast['pvalue OCI to control'])) * (-1)
df_bmi_breast['negative_log_pval_oln'] = np.log10(get_threshold(df_bmi_breast['pvalue OLN to control'])) * (-1)
df_bmi_breast['negative_log_pval_olx'] = np.log10(get_threshold(df_bmi_breast['pvalue OLX to control'])) * (-1)

df_bmi_breast['organism'] = df_bmi_breast['species'].apply(remove_bracket_details)

df_bmi_breast = rename_header(df_bmi_breast)

### Kidney ###
df_bmi_kidney = pd.read_csv(app_dir / "bmi_kidney_match.csv", index_col=0)
df_bmi_kidney.head()

df_bmi_kidney['log2FC_OCI'] = np.log2(df_bmi_kidney['ratio OCI to control'])
df_bmi_kidney['log2FC_OLN'] = np.log2(df_bmi_kidney['ratio OLN to control'])
df_bmi_kidney['log2FC_OLX'] = np.log2(df_bmi_kidney['ratio OLX to control'])

df_bmi_kidney['negative_log_pval_oci'] = np.log10(get_threshold(df_bmi_kidney['pvalue OCI to control'])) * (-1)
df_bmi_kidney['negative_log_pval_oln'] = np.log10(get_threshold(df_bmi_kidney['pvalue OLN to control'])) * (-1)
df_bmi_kidney['negative_log_pval_olx'] = np.log10(get_threshold(df_bmi_kidney['pvalue OLX to control'])) * (-1)

df_bmi_kidney['organism'] = df_bmi_kidney['species'].apply(remove_bracket_details)

df_bmi_kidney = rename_header(df_bmi_kidney)
