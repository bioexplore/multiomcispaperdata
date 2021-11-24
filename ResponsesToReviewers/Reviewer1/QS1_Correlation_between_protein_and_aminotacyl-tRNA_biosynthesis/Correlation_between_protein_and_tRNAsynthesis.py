# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 14:06:53 2021

@author: bioex
"""

import pandas as pd
import numpy as np
from datetime import date
from scipy import stats # For doing the spearman correlation analysis
from statsmodels.stats import multitest # For doing the FDR correction of the p-value

# Read in dataset for protein and tRNA synthesis reaction fluxes seperately
protein_data = pd.read_excel("pvsm_new.xlsx",sheet_name="Protome",index_col=0)
flux_data_bk = pd.read_excel("tRNA_Fluxes_proteome.xlsx",sheet_name="Flux",index_col=3)

#protein_data = protein_data.iloc[:,1:10]
flux_data = flux_data_bk.iloc[:,3:12]
tmp_colnames = protein_data.columns
protein_data.columns = ['D'+name[2:] if 'C_' in name else name for name in tmp_colnames]
flux_data_T = flux_data.T
protein_data_T = protein_data.T

#Form the spearman correlation dataframe by using scipy.stats.spearman function
flux_names=flux_data_T.columns
p_names=protein_data_T.columns
names=[]
corrcoefs=[]
pvalues=[]

for name in flux_names:
    if name in p_names:
        flux = flux_data_T[name]
        protein = protein_data_T[name]        
        available_protein_num = sum(~protein.isna())
        print(name,available_protein_num)
        if available_protein_num<3 or flux[0]==0:
            continue
        #print('name:',name)
        corr,p=stats.spearmanr(flux,protein,nan_policy='omit')
        names.append(name)
        corrcoefs.append(corr)
        pvalues.append(p)
pvsflux_correlation_df = pd.DataFrame({'Name':names,'coef':corrcoefs,'pvalue':pvalues})

# Output the correlation table and do FDR correction of p values
res=multitest.multipletests(pvsflux_correlation_df.pvalue,method='fdr_tsbky')
corrected_pvsflux_correlation_df=pvsflux_correlation_df.iloc[res[0],:]
pvsflux_correlation_df['pcorr']=res[1]
pvsflux_correlation_df['FDR<0.05']=res[0]
today_mark = date.today().strftime("%Y%m%d")
pvsflux_correlation_df.to_csv('corrected_pvsflux_correlation'+today_mark+'.csv')
significant_correlated_df = pvsflux_correlation_df[pvsflux_correlation_df['FDR<0.05']]
pos_correlated_df = significant_correlated_df[significant_correlated_df['coef']>0]
neg_correlated_df = significant_correlated_df[significant_correlated_df['coef']<0]
pos_cor_val=pos_correlated_df.coef.values
neg_cor_val=neg_correlated_df.coef.values