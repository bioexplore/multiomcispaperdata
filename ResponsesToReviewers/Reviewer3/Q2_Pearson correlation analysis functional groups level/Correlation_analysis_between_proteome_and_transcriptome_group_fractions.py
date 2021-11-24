# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 14:49:06 2021

@author: bioex
"""
import os
import pandas as pd
from scipy import stats # For doing the spearman correlation analysis
from statsmodels.stats import multitest # For doing the FDR correction of the p-value

# Read in dataset for protein group fractoin and mRNA group fraction seperately
proteome_group_fraction_df = pd.read_excel('Group_comparison_between_proeome_and_transcriptome.xlsx',sheet_name='Proteins')
transcriptome_group_fraciton_df = pd.read_excel('Group_comparison_between_proeome_and_transcriptome.xlsx',sheet_name='mRNAs')

group_names = proteome_group_fraction_df.columns

names=[]
corrcoefs=[]
pvalues=[]
for name in group_names:
    mRNA_group_fracs = transcriptome_group_fraciton_df[name]
    protein_group_fracs = proteome_group_fraction_df[name]
    corr,p=stats.pearsonr(mRNA_group_fracs,protein_group_fracs)
    names.append(name)
    corrcoefs.append(corr)
    pvalues.append(p)
pvsm_correlation_df = pd.DataFrame({'Name':names,'coef':corrcoefs,'pvalue':pvalues})

#Get the FDR corrected p-values with a critical FDR = 0.05
res=multitest.multipletests(pvsm_correlation_df.pvalue,method='fdr_tsbky')
corrected_pvsm_correlation_df=pvsm_correlation_df.iloc[res[0],:]
pvsm_correlation_df['pcorr']=res[1]
pvsm_correlation_df['FDR<0.05']=res[0]
#form a datestamp
from datetime import date
today = date.today()
today_str=today.strftime('%Y%m%d')
pvsm_correlation_df.to_csv('corrected_pvsm_group_fraction_correlation_new'+today_str+'.csv')
