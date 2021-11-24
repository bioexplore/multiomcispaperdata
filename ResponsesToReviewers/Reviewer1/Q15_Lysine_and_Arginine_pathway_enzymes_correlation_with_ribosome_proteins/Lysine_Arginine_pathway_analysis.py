# -*- coding: utf-8 -*-
"""
Created on Sun Oct  3 08:06:48 2021

This script is designed to analysis the correlation of Lysine/Arginine biosynthesis
pathway proteins to that of Ribosome proeins.

@author: bioex
"""

import pandas as pd
import os
import numpy as np
from scipy import stats
from datetime import date
from statsmodels.stats import multitest

protein_data = pd.read_excel("pvsm_new.xlsx",sheet_name='Protome',index_col=0)
protein_data = protein_data.iloc[:,0:9]

# Readin ribosome protein list
ribosome_protein_df = pd.read_excel("eLIFE-ProteinCategoriesModify20200527.xlsx",sheet_name='Ribosomes')
ribosome_genes = ribosome_protein_df.ORF# Get all ribosome genes
ribosome_proteins_detected = [protein for protein in ribosome_genes if protein in protein_data.index]# filter out genes whose protein did not detected in the proteom dataset

ribosome_proteome_data_df = protein_data.loc[ribosome_proteins_detected]# Get ribosome proteome data
ribosome_proteome_data_df.dropna(inplace=True)# Drop all protein rows with NAs
ribosome_proteome_data_column_form = ribosome_proteome_data_df.T # Transpose the matrix
ribosome_proteome_data_normalized = ribosome_proteome_data_column_form.div(ribosome_proteome_data_column_form.sum(),axis=1)

lysine_biosyn_pathway_genes_df = pd.read_excel("Lysine_Arginine_biosynthesis_genes.xlsx",sheet_name='Lysine')
arginine_biosyn_pathway_genes_df = pd.read_excel("Lysine_Arginine_biosynthesis_genes.xlsx",sheet_name='Arginine')

lysine_genes = lysine_biosyn_pathway_genes_df.Gene
arginine_genes = arginine_biosyn_pathway_genes_df.Gene

# Filter out genes that not detected in the proteome dataset
lysine_proteins_detected = [protein for protein in lysine_genes if protein in protein_data.index]
lysine_proteins_not_detected = [protein for protein in lysine_genes if protein not in protein_data.index]
print("Protein(s) not detected in the proteome data: ",lysine_proteins_not_detected)

arginine_proteins_detected = [protein for protein in arginine_genes if protein in protein_data.index]
arginine_proteins_not_detected = [protein for protein in arginine_genes if protein not in protein_data.index]
print("Protein(s) not detected in the proteome data: ",arginine_proteins_not_detected)


lysine_biosynthesis_pathway_proteome_data_df = protein_data.loc[lysine_proteins_detected]
arginine_biosynthesis_pathway_proteome_data_df = protein_data.loc[arginine_proteins_detected]

lysine_proteome_data_column_form = lysine_biosynthesis_pathway_proteome_data_df.T
arginine_proteome_data_column_form = arginine_biosynthesis_pathway_proteome_data_df.T

# Do normalization of each protein with sum of each protein across all dilution rates
lysine_proteome_data_normalized = lysine_proteome_data_column_form.div(lysine_proteome_data_column_form.sum(),axis=1)
arginine_proteome_data_normalized = arginine_proteome_data_column_form.div(arginine_proteome_data_column_form.sum(),axis=1)

# Do Spearman rank correlation analysis for lysine/arginine biosynthesis enzymes with all ribosome proteins
spearman_coef_prot_lysine_vs_ribosome = [] # Store all pair correlatoin coefficient
protein_pair_lysine_vs_ribosome =[] # corresponding lysine-biosynthesis enzyme vs ribosome protein pair
pval_prot_lysine_vs_ribosome = [] # Corresponding pvalue for each protein pair
spearman_coef_prot_arginine_vs_ribosome = [] # Store all pair correlatoin coefficient
protein_pair_arginine_vs_ribosome =[] # corresponding lysine-biosynthesis enzyme vs ribosome protein pair
pval_prot_arginine_vs_ribosome = [] # Corresponding pvalue for each protein pair
for r_prot in ribosome_proteome_data_normalized.columns:
    if sum(~ribosome_proteome_data_normalized[r_prot].isna())<3:
        continue
    for  l_prot in lysine_proteins_detected:
        # Guarantee it is possible to carry out spearman correlation
        if sum(~lysine_proteome_data_normalized[l_prot].isna())<3:
            continue
        corr,p=stats.spearmanr(lysine_proteome_data_normalized[l_prot],ribosome_proteome_data_normalized[r_prot],nan_policy='omit')
        spearman_coef_prot_lysine_vs_ribosome.append(corr)
        pval_prot_lysine_vs_ribosome.append(p)
        protein_pair_lysine_vs_ribosome.append(l_prot+'_'+r_prot)
    for  a_prot in arginine_proteins_detected:
        if sum(~arginine_proteome_data_normalized[a_prot].isna())<3:
            continue
        corr,p=stats.spearmanr(arginine_proteome_data_normalized[a_prot],ribosome_proteome_data_normalized[r_prot],nan_policy='omit')
        spearman_coef_prot_arginine_vs_ribosome.append(corr)
        pval_prot_arginine_vs_ribosome.append(p)
        protein_pair_arginine_vs_ribosome.append(a_prot+'_'+r_prot)

lysine_spearman_correlation_df = pd.DataFrame({'Protein_pair':protein_pair_lysine_vs_ribosome,
                                               'Spearman_r':spearman_coef_prot_lysine_vs_ribosome,
                                               'Pvalue':pval_prot_lysine_vs_ribosome})
arginine_spearman_correlation_df = pd.DataFrame({'Protein_pair':protein_pair_arginine_vs_ribosome,
                                               'Spearman_r':spearman_coef_prot_arginine_vs_ribosome,
                                               'Pvalue':pval_prot_arginine_vs_ribosome})

# Do FDR correction for multiple test pvalues
res = multitest.multipletests(lysine_spearman_correlation_df.Pvalue,method='fdr_tsbky')
corrected_lysine_correlation_df=lysine_spearman_correlation_df.iloc[res[0],:]
lysine_spearman_correlation_df['pcorr']=res[1]
lysine_spearman_correlation_df['FDR<0.05']=res[0]
today_mark = date.today().strftime("%Y%m%d")
lysine_spearman_correlation_df.to_csv('corrected_lysine_vs_ribosome_spearman_correlation'+today_mark+'.csv')
significant_correlated_df = lysine_spearman_correlation_df[lysine_spearman_correlation_df['FDR<0.05']]
pos_correlated_df = significant_correlated_df[significant_correlated_df['Spearman_r']>0]
neg_correlated_df = significant_correlated_df[significant_correlated_df['Spearman_r']<0]
pos_cor_val=pos_correlated_df.Spearman_r.values
neg_cor_val=neg_correlated_df.Spearman_r.values
print("Siginficant positive correlations:%d\nSignificant negative correlations:%d"%(len(pos_cor_val),len(neg_cor_val)))

res = multitest.multipletests(arginine_spearman_correlation_df.Pvalue,method='fdr_tsbky')
corrected_arginine_correlation_df=arginine_spearman_correlation_df.iloc[res[0],:]
arginine_spearman_correlation_df['pcorr']=res[1]
arginine_spearman_correlation_df['FDR<0.05']=res[0]
today_mark = date.today().strftime("%Y%m%d")
arginine_spearman_correlation_df.to_csv('corrected_arginine_vs_ribosome_spearman_correlation'+today_mark+'.csv')
significant_correlated_df = arginine_spearman_correlation_df[arginine_spearman_correlation_df['FDR<0.05']]
pos_correlated_df = significant_correlated_df[significant_correlated_df['Spearman_r']>0]
neg_correlated_df = significant_correlated_df[significant_correlated_df['Spearman_r']<0]
pos_cor_val=pos_correlated_df.Spearman_r.values
neg_cor_val=neg_correlated_df.Spearman_r.values
print("Siginficant positive correlations:%d\nSignificant negative correlations:%d"%(len(pos_cor_val),len(neg_cor_val)))

# hist plot of the correlation coefficient distribution between protein/mRNA pair
import matplotlib.pyplot as plt
import seaborn as sns

plt.figure()
#Do the hist plot
ax=sns.histplot(lysine_spearman_correlation_df,bins=20,x="Spearman_r",hue="FDR<0.05",kde=True)#,color='black',hist_kws={'linewidth':2,'edgecolor':'black','facecolor':'white'})

#Setting filled density curves
for l in ax.lines:
    x1 = l.get_xydata()[:,0]
    y1 = l.get_xydata()[:,1]
    ax.fill_between(x1,y1,color=l.get_color(),alpha=0.4)
#Set x axis limit
ax.set_xlim(-1,1.1)

#Prepare for annotation in figure
ratio_accepted = sum(lysine_spearman_correlation_df["FDR<0.05"])/len(lysine_spearman_correlation_df["FDR<0.05"])
ratio_rejected = 1 - ratio_accepted
str_accepted = 'Area ratio: %0.2f'%ratio_accepted
str_rejected = 'Area ratio: %0.2f'%ratio_rejected
#Do annotation for the two density distribution curves
ax.annotate(str_accepted,(0.85,40),(-150,50),textcoords='offset points',color='red',fontsize=12,arrowprops=dict(color='black'))
ax.annotate(str_rejected,(0.5,20),(-140,50),textcoords='offset points',color='blue',fontsize=12,arrowprops=dict(color='black'))

#Formatting the figure
ax.set_xlabel('Spearman correlation coefficient',fontsize=12)
ax.set_ylabel('Count',fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax.xaxis.set_tick_params(width=2)
ax.yaxis.set_tick_params(width=2)
ax.spines.top.set_visible(False)
ax.spines.right.set_visible(False)
ax.spines.bottom.set_linewidth(2)
ax.spines.left.set_linewidth(2)
ax.plot(1.01,0,">k",transform=ax.get_yaxis_transform(),clip_on=False)
ax.plot(-1,1.01,"^k",transform=ax.get_xaxis_transform(),clip_on=False)
plt.savefig('Spearman_correlation_dist_lysine_vs_ribosme'+today_mark+'.png',dpi=300,format='png')

plt.figure()
ax2=sns.histplot(arginine_spearman_correlation_df,bins=20,x="Spearman_r",hue="FDR<0.05",kde=True)#,color='black',hist_kws={'linewidth':2,'edgecolor':'black','facecolor':'white'})

#Setting filled density curves
for l in ax2.lines:
    x1 = l.get_xydata()[:,0]
    y1 = l.get_xydata()[:,1]
    ax2.fill_between(x1,y1,color=l.get_color(),alpha=0.4)
#Set x axis limit
ax2.set_xlim(-1,1.1)
#Prepare for annotation in figure
ratio_accepted = sum(arginine_spearman_correlation_df["FDR<0.05"])/len(arginine_spearman_correlation_df["FDR<0.05"])
ratio_rejected = 1 - ratio_accepted
str_accepted = 'Area ratio: %0.2f'%ratio_accepted
str_rejected = 'Area ratio: %0.2f'%ratio_rejected
#Do annotation for the two density distribution curves
ax2.annotate(str_accepted,(0.85,70),(-120,100),textcoords='offset points',color='red',fontsize=12,arrowprops=dict(color='black'))
ax2.annotate(str_rejected,(0.15,40),(-120,100),textcoords='offset points',color='blue',fontsize=12,arrowprops=dict(color='black'))

#Formatting the figure
ax2.set_xlabel('Spearman correlation coefficient',fontsize=12)
ax2.set_ylabel('Count',fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax2.xaxis.set_tick_params(width=2)
ax2.yaxis.set_tick_params(width=2)
ax2.spines.top.set_visible(False)
ax2.spines.right.set_visible(False)
ax2.spines.bottom.set_linewidth(2)
ax2.spines.left.set_linewidth(2)
ax2.plot(1.01,0,">k",transform=ax.get_yaxis_transform(),clip_on=False)
ax2.plot(-1,1,"^k",transform=ax.get_xaxis_transform(),clip_on=False)

plt.savefig('Spearman_correlation_dist_arginine_vs_ribosme'+today_mark+'.png',dpi=300,format='png')