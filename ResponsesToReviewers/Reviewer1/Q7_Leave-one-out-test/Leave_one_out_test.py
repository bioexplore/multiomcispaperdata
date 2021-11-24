# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 20:49:47 2021

这个脚本将执行与R脚本plotProteomeAllocationOfeLIFE11Categories.R相同的功能，用于回答NC审稿人1的Q7
这里用到两个文件：
This script is designed to do leave-one-out test for each category, the logic behind the code 
is as follows:
    for each protein in a specific category, do calculate all other protein fraction sum,
    then do spearman correlation analysis between this leave-one-out fraction sum and the total sum,
    if the correlation is statistical significant then we believe the category based conclusion does
    not influenced by any protein in this category.
There are two depended files:
1. eLIFE-ProteinCategoriesModify20200527.xlsx
2. ProteinAbsConcentration20210924.xlsx

@author: bioex
"""
import pandas as pd
import os
from scipy import stats #This package is imported for the spearmanr correlation analysis for leave-one-out test
from datetime import date # For making date mark for the output file

#Read in data and categories
category_file = "eLIFE-ProteinCategoriesModify20200527.xlsx"
protemo_file = "ProteinAbsConcentration20210924.xlsx"
today_mark = date.today().strftime("%Y%m%d")
leve_one_out_test_result_file = "Leave_One_Out_test_result"+today_mark+".xlsx"

# Helper leave-one-out function
def loo(xlst):
    return [[elem for elem in xlst if elem!=xlst[i]] for i in range(len(xlst))]

category_df_dict = {}
category_name_dict = {}
categories = ['Ribosomes','Translation','Glycolysis','Mitocondria','AA','Energy','Lipids','Chaperons','Nucleotides','Stress','ER','Transcription','Other']

for category in categories:
    category_df_dict[category]=pd.read_excel(category_file,sheet_name=category)
    category_name_dict[category]=category_df_dict[category]['ORF'].to_list()
    
wanted_columns = ['SystemID'] + ['w_S'+str(i) for i in range(1,28)]
proteome_weight_raw_data = pd.read_excel(protemo_file,sheet_name='unit_fmol',usecols=wanted_columns,index_col=0)
proteome_weight_drop_all_na_data = proteome_weight_raw_data.dropna(how='all')#eliminate rows that have its all values to be na

total_weight_of_found_proteins=proteome_weight_drop_all_na_data.sum(axis=0,skipna=True).values # Sum up all proteins for each sample
frac_of_each_protein_for_each_sample = proteome_weight_drop_all_na_data/total_weight_of_found_proteins # Calculate fraction of each protein in the total found protein weight

proteome_name_index = proteome_weight_drop_all_na_data.index

# Do leave-one-out test for each category
found_proteins_in_each_category_dict = {} # Dictionary storing protein name list in each categories, keyword are category names
found_protein_occupy_frac_df_in_each_category_dict = {} # Dictionary storing fractions of each protein in each category, keywords are category names
category_frac_mean_df_dict = {} # Dictionary, storing means for each experiment (triplicate sample mean, here nan is omited for doing the mean)
category_frac_mean_leave_one_out_df_dict = {} # Dictionary, storing total fraction sums with each protein leave, total sum without any protein leave put at the bottom row of the dataframe
exp_names = ['D'+str(i) for i in range(1,10)] # Nine experiment names
samp_names = ['w_S'+str(i) for i in range(1,28)] # 27 sample names
for cat in categories:# loop over all categories
    found_proteins_in_each_category_dict[cat] =  [prot for prot in proteome_name_index if prot in category_name_dict[cat]] # find all proteins in current category represent by cat
    found_protein_occupy_frac_df_in_each_category_dict[cat] = frac_of_each_protein_for_each_sample.loc[found_proteins_in_each_category_dict[cat]] # locate the protein fraction data to the current category represent by cat
    i = 0
    for exp in exp_names:#sum up triplicate samples for each experiment
        found_protein_occupy_frac_df_in_each_category_dict[cat][exp]=found_protein_occupy_frac_df_in_each_category_dict[cat][samp_names[i*3:i*3+3]].mean(axis=1)# sum up horizentally, so use axis = 1
        i+=1
    category_frac_mean_df_dict[cat] = found_protein_occupy_frac_df_in_each_category_dict[cat][exp_names]
    category_frac_mean_leave_one_out_df_dict[cat] = category_frac_mean_df_dict[cat].copy(deep=True) # get a deep copy, means the new dataframe has its own memory loc
    prots_in_current_category = category_frac_mean_df_dict[cat].index
    prots_leave_one_out_in_current_category = loo(prots_in_current_category) # Get a list of protein list that with one protein leave in each list
    for prot,loo_lst in zip(prots_in_current_category,prots_leave_one_out_in_current_category):
        category_frac_mean_leave_one_out_df_dict[cat].loc[prot] = category_frac_mean_df_dict[cat].loc[loo_lst].sum() # Get sum of other proteins fractions
    category_frac_mean_leave_one_out_df_dict[cat].loc['total'] = category_frac_mean_df_dict[cat].sum() # Get sum of all proteins in this category
    
    # Do spearman correlation
    curr_proc_df = category_frac_mean_leave_one_out_df_dict[cat].T # For convenience, do transpose 
    prots = curr_proc_df.columns
    total_sum_vec = curr_proc_df['total']
    coefs = []
    pvals = []
    for p in prots:
        curr_p_leave_one_out_vec = curr_proc_df[p]
        coef,pval = stats.spearmanr(curr_p_leave_one_out_vec,total_sum_vec,nan_policy='omit')
        coefs.append(coef)
        pvals.append(pval)
    
    category_frac_mean_leave_one_out_df_dict[cat]['coef'] = coefs # add one new column represent the spearman correlation coefficients
    category_frac_mean_leave_one_out_df_dict[cat]['pval'] = pvals # and the corresponding p values

# Output all leave-one-out test for each category to an excel file with sheets for individual category
with pd.ExcelWriter(leve_one_out_test_result_file) as writer:
    for cat in categories:
        category_frac_mean_leave_one_out_df_dict[cat].to_excel(writer,sheet_name=cat)
