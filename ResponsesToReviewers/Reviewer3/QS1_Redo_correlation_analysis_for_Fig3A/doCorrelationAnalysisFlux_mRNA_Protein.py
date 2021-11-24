# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 15:16:57 2020
This script is designed for dealing with the correlation analysis between reaction
flux and mRNA level or protein level, only catabolic enzyme coding genes that carry
reaction flux get involved in the analysis. 

These catabolic enzyme coding genes are classified into three categories: one reaction
vs one protein (unique gene), one reaction vs several protein (ioenzymes), and one 
reaction vs multiple proteins (enzyme complex), based on the gene-protein-reaction 
rules defined in the genome scale metabolic model Yeast ver 7.6. 

Correlation coefficient rules for these three type of reactions:
For reactions with only one unique enzyme:
    Correlation analysis is carried out between reaction fluxes and protein levels 
    or mRNA levels of the unique enzyme.
For reactions with multiple isoenzymes:
    Correlation analysis is carried out between reaction fluxes and each ioenzyme, but
    the maximum correaltion coefficient is selected out.
For reactions with enzyme complex:
    Correlation analysis is carried out between reaction flux and complex component
    protein, which is detected and with the minimum level compared to other components,
    the same rule is applied for mNRA.
    
Modified on 2021/10/08
Codes doing FDR correction of p-values were added 


@author: bioex
"""
import pandas as pd
import re
import os
import numpy as np
from scipy import stats # for pearson r calculation
from collections import OrderedDict
import itertools
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats import multitest # For FDR corrections of p-values

# Help function for binary string to nums
def booleanarray2num(array):
    '''
    Convert a boolean array to an integer number.

    Parameters
    ----------
    array : np.array<boolean>
        Boolean np array, e.g. np.array([True,True,False,True]).

    Returns
    -------
    int 
        Return the corresponding integer number of the binary array.

    '''
    intarray = array.astype(int)
    return intarray.dot(2**np.arange(len(intarray))[::-1])
    
def r_gcomplex(gprlist):
    '''
    This function will generate the protein complex list that encoded in the GEM model.
    :param gprlist: the coded gene-protein-reaction list coded in the model
    :return:gene_complex as a list
    '''
    gene_complex=[]
    for gpr in gprlist:
        #print(gpr)
        
        # CONDITION 1: no gene assigned
        if not isinstance(gpr,str):
            gene_complex.append(np.nan)
            continue
        
        isozymes=re.findall(r' or ', gpr, flags=re.IGNORECASE)#find all the ' or ' substring in gpr
        isozyme_num=len(isozymes)+1
        
        # CONDITION 2: One enzyme or One complex
        if (isozyme_num==1):#only one enzyme or enzyme complex
            complex_num=len(re.findall(r' and ', gpr, flags=re.IGNORECASE))+1
            if(complex_num==1):#only one enzyme
                gene_complex.append([gpr])
            else:
                gpr_complex=re.findall(r'[YQ][A-Z0-9]+[WC]?-?[A-Z]{0,1}',gpr)# get all genes
                gene_complex.append(gpr_complex)
        
        # CONDITION 3: Multi- isoenzymes
        else:#more isozymes
            complex_str_list=re.findall(r'\([^(^)]+\)', gpr)# find all () pairs 
                                                            #Explaination: \( match the begining of (
                                                            #              \) means the end of )
                                                            #              [^(^) means this match donot contain any ( or )
            complex_num=len(complex_str_list) # See how many complexes we find
            # CASE 1: No complexes
            if(complex_num==0):
                gpr_complex=re.findall(r'[YQ][A-Z0-9]+[WC]?-?[A-Z]{0,1}',gpr)
                gpr_complex_list = [[comp] for comp in gpr_complex]
                gene_complex.append(gpr_complex_list)
            
            # CASE 2: More than 1 complex
            else:
                
                # C1:All isoenzymes are complex 
                if(isozyme_num==complex_num): # All isozymes are complex (contains more than one protein)
                    gpr_complex = [re.findall(r'[YQ][A-Z0-9]+[WC]?-?[A-Z]{0,1}',g) for g in complex_str_list]
                
                # C2:Some isoenzymes are note complex
                #   e.g. [['YDR353W', 'YGR209C'], ['YDR353W', 'YLR043C'], 'YDR353W']
                elif(isozyme_num>complex_num):
                    
                    if (complex_num==1):
                        isozyme_str_list = gpr.split(sep='or')
                        sub_complex_list = [re.findall(r'[YQ][A-Z0-9]+[WC]?-?[A-Z]{0,1}',g) for g in isozyme_str_list]
                        gpr_complex = sub_complex_list
                    else:
                        sub_complex_list = [re.findall(r'[YQ][A-Z0-9]+[WC]?-?[A-Z]{0,1}',g) for g in complex_str_list]
                        nocomplex_str = re.sub(r'\([^(^)]+\)','', gpr)
                        other_isozyme = re.findall(r'[YQ][A-Z0-9]+[WC]?-?[A-Z]{0,1}',nocomplex_str)
#                        if len(other_isozyme) == 1:
#                            other_isozyme_list = other_isozyme
#                        else:
                        other_isozyme_list = [[comp] for comp in other_isozyme]
                        gpr_complex = sub_complex_list + other_isozyme_list
                gene_complex.append(gpr_complex)
    return gene_complex

#### Read in data ####
data_file = "FluxTranscriptomeAndProteomeCorrelationData.xlsx"
with pd.ExcelFile(data_file) as workbook:
    yeast_7_6_df = pd.read_excel(workbook,'Yeast_7_6',na_values=['NA'])
    enzconv_df = pd.read_excel(workbook,'EnzymeNameConversion',na_values=['NA'])
    enzdata_df = pd.read_excel(workbook,'Enzyme',na_values=['NA'])
    mRNAdata_df = pd.read_excel(workbook,'mRNA',na_values=['NA'])
    flxdata_fba_df = pd.read_excel(workbook,'Fluxes',na_values=['NA'])

enz_dict = dict((key,value) for key,value in zip(enzconv_df.GeneName,enzconv_df.EnzymeOrderId))# Get enzyme name conversion dictionary
enz_abs_df = enzdata_df.drop(enzdata_df.columns[[1,2,3,4]],axis=1) # leave only the absolute values 

rxn_ids = yeast_7_6_df.Abbreviation
rxn_genes = yeast_7_6_df.GPR
enzyme_complex_list = r_gcomplex(rxn_genes)

# 将所有参与代谢反应的酶基因都列出来
list_all_catabolic_enzymes = []
i=0
unpredicted_count=0
unpredicted_list=[]
for complex_e in enzyme_complex_list:
    #if i>0: 
    #    break
    if isinstance(complex_e,list):
        for item in complex_e:
            if isinstance(item, list):
                list_all_catabolic_enzymes += item
            elif isinstance(item,str):
                list_all_catabolic_enzymes.append(item)
            else:
                print("Unpredicted:",item)
                if item not in unpredicted_list:
                    unpredicted_list.append(item)
    else:
        unpredicted_count+=1
    i+=1
# Remove duplicates
unique_catabolic_enzymes = []
[unique_catabolic_enzymes.append(x) for x in list_all_catabolic_enzymes if x not in unique_catabolic_enzymes]

# catabolic_enzyme_absvalue_df,catabolic_mRNA_absvalue_df
# are used to store unique_catabolic_enzymes levels under 9 dilution rates
valid_keys = [True if item in enz_dict.keys() else False for item in unique_catabolic_enzymes]
from itertools import compress
detected_protein_list = list(compress(unique_catabolic_enzymes,valid_keys))#Detected proteins in the proteome data 蛋白组中检测到的蛋白
from operator import itemgetter
# itemgetter can get dictionary values by the given key list
detected_proteinid_list = itemgetter(*detected_protein_list)(enz_dict)
catabolic_enzyme_absvalue_df = enz_abs_df[enz_abs_df.EnzymeOrderId.isin(detected_proteinid_list)]#Need to average这里的蛋白数据是每一个独立样品的，需要进行平均
catabolic_mRNA_absvalue_df = mRNAdata_df[mRNAdata_df.gene.isin(unique_catabolic_enzymes)]

catabolic_enzyme_absvalue_df1 = pd.DataFrame(columns=['gene','d1','d2','d3','d4','d5','d6','d7','d8','d9'])
i=0
for protein,pid in zip(detected_protein_list,detected_proteinid_list):
    tmp_vals = enz_abs_df[enz_abs_df.EnzymeOrderId==pid].iloc[:,1:].values.reshape((9,3))
    tmp_val_list = np.nanmean(tmp_vals,axis=1)
    data_series = pd.Series([protein]+tmp_val_list.tolist(),index=catabolic_enzyme_absvalue_df1.columns)
    catabolic_enzyme_absvalue_df1 = catabolic_enzyme_absvalue_df1.append(data_series,ignore_index=True)
    if i<2:
        print('tmp_vals:',tmp_vals,'\ntmp_val_list:',tmp_val_list,'\ndata_series:',data_series)
    i+=1

catabolic_enzymes_detected_under_9conds_df = catabolic_enzyme_absvalue_df1.dropna() #Get all catabolic enzymes detected获得所有9个稀释速率下都检测到的催化功能的蛋白数据
# Ok, there are 496 catabolic proteins which are detected under all 9 conditions, 496 = row numbers of catabolic_enzymes_detected_under_9conds_df

# Construct a dictionary, which links reaction name with enzyme complex
# enz_cmplx_dict[rxn_name] will give an enzyme complex value
# enzyme complex's format:
#   three kinds of value types:
#       Type 1: 'YCL040W' only one protein
#       Type 2: ['YCL040W','YFR053C','YGL253W'] complex that contains three proteins
#       Type 3: [['YCL040W'],['YFR053C'],['YGL253W']] three isozymes all of which can catalyze the reaction
enz_cmplx_dict = dict()
for rxn,cmplx in zip(yeast_7_6_df.Abbreviation,enzyme_complex_list):
    if not isinstance(cmplx,str) and not isinstance(cmplx,list):
        continue
    else:
        enz_cmplx_dict[rxn] = cmplx

# get flux values for each reaction
flux_vals = OrderedDict()
for rxn in enz_cmplx_dict.keys():
    flx_val = flxdata_fba_df[flxdata_fba_df.Rxn_name==rxn].drop('Rxn_name',axis=1).values.reshape((9,))
    flux_vals[rxn] = flx_val
    
rxn_num = len(enz_cmplx_dict)
l_num = 0
for rxn in enz_cmplx_dict.keys():
    cplx = enz_cmplx_dict[rxn]
    if isinstance(cplx,list):
        l_num += 1

# Check all is list
if rxn_num == l_num:
    print("check passed!")
    
# Parse the reaction-isozyme-genename triplet from enz_cmplx_dict
rxn_iso_gene_df = pd.DataFrame(columns=['rxn','isozyme','gene'])
for rxn in enz_cmplx_dict.keys():
    cplx = enz_cmplx_dict[rxn]
    iso_nu = 0

    for isozyme in cplx:
        if isinstance(isozyme,list):
            for prot in isozyme:
                row_info = pd.Series([rxn,iso_nu,prot],index=['rxn','isozyme','gene'])
                rxn_iso_gene_df = rxn_iso_gene_df.append(row_info,ignore_index=True)
            iso_nu += 1
        else:
            row_info = pd.Series([rxn,iso_nu,isozyme],index=['rxn','isozyme','gene'])
            rxn_iso_gene_df = rxn_iso_gene_df.append(row_info,ignore_index=True)

# rxn_iso_gene_vals_df is used to generate orignal data for 'Hierarchical_analysis_results.xlsx'
rxn_iso_gene_vals_df = rxn_iso_gene_df.copy()
rxn_iso_gene_vals_df = rxn_iso_gene_vals_df.reindex(columns=['rxn','isozyme','gene','e1','e2','e3','e4','e5','e6','e7','e8','e9','flx1','flx2','flx3','flx4','flx5','flx6','flx7','flx8','flx9'],fill_value = 0)
rownus = len(rxn_iso_gene_vals_df.gene)
for rowidx,gene,rxn in zip(range(rownus),rxn_iso_gene_vals_df.gene,rxn_iso_gene_vals_df.rxn):
    if gene in catabolic_enzyme_absvalue_df1.gene.values:
        rxn_iso_gene_vals_df.iloc[rowidx,3:12] = catabolic_enzyme_absvalue_df1[catabolic_enzyme_absvalue_df1.gene==gene].iloc[0,-9:].values
    rxn_iso_gene_vals_df.iloc[rowidx,12:21] = flux_vals[rxn]


# Count how many isozymes contain more than one genes
tt=rxn_iso_gene_df.groupby(['rxn','isozyme']).count().sort_values(by='gene',ascending=False)
sum(tt.gene>1)# numbers of isozymes that contain more than one gene

rxn_complex_num_df = rxn_iso_gene_df.groupby('rxn')['isozyme'].nunique()# how many complexes for each reaction
sorted_rxn_cplnu_df = rxn_complex_num_df.sort_values(ascending=False)# sorted version of reaction complex pairt dataframe
# See how many reaction catalyzed by more than one isozymes
sum(sorted_rxn_cplnu_df>1)# 601 reactions have more than one isozymes

    
mRNA_vals,protein_vals = OrderedDict(),OrderedDict() 
for gene in mRNAdata_df.gene:
    mRNA_val = mRNAdata_df[mRNAdata_df.gene == gene].drop('gene',axis=1).values.reshape((9,))
    mRNA_vals[gene] = mRNA_val
for gene in enz_dict.keys():
    idx = enz_dict[gene]
    raw_data = enzdata_df[enzdata_df.EnzymeOrderId==idx]
    raw_data = raw_data.drop(raw_data.columns[0:5],axis=1)
    raw_data = np.nanmean(raw_data.values.reshape((9,3)),axis=1)
    protein_val = np.nan_to_num(raw_data)
    protein_vals[gene] = protein_val
    
# For each reaction using an integer number to represent the reaction carries flux under which diluton rates
# For example: 24 can be converted to '000011000' which denotes under d5,d6 conditions the reaction carry flux
# The same info are parsed to store in mRNA_detectstr_dict and protein_detectstr_dict
# e.g. if mRNA_detctstr_dict['YDR012W']=24 it means mRNA of gene 'YDR012W' were detected under d5, d6
flx_detectstr_dict,mRNA_detectstr_dict,protein_detectstr_dict=OrderedDict(),OrderedDict(),OrderedDict()
unique_rxns = rxn_iso_gene_df.rxn.drop_duplicates()# Get unique reactions in the yeast 7.6 mdoel that have enzyme assigned, 2302 reactions included
unique_genes = rxn_iso_gene_df.gene.drop_duplicates()# Get unique genes in the yeast7.6 model, 909 genes included
for rxn in unique_rxns:
    flx_detectstr_dict[rxn]=booleanarray2num(flux_vals[rxn]!=0)
for gene in unique_genes:
    if gene in mRNA_vals.keys():
        mRNA_detectstr_dict[gene]=booleanarray2num(mRNA_vals[gene]!=0)
    else:
        mRNA_detectstr_dict[gene]=0
    if gene in protein_vals.keys():
        protein_detectstr_dict[gene]=booleanarray2num(protein_vals[gene]!=0)
    else:
        protein_detectstr_dict[gene]=0

rxn_iso_gene_detected_info_df = rxn_iso_gene_df.copy()
rxn_iso_gene_detected_info_df['fluxdstr']=0#Add one column as flux detected info num, as stored in flx_detectstr_dict
rxn_iso_gene_detected_info_df['mRNAdstr']=0#Add one column as mRNA detected info num, as stored in mRNA_detctstr_dict
rxn_iso_gene_detected_info_df['protdstr']=0#Add one column as protein detected info num, ast stored in protein_detectstr_dict
t_rownu = len(rxn_iso_gene_detected_info_df)
for i in range(t_rownu):
    rxn = rxn_iso_gene_detected_info_df.loc[i,'rxn']
    gene = rxn_iso_gene_detected_info_df.loc[i,'gene']
    rxn_iso_gene_detected_info_df.loc[i,'fluxdstr'] = flx_detectstr_dict[rxn]
    if gene in mRNA_detectstr_dict.keys():
        rxn_iso_gene_detected_info_df.loc[i,'mRNAdstr'] = mRNA_detectstr_dict[gene]
    if gene in protein_detectstr_dict.keys():
        rxn_iso_gene_detected_info_df.loc[i,'protdstr'] = protein_detectstr_dict[gene]

rxn_iso_gen_carry_flux_df = rxn_iso_gene_detected_info_df[rxn_iso_gene_detected_info_df.fluxdstr!=0]


# find all reactions that do carry fluxes 
reactions_carry_flux = []
reactions_fluxes_num = []
for rxn in flux_vals.keys():
    if flux_vals[rxn].mean() == 0:
        continue
    reactions_carry_flux.append(rxn)
    reactions_fluxes_num.append(np.count_nonzero(flux_vals[rxn]))


# Do correaltion analysis    
corr_coef_df = pd.DataFrame(columns=['rxn','gene','coef_prot','p_prot','coef_mRNA','p_mRNA'])
i=0;midx=0
for rxn,enzyme_complex in zip(rxn_ids,enzyme_complex_list):
    
    if isinstance(enzyme_complex,list): # this is a complex need find the smallest component
        # two kinds of list
        # ['',''] this is real complex, select the detectable minimum component
        # [[''],['']] this is isozymes, sum all component up
        flux_val = flux_vals[rxn]
        # 1. deal with single proteins
        if len(enzyme_complex)==1 and isinstance(enzyme_complex[0],str):# A single protein
            
            # Get the abundance data of the protein
            if (enzyme_complex[0] in enz_dict.keys()):
                idx = enz_dict[enzyme_complex[0]]
                raw_data = enzdata_df[enzdata_df.EnzymeOrderId==idx]
                raw_data = raw_data.drop(raw_data.columns[0:5],axis=1)
                raw_data = np.nanmean(raw_data.values.reshape((9,3)),axis=1)
                enz_data = np.nan_to_num(raw_data)
                prot_corr_coef_and_p = stats.pearsonr(flux_val,enz_data)
                coef_prot = prot_corr_coef_and_p[0]
                p_prot = prot_corr_coef_and_p[1]
            else:
                coef_prot = None
                p_prot = None
            if (enzyme_complex[0] in mRNA_vals.keys()):
                mRNA_corr_coef_and_p = stats.pearsonr(flux_val,mRNA_vals[enzyme_complex[0]])
                coef_mRNA = mRNA_corr_coef_and_p[0]
                p_mRNA = mRNA_corr_coef_and_p[1]  
            else:
                coef_mRNA = None
                p_mRNA = None
            row_values=[rxn,enzyme_complex[0],coef_prot,p_prot,coef_mRNA,p_mRNA]
            corr_coef_df.loc[i] = row_values
            i+=1
        # 2. deal with multiple proteins
        else:
            # Flatten the complex and make unique
            if isinstance(enzyme_complex[0],str):
                flatten_unique_genes = enzyme_complex
            else:
                flatten_unique_genes = np.unique(list(itertools.chain.from_iterable(enzyme_complex)))
            prot_coef_p_dict = OrderedDict()
            mRNA_coef_p_dict = OrderedDict()
            for gene in flatten_unique_genes:# deal with individual gene
                # Get the abundance of protein
                if(gene in enz_dict.keys()):
                    idx = enz_dict[gene]
                    raw_data = enzdata_df[enzdata_df.EnzymeOrderId==idx]
                    raw_data = raw_data.drop(raw_data.columns[0:5],axis=1)
                    raw_data = np.nanmean(raw_data.values.reshape((9,3)),axis=1)
                    enz_data = np.nan_to_num(raw_data)
                    prot_corr_coef_and_p = stats.pearsonr(flux_val,enz_data)
                    coef_prot = prot_corr_coef_and_p[0]
                    p_prot = prot_corr_coef_and_p[1]
                else:
                    coef_prot = None
                    p_prot = None
                if (gene in mRNA_vals.keys()):
                    mRNA_corr_coef_and_p = stats.pearsonr(flux_val,mRNA_vals[gene])
                    coef_mRNA = mRNA_corr_coef_and_p[0]
                    p_mRNA = mRNA_corr_coef_and_p[1]  
                else:
                    coef_mRNA = None
                    p_mRNA = None
                row_values=[rxn,gene,coef_prot,p_prot,coef_mRNA,p_mRNA]
                corr_coef_df.loc[i] = row_values
                i+=1            
    else:
        midx+=1
        print(midx,rxn,"No enzyme complex assigned!")

output_corr_coef_df=corr_coef_df.dropna()#不含nan的行

corr_coef_df1 = pd.DataFrame(columns=['rxn','gene','coef_prot','p_prot','coef_mRNA','p_mRNA'])#only flux-carring reactions
i=0;midx=0
for rxn,enzyme_complex in zip(rxn_ids,enzyme_complex_list):
    if(rxn not in reactions_carry_flux):
        continue
    if isinstance(enzyme_complex,list): # this is a complex need find the smallest component
        # two kinds of list
        # ['',''] this is real complex, select the detectable minimum component
        # [[''],['']] this is isozymes, sum all component up
        flux_val = flux_vals[rxn]
        # 1. deal with single proteins
        if len(enzyme_complex)==1 and isinstance(enzyme_complex[0],str):# A single protein
            
            # Get the abundance data of the protein
            if (enzyme_complex[0] in enz_dict.keys()):
                idx = enz_dict[enzyme_complex[0]]
                raw_data = enzdata_df[enzdata_df.EnzymeOrderId==idx]
                raw_data = raw_data.drop(raw_data.columns[0:5],axis=1)
                raw_data = np.nanmean(raw_data.values.reshape((9,3)),axis=1)
                enz_data = np.nan_to_num(raw_data)
                prot_corr_coef_and_p = stats.pearsonr(flux_val,enz_data)
                coef_prot = prot_corr_coef_and_p[0]
                p_prot = prot_corr_coef_and_p[1]
            else:
                coef_prot = None
                p_prot = None
            if (enzyme_complex[0] in mRNA_vals.keys()):
                mRNA_corr_coef_and_p = stats.pearsonr(flux_val,mRNA_vals[enzyme_complex[0]])
                coef_mRNA = mRNA_corr_coef_and_p[0]
                p_mRNA = mRNA_corr_coef_and_p[1]  
            else:
                coef_mRNA = None
                p_mRNA = None
            row_values=[rxn,enzyme_complex[0],coef_prot,p_prot,coef_mRNA,p_mRNA]
            corr_coef_df1.loc[i] = row_values
            i+=1
        # 2. deal with multiple proteins
        else:
            # Flatten the complex and make unique
            if isinstance(enzyme_complex[0],str):
                flatten_unique_genes = enzyme_complex
            else:
                flatten_unique_genes = np.unique(list(itertools.chain.from_iterable(enzyme_complex)))
            prot_coef_p_dict = OrderedDict()
            mRNA_coef_p_dict = OrderedDict()
            for gene in flatten_unique_genes:# deal with individual gene
                # Get the abundance of protein
                if(gene in enz_dict.keys()):
                    idx = enz_dict[gene]
                    raw_data = enzdata_df[enzdata_df.EnzymeOrderId==idx]
                    raw_data = raw_data.drop(raw_data.columns[0:5],axis=1)
                    raw_data = np.nanmean(raw_data.values.reshape((9,3)),axis=1)
                    enz_data = np.nan_to_num(raw_data)
                    prot_corr_coef_and_p = stats.pearsonr(flux_val,enz_data)
                    coef_prot = prot_corr_coef_and_p[0]
                    p_prot = prot_corr_coef_and_p[1]
                else:
                    coef_prot = None
                    p_prot = None
                if (gene in mRNA_vals.keys()):
                    mRNA_corr_coef_and_p = stats.pearsonr(flux_val,mRNA_vals[gene])
                    coef_mRNA = mRNA_corr_coef_and_p[0]
                    p_mRNA = mRNA_corr_coef_and_p[1]  
                else:
                    coef_mRNA = None
                    p_mRNA = None
                row_values=[rxn,gene,coef_prot,p_prot,coef_mRNA,p_mRNA]
                corr_coef_df1.loc[i] = row_values
                i+=1            
    else:
        midx+=1
        print(midx,rxn,"No enzyme complex assigned!")  

# Construct tidy data of the corr_coef_df1
tt = pd.melt(corr_coef_df1,id_vars=['rxn','gene'],value_vars=['coef_prot','coef_mRNA'])
tt.loc[tt.variable=='coef_prot','variable']='prot'
tt.loc[tt.variable=='coef_mRNA','variable']='mRNA'
tt.rename(columns={"variable":"type","value":"PearsonR"},inplace=True)
ttp = pd.melt(corr_coef_df1,id_vars=['rxn','gene'],value_vars=['p_prot','p_mRNA'])
ttp.loc[ttp.variable=='p_prot','variable']='prot'
ttp.loc[ttp.variable=='p_mRNA','variable']='mRNA'
ttp.rename(columns={"variable":"type","value":"p_value"},inplace=True)

tidy_corr_coef_df1 = tt.copy()
tidy_corr_coef_df1['p_value']=ttp['p_value']

#check consistent
for idx in tt.index:
    if tt.loc[idx,'rxn']==ttp.loc[idx,'rxn'] and tt.loc[idx,'gene']==ttp.loc[idx,'gene'] and tt.loc[idx,'type']==ttp.loc[idx,'type']:
        if tidy_corr_coef_df1.loc[idx,'p_value'] == ttp.loc[idx,'p_value']:
            print(idx)


# Do plot pearson correlation between mRNA and reaction flux
g=sns.catplot(x='rxn',y='coef_mRNA',data=corr_coef_df1.sort_values(by='coef_mRNA'))
plt.setp(g.ax.get_yticklabels(),rotation=30)

        
output_corr_coef_df1 = corr_coef_df1.dropna()#drop rows containing nan value不含nan的行

# Protein list with their corresponding reaction carrying flux带有流量的这些反应对应的蛋白列表
plist_carry_flux = corr_coef_df1.drop_duplicates('gene').gene.values.astype(str).tolist()
plist_detected_under_9d = catabolic_enzymes_detected_under_9conds_df.gene.values.astype(str).tolist()
# to see how many of the proteins in the above list has been detected in all 
plist_carry_flux_and_detetected_9d = [ item for item in plist_carry_flux if item in plist_detected_under_9d]

# Sorted data table排序后的数据表格
sorted_output_corr_coef_df = output_corr_coef_df.sort_values(by='coef_prot')
# Select out the maximum correlation coefficient among isozymes 选出每个反应中相关系数最大的结果
sorted_output_corr_coef_df_abs = sorted_output_corr_coef_df.copy()
sorted_output_corr_coef_df_abs['coef_prot_abs'] = sorted_output_corr_coef_df.coef_prot.abs()

selected_rxn_df = pd.DataFrame(columns=sorted_output_corr_coef_df_abs.columns)
unique_rxns = sorted_output_corr_coef_df_abs.drop_duplicates('rxn').rxn.values.astype(str).tolist()
for rxn in unique_rxns:
    get_row =sorted_output_corr_coef_df_abs[sorted_output_corr_coef_df_abs.rxn==rxn].sort_values(by='coef_prot_abs').iloc[-1]
    selected_rxn_df = selected_rxn_df.append(get_row,ignore_index=True)

sort_selected_rxn_df = selected_rxn_df.sort_values(by='coef_prot')

fig,axes =plt.subplots(1,2,figsize=(10,3))
sc=axes[0].scatter(sort_selected_rxn_df.coef_prot,sort_selected_rxn_df.rxn,c=sort_selected_rxn_df.p_prot)
#ax.scatter(sort_selected_rxn_df.rxn,sort_selected_rxn_df.coef_mRNA,c=sort_selected_rxn_df.p_mRNA)
cbar=fig.colorbar(sc)
cbar.set_label('p_val')

yticks = list(range(0,len(sort_selected_rxn_df.rxn),10))
ylabels = [sort_selected_rxn_df.rxn[x] for x in yticks]
axes[0].set_yticks(yticks)
axes[0].set_yticklabels(ylabels)
axes[0].set_xlabel('Pearson correlation coefficient')
axes[0].set_ylabel('Reaction ID')
fig.tight_layout()

# Do FDR correction for multiple test pvalues
res_prot = multitest.multipletests(sort_selected_rxn_df.p_prot,method='fdr_tsbky')
sort_selected_rxn_df['pcorr_prot']=res_prot[1]
sort_selected_rxn_df['Prot_FDR<0.05']=res_prot[0]
res_mRNA = multitest.multipletests(sort_selected_rxn_df.p_mRNA,method='fdr_tsbky')
sort_selected_rxn_df['pcorr_mRNA']=res_mRNA[1]
sort_selected_rxn_df['mRNA_FDR<0.05']=res_mRNA[0]

# Do FDR correction for multiple test pvalues
res_prot = multitest.multipletests(output_corr_coef_df1.p_prot,method='fdr_tsbky')
res_mRNA = multitest.multipletests(output_corr_coef_df1.p_mRNA,method='fdr_tsbky')
output_corr_coef_df1['pcorr_prot']=res_prot[1]
output_corr_coef_df1['Prot_FDR<0.05']=res_prot[0]
output_corr_coef_df1['pcorr_mRNA']=res_mRNA[1]
output_corr_coef_df1['mRNA_FDR<0.05']=res_mRNA[0]

# Use seaborn to plot
#sns.set(style='whitegrid')
g = sns.jointplot('coef_mRNA','coef_prot',data=sort_selected_rxn_df,
                  marginal_kws=dict(bins=15))
g.ax_joint.cla() # clear the scatter plot
plt.sca(g.ax_joint)
plt.scatter(sort_selected_rxn_df.coef_mRNA,sort_selected_rxn_df.coef_prot,c=sort_selected_rxn_df.p_prot)

sns.distplot(sort_selected_rxn_df.coef_prot,hist=True,kde=True)

sns.kdeplot(sort_selected_rxn_df.coef_prot,shade=True,cut=0)
sns.rugplot(sort_selected_rxn_df.coef_prot)
sns.jointplot(x='coef_prot',y='coef_mRNA',data=sort_selected_rxn_df,kind='kde')

sns.pairplot(data=output_corr_coef_df)

 
plt.scatter(sorted_output_corr_coef_df.rxn,sorted_output_corr_coef_df.coef_prot,c=sorted_output_corr_coef_df.p_prot)
plt.plot(list(sorted_output_corr_coef_df.coef_prot.values))

coef_p = output_corr_coef_df1.coef_prot
coef_m = output_corr_coef_df1.coef_mRNA
p_p = output_corr_coef_df1['Prot_FDR<0.05']
p_m = output_corr_coef_df1['mRNA_FDR<0.05']
significant_sym = p_p.copy()
significant_sym[:]=2
significant_sym.loc[(coef_p>0.5)&(p_p)&(coef_m>0.5)&(p_m)]=1
significant_sym[(coef_p<-0.5)&(p_p)&(coef_m<-0.5)&(p_m)]=3

import matplotlib.colors as col
import matplotlib.cm as cm
from pylab import tick_params
startcolor='#ff726f';midcolor = '#c0c0c0';endcolor='#008d90'
cmap2 = col.LinearSegmentedColormap.from_list('rgb_2',[startcolor,midcolor,endcolor])
cm.register_cmap(cmap=cmap2)


sns.set_style('ticks')
font1 = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 20,
}
fig,ax=plt.subplots(figsize=(12,10));ax.plot([-1,1],[-1,1],color='#c0c0c0',linewidth=10,alpha=0.9,zorder=0)
scatter= ax.scatter(coef_p,coef_m,c=significant_sym,cmap="rgb_2",alpha=0.5,s=100,zorder=1)
#from matplotlib.lines import Line2D
#legend_elements = [Line2D([0],[0],marker='o',linestyle='None',color='#ff726f',label='Positive'),
#                          Line2D([0],[0],marker='o',linestyle='None',color='#c0c0c0',label='NoSignificant'),
#                                 Line2D([0],[0],marker='o',linestyle='None',color='#008d90',label='Negative')]
#legend1=ax.legend(handles=legend_elements,loc='upper left',prop=font1)
legend1=ax.legend(*scatter.legend_elements(),loc='upper left',prop=font1,framealpha=0.5,markerscale=3,borderpad=0.2,handletextpad=0.05)#自动获取scatter的legend,markerscale是marker的大小，handletextpad是marker和text之间的距离
legend1.get_texts()[0].set_text('Positive')#修改legend的label
legend1.get_texts()[0].set_fontsize(20)
legend1.get_texts()[0].set_fontweight('bold')

legend1.get_texts()[1].set_text('Not significant')
legend1.get_texts()[1].set_fontsize(20)
legend1.get_texts()[1].set_fontweight('bold')

legend1.get_texts()[2].set_text('Negative')
legend1.get_texts()[2].set_fontsize(20)
legend1.get_texts()[2].set_fontweight('bold')

ax.add_artist(legend1);ax.set_xlim(-1,1);ax.set_ylim(-1,1)
ax.set_ylabel('Pearson correlation coefficient for mRNA',fontsize=20,fontweight='bold');ax.set_xlabel('Pearson correlation coefficient for protein',fontsize=20,fontweight='bold')
x=coef_p[significant_sym==3].min();w=coef_p[significant_sym==3].max()-coef_p[significant_sym==3].min();
y=coef_m[significant_sym==3].min();h=coef_m[significant_sym==3].max()-coef_m[significant_sym==3].min();
ax.add_patch(plt.Rectangle((x,y),w,h,ls='--',lw=2,ec='#008d90',fc='None',alpha=0.9))#绘制一个方框将positive的点框起来
x=coef_p[significant_sym==1].min();w=coef_p[significant_sym==1].max()-coef_p[significant_sym==1].min();
y=coef_m[significant_sym==1].min();h=coef_m[significant_sym==1].max()-coef_m[significant_sym==1].min();
ax.add_patch(plt.Rectangle((x,y),w,h,ls='--',lw=2,ec='#ff726f',fc='None',alpha=0.9))#绘制一个方框将negative的点框起来
ax.tick_params(which='major',width=2)
#设置字体大小
for label in ax.xaxis.get_ticklabels()+ax.yaxis.get_ticklabels():
    label.set_fontsize(20)
    label.set_fontweight('bold')
ax.spines['bottom'].set_linewidth(2)#Setting line width of bottom line
ax.spines['left'].set_linewidth(2)#Setting line width of left line
ax.spines['right'].set_linewidth(2)#Setting line width of right line
ax.spines['top'].set_linewidth(2)#Setting line width of top line

for tick in ax.get_xticklines():
    tick.set_linewidth(2)
for tick in ax.get_yticklines():
    tick.set_linewidth(2)
    
from mpl_toolkits.axes_grid1 import make_axes_locatable
#绘制边上的hist统计分布图
#create new axes on the right and on the top of the current axes
divider = make_axes_locatable(ax)
#below height and pad are in inches
ax_histx = divider.append_axes("top",1.2,pad=0.1,sharex=ax)
ax_histx.tick_params(which='major',width=2)
ax_histy = divider.append_axes("right",1.2,pad=0.1,sharey=ax)
ax_histy.tick_params(which='major',width=2)
#Make labels invisible
#ax_histx.set_xticks([])
ax_histx.set_yticks([])
#ax_histx.set_xticklabels('')
#ax_histx.set_yticklabels('')

ax_histy.set_xticks([])
#ax_histy.set_yticks([])
#ax_histy.set_xticklabels('')
#ax_histy.set_yticklabels('')
ax_histx.tick_params(labelbottom=False)
ax_histy.tick_params(labelleft=False)
ax_histx.spines.bottom.set_linewidth(2)
ax_histx.spines.left.set_linewidth(2)
ax_histx.spines.right.set_linewidth(2)
ax_histx.spines.top.set_linewidth(2)
ax_histy.spines.bottom.set_linewidth(2)
ax_histy.spines.left.set_linewidth(2)
ax_histy.spines.right.set_linewidth(2)
ax_histy.spines.top.set_linewidth(2)
#now determine nice limits by hand
binwidth=0.15
histx=coef_p.values
histy=coef_m.values
xymax=max(np.max(np.abs(histx)),np.max(np.abs(histy)))
lim = (int(xymax/binwidth)+1)*binwidth
bins = np.arange(-lim,lim+binwidth,binwidth)
ax_histx.hist(histx,bins=bins)
ax_histy.hist(histy,bins=bins,orientation='horizontal')

fig.tight_layout()
from datetime import date
today_mark = date.today().strftime("%Y%m%d")
plt.savefig('Correlation_plot_mRNAFlux_vs_protFlux1'+today_mark+'.png',dpi=300,format='png')

#列出positive和negtive 部分的rxn-gene对
tt3=output_corr_coef_df1.loc[significant_sym==3,['rxn','gene']]# Negative correlated part rxn-gene pairs
tt1=output_corr_coef_df1.loc[significant_sym==1,['rxn','gene']]# Positive correlated part rxn-gene pairs
tt=output_corr_coef_df1[['rxn','gene']] # all rxn-gene pairs

# show some statistics info about the data
# 1. For positively correlation pairs
len(tt1.groupby('rxn').count())# How many reactions in total, 77
len(tt1.groupby('gene').count())# how many genes in total, 86
sum(tt1.groupby('rxn').count().gene>1)#how many rxns have more than one gene positively correlated with flux, 8
sum(tt1.groupby('gene').count().rxn>1)#how many genes have more than one rxn positively correlated with its abundance, 8
tt1.groupby('rxn').count()[tt1.groupby('rxn').count().gene>1].sort_values(by='gene',ascending=False) # see what are the reactions that pair with multigenes
tt1.groupby('gene').count()[tt1.groupby('gene').count().rxn>1].sort_values(by='rxn',ascending=False)

# 2. For negatively correlation pairs
len(tt3.groupby('rxn').count())# How many reactions in total, 47
len(tt3.groupby('gene').count())# how many genes in total, 47
sum(tt3.groupby('rxn').count().gene>1)#how many rxns have more than one gene positively correlated with flux, 9
sum(tt3.groupby('gene').count().rxn>1)#how many genes have more than one rxn positively correlated with its abundance, 8
tt3.groupby('rxn').count()[tt3.groupby('rxn').count().gene>1].sort_values(by='gene',ascending=False) # see what are the reactions that pair with multigenes
tt3.groupby('gene').count()[tt3.groupby('gene').count().rxn>1].sort_values(by='rxn',ascending=False)

# 探索关于所有反应中的一些complex相关信息：
reaction_flux_info_df = pd.DataFrame({'rxn':reactions_carry_flux,'flx_num':reactions_fluxes_num}, columns=['rxn','flx_num'])
# how many flux carry 9 fluxes with it
len(reaction_flux_info_df.groupby('rxn').count())# in total 363 reactions contain flux under at least one condition
sum(reaction_flux_info_df.flx_num==9)# see how many rxns carry flux under all 9 conditions, 290
sum(reaction_flux_info_df.flx_num==8) # 10
sum(reaction_flux_info_df.flx_num==7) # 14
sum(reaction_flux_info_df.flx_num==6) # 4
sum(reaction_flux_info_df.flx_num==5) # 2
sum(reaction_flux_info_df.flx_num==4) # 4
sum(reaction_flux_info_df.flx_num==3) # 2
sum(reaction_flux_info_df.flx_num==2) # 18
sum(reaction_flux_info_df.flx_num==1) # 19
reaction_flux_info_df[reaction_flux_info_df.flx_num==8]


