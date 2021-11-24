# -*- coding: utf-8 -*-
"""
Created on Sun Oct 03 12:31:24 2021
This script use the biopython package to read and do statistic study on the 
yeast ribosome protein sequence analysis.
@author: Jianye Xia
example code: download from https://www.biostars.org/p/62874/
from Bio import SeqIO
from Bio.SeqUtils import ProtParam

handle = open("example.fasta") 
for record in SeqIO.parse(handle, "fasta"): 
    seq = str(record.seq)
    X = ProtParam.ProteinAnalysis(seq)
    print X.count_amino_acids() 
    print X.get_amino_acids_percent() 
    print X.molecular_weight() 
    print X.aromaticity() 
    print X.instability_index() 
    print X.flexibility() 
    print X.isoelectric_point() 
    print X.secondary_structure_fraction()
"""
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import gzip
import os

import pandas as pd
from datetime import date

fast_file='sgd.fasta.gz'# fast file contain all protein amino acid sequences
sgd2swissprot_map_file = "sgdSystemId2SWISSPROT.xlsx"
seq_file_handler=gzip.open(fast_file,'rt')# use gzip to open the library
sgd2swissprot_map_df = pd.read_excel(sgd2swissprot_map_file,index_col=1)# Readin the systemid to swiss prot id map file
# usage: sgd2swissprot_map_df.loc['P23255',0] will return 'YCR042C'
ribosome_protein_df = pd.read_excel("eLIFE-ProteinCategoriesModify20200527.xlsx",sheet_name='Ribosomes',index_col=0)
ribosome_genes = ribosome_protein_df.index# Get all ribosome genes

pname_list=[]# protein name list
seq_list=[]# protein sequence list
stat_df_list=[] # statistic dataframe list
aa_count_dict_list=[]# aa count for each protein
stats_df = pd.DataFrame()
for record in SeqIO.parse(seq_file_handler,'fasta'):
    pname = record.name.split('|')
    assert len(pname)>=2
    if pname[1] in sgd2swissprot_map_df.index:
        #print(pname[1])
        if len(sgd2swissprot_map_df.loc[pname[1]].shape)>1:
            print(pname[1],': more than one sgdname mapped!')
            sgdname = sgd2swissprot_map_df.loc[pname[1]].ORF[0]
            #print('sgdname: ',sgdname)
        else:
            sgdname = sgd2swissprot_map_df.loc[pname[1]][0]
            #print('sgdname: ',sgdname)
        
        if sgdname in ribosome_genes:
            pname_list.append(sgdname)
            seq_list.append(str(record.seq))
            prot_analyze_stats = ProteinAnalysis(str(record.seq))
            aa_percent_dict = prot_analyze_stats.get_amino_acids_percent()
            aa_count_dict = prot_analyze_stats.count_amino_acids()
            aa_count_dict_list.append(aa_count_dict)
            aa_stats_df = pd.DataFrame.from_dict(aa_percent_dict,orient='index',columns=[sgdname])
            stat_df_list.append(aa_stats_df)
            stats_df=pd.concat([stats_df,aa_stats_df],axis=1)
    else:
        print(pname[1], 'is not found in the name database!')
today_mark = date.today().strftime("%Y%m%d")
stats_df.to_excel("ribosome_amino_acid_statistics"+today_mark+'.xlsx')