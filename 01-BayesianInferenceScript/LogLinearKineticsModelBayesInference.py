# -*- coding: utf-8 -*-
"""
This script is used to do MCMC Bayes inference on the loglinear kinetics model simulation, by which multi-omics integration was carried out!
@author: Jianye Xia
@email: xiajy@tib.cas.cn
"""
#%matplotlib inline
import pandas as pd
import numpy as np
import pymc3 as pm
from scipy import optimize
import matplotlib.pyplot as plt
from collections import OrderedDict
import re

def r_gcomplex(gprlist):
    '''
    This function will generate the protein complex list that encoded in the GEM model.
    :param gprlist: the coded gene-protein-reaction list coded in the model
    :return:gene_complex as a list
    '''
    gene_complex=[]
    for gpr in gprlist:
        if (np.isnan(gpr)):
            gene_complex.append(np.nan)
            continue
        isozymes=re.findall(r' or ', gpr, flags=re.IGNORECASE)#find all the ' or ' substring in gpr
        isozyme_num=len(isozymes)+1
        if (isozyme_num==1):#only one enzyme or enzyme complex
            complex_num=len(re.findall(r' and ', gpr, flags=re.IGNORECASE))+1
            if(complex_num==1):#only one enzyme
                gene_complex.append(gpr)
            else:
                gpr_complex=re.findall(r'[YQ][A-Z0-9]+[WC]?-?[A-Z]{0,1}',gpr)# get all genes
                gene_complex.append(gpr_complex)
        else:#more isozymes
            complex_str_list=re.findall(r'\([^)]+\)', gpr)# find all () pairs
            complex_num=len(complex_str_list)
            if(complex_num==0):
                gpr_complex=re.findall(r'[YQ][A-Z0-9]+[WC]?-?[A-Z]{0,1}',gpr)
                gpr_complex_list = [[comp] for comp in gpr_complex]
                gene_complex.append(gpr_complex_list)
            else:
                if(isozyme_num==complex_num):
                    gpr_complex = [re.findall(r'[YQ][A-Z0-9]+[WC]?-?[A-Z]{0,1}',g) for g in complex_str_list]
                elif(isozyme_num>complex_num):
                    sub_complex_list = [re.findall(r'[YQ][A-Z0-9]+[WC]?-?[A-Z]{0,1}',g) for g in complex_str_list]
                    nocomplex_str = re.sub(r'\([^)]+\)','', gpr)
                    other_isozyme = re.findall(r'[YQ][A-Z0-9]+[WC]?-?[A-Z]{0,1}',nocomplex_str)
                    other_isozyme_list = [[comp] for comp in other_isozyme]
                    gpr_complex = sub_complex_list + other_isozyme_list
                gene_complex.append(gpr_complex)
    return gene_complex

def get_absval_enzyme(enzyme_complex,enzyme_data,enzyme_dict):
    '''
    This function will get the absolute protein abunadance for the given
    enzyme complex, the rule is like as follows:
        isoenzymes: all the detected isoenzyme will sum up
        complex components: the minimum abundance among all detected compoment is used
    [INPUT]
        enzyme_complex: the given enzyme complex 
            e.g. ['YKL001C','YBL002W'] means this complex has two components 
                [['KL001C'],['YBL002W']] means two isoenzymes
        enzyme_data: contain the dataframe storying the absolute quantification data
            note: because we have triplicate measurments, so both mean and std is calculated
        enzyme_dict: give the dictionary where to find the data
    [OUTPUT]
        absval_enzyme_list, [list], list of absolute abundance values for 9 dilution rates
            note: mean values list for each chemostat
        std_enzyme_list: [list], list of std measurement err for 9 dilution rates
            note: sample std err list for each chemostat
        Actually, a tuple (absval_enzyme_list,std_enzyme_list) will return
    '''
    absval_enzyme_list = np.zeros(27,dtype=float) # Make empty mean value list
    
    if isinstance(enzyme_complex,list): # this is a complex need find the smallest component
        #two kinds of list
        # ['',''] this is real complex, select the detectable minimum component
        # [[''],['']] this is isozymes, sum all component up
        
        # detect isoenzymes
        isoenzymes_identy = [isinstance(element,list) for element in enzyme_complex]
        if (sum(isoenzymes_identy)==0): #This is a real complex
            for component in enzyme_complex:
                abs_val = get_absval_enzyme(component,enzyme_data,enzyme_dict)
                for idx in range(27):
                    if not np.isnan(abs_val[idx]):
                        if absval_enzyme_list[idx]==0.0:
                            absval_enzyme_list[idx]=abs_val[idx]
                        else:
                            if abs_val[idx]!=0:
                                absval_enzyme_list[idx]=np.minimum(absval_enzyme_list[idx],abs_val[idx])
            return absval_enzyme_list
        else: #This is isoenzymes
            for component in enzyme_complex:
                abs_val = get_absval_enzyme(component,enzyme_data,enzyme_dict)
                absval_enzyme_list += abs_val
            return absval_enzyme_list
    elif isinstance(enzyme_complex,str):# A single protein enzyme
        if(enzyme_complex in enzyme_dict.keys()):
            idx = enzyme_dict[enzyme_complex]
            raw_data = enzyme_data[enzyme_data.EnzymeOrderId==idx]
            raw_data = raw_data.drop(raw_data.columns[0],axis=1)
            raw_data = raw_data.values.reshape((27,))
            absval_enzyme_list = np.array([a if np.isnan(m) else m+a for m,a in zip(raw_data,absval_enzyme_list)]) #list comprehension syntax [[cond true value,cond flase value][condition] for elem in list]
            return absval_enzyme_list
    else:
        print("Some error happened!")
       
# Set relative directory of data file and read it to memory
data_file = "./CChemostatsBayesModelData-withmacro.xlsm"
with pd.ExcelFile(data_file) as workbook:
    model_df = pd.read_excel(workbook,'Model',na_values=['NA'])
    metconv_df = pd.read_excel(workbook,'MetaboliteNameConversion',na_values=['NA'])
    enzconv_df = pd.read_excel(workbook,'EnzymeNameConversion',na_values=['NA'])
    metdata_df = pd.read_excel(workbook,'Metabolites',na_values=['NA'])
    enzdata_df = pd.read_excel(workbook,'Enzyme',na_values=['NA'])
    flxdata_fva_df = pd.read_excel(workbook,'FVA_Fluxes',na_values=['NA'])
    flxdata_fba_df = pd.read_excel(workbook,'Fluxes',na_values=['NA'])

met_dict = dict((key,value) for key,value in zip(metconv_df.MetId,metconv_df.MetOrderId))# Get metname conversion dictionary
enz_dict = dict((key,value) for key,value in zip(enzconv_df.GeneName,enzconv_df.EnzymeOrderId))# Get enzyme name conversion dictionary

enz_abs_df = enzdata_df.drop(enzdata_df.columns[[1,2,3,4]],axis=1) # leave only the absolute values 
met_rel_df = metdata_df.drop(metdata_df.columns[np.arange(1,13)],axis=1).replace(np.nan,0) # leave only the quantification values

# Parse the model in the reaction order
rxn_ids = model_df.Rxn_name
rxns = model_df.Formula
rxn_genes = model_df.Genes
rxn_pathway = model_df.Pathway

# Step 1. split the reactions formula to get substrates and products
split_rxns=[r.split() for r in rxns]#split all reactions
split_rxns_without_plus=[[value for value in r if value!='+'] for r in split_rxns]#eliminate plus symbol
reactants=[r[0:r.index('<=>')] for r in split_rxns_without_plus] #get reactants
products=[r[r.index('<=>')+1:len(r)] for r in split_rxns_without_plus]#get products

# Step 2. split the genes to get enzyme complexes
enzyme_complex_list = r_gcomplex(rxn_genes)
# Fill in the absolute proteome data for each enzyme complex
# For isoenzymes sum up all detected proteins value, while enzyme complex select the nonzero miminum component's abs value
# Enzyme unit is [molecules/cell], can be converted to unit [micromol/gDW] by dividing the value with 13 [pgDW/cell]= 13e-12 [gDW/cell], then divide by Avogadro constant 6.02e23 [molecules/mol], then multiply by 1e6 change mol to micromol, in a whole divide by (1.3*6.02*1e6) to change to unit of [umol/gDW]
enzyme_vals = []# enzyme values under all 27 chemostats
for elment in enzyme_complex_list:
    enzyme_val =get_absval_enzyme(enzyme_complex=elment, enzyme_data=enz_abs_df,enzyme_dict=enz_dict) 
    enzyme_vals.append(enzyme_val)
enzyme_vals_mean = [np.nanmean(element.reshape((9,3)),axis=1)/(1.3*6.02e6) for element in enzyme_vals] # Calculate average , in unit of [umol/gDW]
enzyme_vals_std = [np.nanstd(element.reshape((9,3)),axis=1,ddof=1)/(1.3*6.02e6) for element in enzyme_vals] # unit of [umol/gDW]

# Step 3. Get flux value for each testing reaction
# flux unit is [mmol/gDW/h]
flux_vals_up = []# up range of the flux calculated by FVA
flux_vals_low = []# lower range of the flux calculated by FVA
flux_vals = [] # FBA determined flux 
for rec in rxn_ids:
    flx_val_low = flxdata_fva_df[flxdata_fva_df.Rxn_name==rec].drop('Rxn_name',axis=1).values[:,::2].reshape((9,)) # odd columns corresponding to lower limit of flux
    flx_val_up = flxdata_fva_df[flxdata_fva_df.Rxn_name==rec].drop('Rxn_name',axis=1).values[:,1::2].reshape((9,)) # even columns corresponding to upper limit of flux
    flx_val = flxdata_fba_df[flxdata_fba_df.Rxn_name==rec].drop('Rxn_name',axis=1).values.reshape((9,))
    flux_vals_up.append(flx_val_up)
    flux_vals_low.append(flx_val_low)
    flux_vals.append(flx_val)
# Split the flux to up and low limit

# Step 4. Calculate the specific flux as respect to enzyme amount
specific_flux_vals_up = [flx/enz for flx,enz in zip(flux_vals_up,enzyme_vals_mean)]
specific_flux_vals_low = [flx/enz for flx,enz in zip(flux_vals_low,enzyme_vals_mean)]
specific_flux_vals = [flx/enz for flx,enz in zip(flux_vals,enzyme_vals_mean)]


# Step 5. Get metabolites values for each reaction's substrates and products
reactant_vals = []
product_vals = []
# Get normalized data. Since all data of metabolome here are raw MS ion intensity, MS intensity ratio between different metabolites do not means the same as their abundance
# Here we first normalize for each metabolites their raw MS intensity with its mean value over all conditions, after that the normalized intensity signals were then normalized 
# as respect to biomass, Here totalProteinAmount means total protein amount measured in each sample, this was done by Metabolon co ltd., totalProtenPercent means the percentage of
# total protein in mass in biomass for each dilution rate. normalization_constant has the unit of [gProt/gProt/gDW]=[gDW]  
totalProteinAmount = metdata_df[metdata_df.BIOCHEMICAL=='PROTEIN'].drop(metdata_df.columns[np.arange(13)],axis=1).values.reshape((26,))# Because the first sample S1 was contaminated, so we only have 2 samples for the first dilution rate (D=0.025h-1)
totalProteinPercent = metdata_df[metdata_df.BIOCHEMICAL=='TOTAL_PROTEIN_CONTENT'].drop(metdata_df.columns[np.arange(13)],axis=1).values.reshape((26,))
normalization_constant = totalProteinAmount/totalProteinPercent
ridx=0
for racs,prds in zip(reactants,products):
    subs = np.zeros((len(racs),26),dtype=float)# Each substrate occupy one row,for first dilution rate there are only 2 samples 
    prods = np.zeros((len(prds),26),dtype=float)# Each product occupy one row
    i = 0
    for rac in racs:
        if(rac in met_dict.keys()):
            idx = met_dict[rac]   # found the metabolite's correspnding idx
            if (idx != 0): # idx == 0 means this metabolite was not detected in our metabolome dataset
                raw_data = metdata_df[metdata_df.MetOderId==idx] # get the raw MS ion intensity values for the metabolite
                raw_data = raw_data.drop(raw_data.columns[np.arange(13)],axis=1).values.reshape((26,))#drop the first 13 columns
                if idx!=421:#omit ATP
                    # First do internal normalization
                    raw_data_int_norm = raw_data/np.nanmedian(raw_data) 
                    # then do normalization to biomass
                    raw_data_norm = raw_data_int_norm/normalization_constant # After this step, the unit of metabolite should be [a.u./gDW]
                    subs[i,:] = raw_data_norm
                else:
                    subs[i,:] = raw_data# ATP data is filled manually, which is already normalized to biomass, and has the unit of [mmol/gDW]
            else:# print out the not detected metabolites
                MetInfo = metconv_df[metconv_df.MetId==rac].MetName
                
                print('undetected substrate:',MetInfo.values[0],'reaction:',model_df.Rxn_name[ridx])
        i += 1
    j = 0
    for prd in prds:
        if(prd in met_dict.keys()):
            idx = met_dict[prd]   # found the key
            if (idx != 0):
                raw_data = metdata_df[metdata_df.MetOderId==idx]
                raw_data = raw_data.drop(raw_data.columns[np.arange(13)],axis=1).values.reshape((26,))
                if idx!=421:#omit ATP
                    # First do internal normalization
                    raw_data_int_norm = raw_data/np.nanmedian(raw_data)
                    # then do normalization to biomass
                    raw_data_norm = raw_data_int_norm/normalization_constant
                    prods[j,:] = raw_data_norm
                else:
                    prods[j,:] = raw_data
            else:
                MetInfo = metconv_df[metconv_df.MetId==prd].MetName
                print('undetected product:',MetInfo.values[0],'reaction:',model_df.Rxn_name[ridx])
        j += 1
    reactant_vals.append(subs)
    product_vals.append(prods)
    ridx += 1

reactant_vals_mean = []# Store the reactants relative concentration values
reactant_vals_std = []# Standard error
product_vals_mean = []# Store the products relative concentration values
product_vals_std = []# Standard error
# Get avergae and std of metabolite values for 9 chemostat
for reactants,products in zip(reactant_vals,product_vals):
    sub_nu = reactants.shape[0] # Get how many substrates take part in this reaction
    prd_nu = products.shape[0]
    exp_nu = (reactants.shape[1]+1)/3 # Condition numbers with triplicate samples for each dilution rate except the first dilution rate has duplicate samples
    sub_nul_col = np.array([np.nan]*sub_nu).reshape((sub_nu,1))# Add one null column to simulate the S1, so that make the reshape be easier
    prd_nul_col = np.array([np.nan]*prd_nu).reshape((prd_nu,1))
    
    tt_subs = np.hstack((sub_nul_col,reactants)).flatten().reshape((-1,3))
    tt_subs_mean = np.nanmean(tt_subs,axis=1).reshape((sub_nu,int(exp_nu)))
    tt_subs_std = np.nanstd(tt_subs,axis=1).reshape((sub_nu,int(exp_nu)))
    reactant_vals_mean.append(tt_subs_mean)
    reactant_vals_std.append(tt_subs_std)
    
    tt_prds = np.hstack((prd_nul_col,products)).flatten().reshape((-1,3))
    tt_prds_mean = np.nanmean(tt_prds,axis=1).reshape((prd_nu,int(exp_nu)))
    tt_prds_std = np.nanstd(tt_prds,axis=1).reshape((prd_nu,int(exp_nu)))
    product_vals_mean.append(tt_prds_mean)
    product_vals_std.append(tt_prds_std)


# Define the log likelihood function for being used in the bayes inference
def logp(ux,lx,mu,sigma):
    cdf_up = pm.math.exp(pm.Normal.dist(mu,sigma).logcdf(ux))
    cdf_low = pm.math.exp(pm.Normal.dist(mu,sigma).logcdf(lx))
    return pm.math.log(cdf_up-cdf_low)-pm.math.log(ux-lx)

# Define a recursive function to seek for the best model that has the lowest LOO value
def findMinBestModel(optionalNM,rltv_met_vals_mean,j0,j_obs_up,j_obs_low,j_obs_val,mets_list,best_nm=None, lowerModel=None):
    sim_type = 0 # simulation type 0: use average kcat as mu, average std as sigma
                 #                 1: use 1 as mean and 10 as sigma
    
    if lowerModel is None:#it means this is the first round
        models,traces,loos=OrderedDict(),OrderedDict(),OrderedDict()
        compareDict,nameConvDict = dict(),dict()
        j_obs_mean = j_obs_val.mean()
        j_obs_std = j_obs_val.std()
        met_idx=0
        for nm in optionalNM:
            try:
                with pm.Model() as models[nm]:
                    mu_par_coef  = np.log(rltv_met_vals_mean[met_idx]).max()
                    if sim_type == 0 :
                        if(mu_par_coef!=0):
                            a = pm.Normal(nm,mu=j_obs_mean/mu_par_coef,sigma=j_obs_std)
                        else:
                            a = pm.Normal(nm,mu=1,sigma=10)
                    else:
                        a = pm.Normal(nm,mu=1,sigma=10)
                    j_prd = pm.Deterministic('j_prd',j0 + a.dot(pm.math.log(rltv_met_vals_mean[met_idx])))
                    j_std = pm.Deterministic('j_std',pm.math.sqrt(((j_prd-j_obs_val)**2).sum()/(8-1)))
                    normal_dist = pm.Normal.dist(j_prd,j_std)
                    j_obs = pm.DensityDist('j_obs',logp,observed={'ux':j_obs_up,'lx':j_obs_low,'mu':j_prd,'sigma':j_std},random=normal_dist.random)
                    traces[nm] = pm.sample(cores=4,start=pm.find_MAP(fmin=optimize.fmin_powell),progressbar=True)
            except RuntimeError:
                print('\033[1;31;43mCatch error!\033[0m, move to next metabolite!')
                continue
            loos[nm] = pm.loo(traces[nm],models[nm])
            kvals = loos[nm].pareto_k.values
            if len(kvals) > 2* sum(kvals>0.7):
                compareDict[models[nm]] = traces[nm]
                nameConvDict[models[nm]] = nm
            met_idx += 1
        if not compareDict:
            print('\033[1;31;43m Warning \033[0m: No suitable model found for reaction',rxn_idx,'!')
            return None,None
        compRst = pm.compare(compareDict)
        print('Now among\033[1;31;43m',optionalNM,'\033[0mvalid models are\033[1;31;43m',list(nameConvDict.values()),'\033[0m\n\tNow do comparison between these models!')
        print(compRst)
        best_md_loc = compRst.index[compRst['rank']==0][0]
        best_nm.append(nameConvDict[best_md_loc])
        best_tc_loc = traces[nameConvDict[best_md_loc]]
        best_md = (best_md_loc,best_tc_loc)
        print('Until now the best model is\033[1;31;43m',best_nm,'\033[0m')
        optionalNM.remove(nameConvDict[best_md_loc])
        return findMinBestModel(optionalNM,rltv_met_vals_mean,j0,j_obs_up,j_obs_low,j_obs_val,mets_list,best_nm,best_md)
    else:
        assert best_nm
        models,traces,loos=OrderedDict(),OrderedDict(),OrderedDict()
        compareDict,nameConvDict = dict(),dict()
        j_obs_mean = j_obs_val.mean()
        j_obs_std = j_obs_val.std()
        for nm in optionalNM:
            to_be_evaluated_list = [int(mets_list.index(m)) for m in best_nm]            
            to_be_evaluated_list.append(int(mets_list.index(nm)))
            r_m_v_m_matrix = rltv_met_vals_mean[to_be_evaluated_list]
            para_nums = len(to_be_evaluated_list)
            evaluate_log_mets = np.log(r_m_v_m_matrix).mean(axis=1).max()
            mu_par_coef = para_nums * evaluate_log_mets
            print('to be evaluated mets:\033[1;31;43m',best_nm+[nm],'\033[0m')
            with pm.Model() as models[nm]:
                if sim_type == 0:
                    a = pm.Normal('a', mu=j_obs_mean/mu_par_coef, sigma=j_obs_std, shape=len(to_be_evaluated_list))
                else:
                    a = pm.Normal('a', mu=1, sigma=10, shape=len(to_be_evaluated_list))
                j_prd = pm.Deterministic('j_prd', j0 + a.dot(pm.math.log(r_m_v_m_matrix)))
                j_std = pm.Deterministic('j_std', pm.math.sqrt(((j_prd-j_obs_val)**2).sum()/(8-len(to_be_evaluated_list))))
                normal_dist = pm.Normal.dist(j_prd,j_std)
                j_obs = pm.DensityDist('j_obs',logp,observed={'ux':j_obs_up,'lx':j_obs_low,'mu':j_prd,'sigma':j_std},random=normal_dist.random)
                traces[nm] = pm.sample(cores=4,start=pm.find_MAP(fmin=optimize.fmin_powell),progressbar=True)
            loos[nm] = pm.loo(traces[nm],models[nm])
            kvals = loos[nm].pareto_k.values
            if len(kvals) > 2* sum(kvals>0.7):
                compareDict[models[nm]] = traces[nm]
                nameConvDict[models[nm]] = nm
        compareDict[lowerModel[0]] = lowerModel[1]
        assert compareDict
        compRst = pm.compare(compareDict)
        print('Except\033[1;31;43m',best_nm,' \033[0mother valid optional combinations will be among\033[1;31;43m',list(nameConvDict.values()),'\033[0m')
        print(compRst)
 
        best_md_loc = compRst.index[compRst['rank']==0][0]
        if best_md_loc == lowerModel[0]:
            print('Finally, found the best model is\033[1;31;43m',best_nm,'\033[0m')
            return best_nm,lowerModel
        else:
            best_tc_loc = traces[nameConvDict[best_md_loc]]
            best_md = (best_md_loc,best_tc_loc)
            best_nm.append(nameConvDict[best_md_loc])
            print('Until now the best model is\033[1;31;43m',best_nm,'\033[0m')
            optionalNM.remove(nameConvDict[best_md_loc])
            return findMinBestModel(optionalNM,rltv_met_vals_mean,j0,j_obs_up,j_obs_low,j_obs_val,mets_list,best_nm,best_md)

# Define a recursive function to seek for the best model that has the lowest LOO value
def findLevelBestModel(rltv_met_vals_mean,j0,j_obs_up,j_obs_low,j_obs_val,condidates,mets_list,level):
    sim_type = 0 # simulation type 0: use average kcat as mu, average std as sigma
                 #                 1: use 1 as mean and 10 as sigma
    # level means how many metabolites was eliminated from the participates
    # level=0 means all n take part in metabolites
    # level=1 means all n-1 participated metabolites
    # level=n-1 means all individual metabolite
    
    # First, get all candiate model names
    from itertools import combinations
    met_nums = len(condidates)
    condidate_list = list(combinations(condidates,met_nums-level))
    models,traces,loos = OrderedDict(),OrderedDict(),OrderedDict() # dict to save the model, trace and loo
    compareDict,nameConvDict = dict(),dict() # dictionary used to save valid models/traces and name convertion dict
    # Output some information
    print('\tNow we are in level:\033[1;31;43m',level,'\033[0m, there are \033[1;31;43m',met_nums-level,'\033[0m pars!')
    
    for condidate in condidate_list:#loop over all condidate
        # Output some information
        print('\t\tNow processing the condidate:\033[1;31;43m',condidate,'\033[0m!')
        cd_nm = '-'.join(condidate)#Name the condidate
        to_be_evaluated_list = [int(mets_list.index(m)) for m in condidate]
        # Get parameter estimation for the prioor
        j_obs_mean = j_obs_val.mean()
        j_obs_std = j_obs_val.std()
        
        par_nums = met_nums - level # Determine how many parameters are there in the model
        if par_nums > 1:
            max_val = np.log(rltv_met_vals_mean[to_be_evaluated_list]).mean(axis=1).max()# Get the log value of relative metabolite value
            min_val = np.log(rltv_met_vals_mean[to_be_evaluated_list]).mean(axis=1).min()
            max_log_rltv_met = max_val if np.abs(max_val)>np.abs(min_val) else min_val
            if not max_log_rltv_met :
                max_log_rltv_met = 1
        elif par_nums == 1:
            max_log_rltv_met = np.log(rltv_met_vals_mean[to_be_evaluated_list]).mean().max()
            if not max_log_rltv_met:
                max_log_rltv_met = 1
        if sim_type == 0:
            mu_prior = j_obs_mean/par_nums/max_log_rltv_met
            sigma_prior = j_obs_std
        else:
            mu_prior = 1
            sigma_prior = 10
        
        #DEBUG INFO
        if max_log_rltv_met == 0:
            print('Got it!')
        print('par_nums:',par_nums,'max_log_rltv_met:',max_log_rltv_met,'mu_prior:',mu_prior)
        
        # Do the Bayes inference for current condidate
        with pm.Model() as models[cd_nm]:
            assert mu_prior,sigma_prior
            a = pm.Normal('a',mu=mu_prior,sigma=sigma_prior,shape=par_nums)
            j_prd = pm.Deterministic('j_prd', j0 + a.dot(pm.math.log(rltv_met_vals_mean[to_be_evaluated_list])))
            j_std = pm.Deterministic('j_std', pm.math.sqrt(((j_prd-j_obs_val)**2).sum()/(8-par_nums)))
            normal_dist = pm.Normal.dist(j_prd,j_std)
            j_obs = pm.DensityDist('j_obs',logp,observed={'ux':j_obs_up,'lx':j_obs_low,'mu':j_prd,'sigma':j_std},random=normal_dist.random)
            traces[cd_nm] = pm.sample(draws=20000,tune=20000,cores=4,start=pm.find_MAP(fmin=optimize.fmin_powell),progressbar=True)
        loos[cd_nm] = pm.loo(traces[cd_nm],models[cd_nm])
        #DEBUG INFO
        print('loo[',cd_nm,']:\t',loos[cd_nm])
        kvals = loos[cd_nm].pareto_k.values # Extract k values from the loo result
        if len(kvals) > 2*sum(kvals>0.7): # drop all models that its kvalues > 0.7 numbers larger than half of the parameters
            compareDict[models[cd_nm]] = traces[cd_nm] # Save valid model and trace for candidate metabolite list cd_nm
            nameConvDict[models[cd_nm]] = cd_nm # Get memory of the model corresponding candiate metalist 
    # After finish the loop, need get the best 
    if not compareDict:
        return None,None,None
    compRst = pm.compare(compareDict) # Do comparison between all models
    # Show the whole table
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(compRst)
        
    best_level_model = compRst.index[compRst['rank']==0][0] # Find model in the 1st rank 
    best_level_model_name = nameConvDict[best_level_model]
    best_level_trace = traces[best_level_model_name]
    
    print('Caution: find the best model for level\033[1;31;43m',level,'\033[0m is \033[1;31;43m',best_level_model_name,'\033[0m')
    return best_level_model,best_level_trace,best_level_model_name

def findBestModel(rltv_met_vals_mean,j0,j_obs_up,j_obs_low,j_obs_val,candidates,mets_list,rxn_nm):
    compareDict,nameConvDict = dict(),dict()
    met_nums = len(candidates)
    for level in range(met_nums):
        mod,trc,nm = findLevelBestModel(rltv_met_vals_mean,j0,j_obs_up,j_obs_low,j_obs_val,candidates,mets_list,level)
        if not mod:
            continue
        compareDict[mod] = trc
        nameConvDict[mod] = nm
    #DEBUG
    print('\nWill do compare for reaction:\033[1;31;43m',rxn_nm,'\033[0m\n')
    if not compareDict: # if compareDict is empty
        return None,None,None
    compRst = pm.compare(compareDict) # Do model comparison
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        print(compRst)
    best_model = compRst.index[compRst['rank']==0][0]
    best_trace = compareDict[best_model]
    best_model_nm = nameConvDict[best_model]
    print('In function: \033[1;31;43mfindBestModel\033[0m, we find the best model is \033[1;31;43m',best_model_nm,'\033[0m')
    return best_model_nm,best_model,best_trace

# Get ready for the MCMC Bayes inference
# The loglinear kinetics is as follows:
# j = j_0 + sum_i(a_i*sub_i_c/sub_i_0)
# subscript 0 is the reference state, subscript i represent the ith substrate
# subscript j repreents the jth product, subscript c means c condition, subscript 0 means reference state
# New strategy for model construction:
#   1. Add only one model coefficient, do sampling, after converged compare and choose the best fit one
#   2. Start from the above determined intrinsic turn over time coefficient, add one extra metabolite into the model, and using likelihood ratio test to determine whether the fit get better
#       if the new model performs better (through kappa-square test), update the model structure to include the corresponding coefficient, if not stop the construction
#   3. Do step 2 for one extra coefficient, until no better model structure is found
model_trace,model_coeff_nm = OrderedDict(),OrderedDict()
for rxn_idx in range(len(rxn_ids)): #choose one reaction
    print('\nNow processing reaction:\033[1;31;43m',rxn_idx,'\033[0mreaction name:\033[1;31;43m',rxn_ids[rxn_idx],'\033[0m')
    # Get all needed infos for doing sampling from this reaction
    # j_obs_up              : Up limit for enzymatic specific flux, vector, contains 8 values (corresponding to 8 dilution rates a.r.t. reference dilution rate)
    # j_obs_low             : Lower limit for enzymatic specific flux, vector, contains 8 values
    # j_obs_val             : FBA determined fix point enzymatic specific flux, vector, contains 8 values
    # j0                    : Enzymetic specific flux under the reference state, here we choose D = 0.1 h-1
    # rltv_subs_vals_mean   : Relative abundance of substrates as respect to reference condition, i.e. (X^c/X^0), vector, contains 8 values
    # rltv_prds_vals_mean   : Relative abundance of products a.r.t. reference condtion, i.e. (X^c/X^0), vector, contains 8 values
    
    j0 = specific_flux_vals[rxn_idx][1]
    j_obs_up = np.delete(specific_flux_vals_up[rxn_idx],1)
    j_obs_low = np.delete(specific_flux_vals_low[rxn_idx],1)
    j_obs_val = np.delete(specific_flux_vals[rxn_idx],1)
    if(j_obs_val.sum()==0):
        cprint('\nWarning: all flux of this reaction if zero! Escape it!','blue','on_yellow')
        continue
    #elif(rxn_ids[rxn_idx] in ['r_0984','r_1021']):
    #    print('\033[1;31;43mWarning\033[0m: Reaction r_0984 was skipped as its all metabolites were not deteced')
    #    continue
                                        
    # In the following block, all columns of the matrix of substrate values wad divided by the 2nd column and remove the 2nd column after that
    # e.g. tt = np.arange(15).reshape((3,5))
    # tt: 
    #array([[ 0,  1,  2,  3,  4],
    #       [ 5,  6,  7,  8,  9],
    #       [10, 11, 12, 13, 14]])
    # divide each column by the 2nd column will give, (tt.T/tt[:,1]).T
    #array([[ 0/1,  1/1,  2/1,  3/1,  4/1],
    #       [ 5/6,  6/6,  7/6,  8/6,  9/6],
    #       [10/11, 11/11, 12/11, 13/11, 14/11]])
    # which equals
    #array([[0.        , 1.        , 2.        , 3.        , 4.        ],
    #       [0.83333333, 1.        , 1.16666667, 1.33333333, 1.5       ],
    #       [0.90909091, 1.        , 1.09090909, 1.18181818, 1.27272727]])
    # remove the 2nd column, np.delete((tt.T/tt[:,1]).T,1,axis=1)
    #array([[0.        , 2.        , 3.        , 4.        ],
    #       [0.83333333, 1.16666667, 1.33333333, 1.5       ],
    #       [0.90909091, 1.09090909, 1.18181818, 1.27272727]])
    rltv_subs_vals_mean = reactant_vals_mean[rxn_idx].T/reactant_vals_mean[rxn_idx][:,1] # This is a matrix, e.g. two substrates, it will be 2*9, divide all other columns a.r.t. the reference column
    rltv_subs_vals_mean = rltv_subs_vals_mean.T
    rltv_subs_vals_mean = np.delete(rltv_subs_vals_mean,1,axis=1)
    
    rltv_prds_vals_mean = product_vals_mean[rxn_idx].T/product_vals_mean[rxn_idx][:,1] # This is a matrix, e.g. two substrates, it will be 2*9, divide all other columns a.r.t. the reference column
    rltv_prds_vals_mean = rltv_prds_vals_mean.T
    rltv_prds_vals_mean = np.delete(rltv_prds_vals_mean,1,axis=1)
    
    # vertical stack the substrates and products data together
    rltv_met_vals_mean = np.vstack((rltv_subs_vals_mean,rltv_prds_vals_mean))
    rltv_met_vals_mean = np.nan_to_num(rltv_met_vals_mean,nan=1) # Replace np.nan to 1
    
    # Do the best model searching
    met_nums = rltv_met_vals_mean.shape[0]
    optional_nm = ['m'+str(i) for i in range(1,met_nums+1) if rltv_met_vals_mean[i-1].mean()!=1.0] # Skip not detected metabolites
    mets_list = ['m'+str(i) for i in range(1,met_nums+1)]
    if len(optional_nm) == 0:
        continue
    best_nm = []
    # model_trace       Store the best model and corresponding trace
    # model_coeff_nm    Store the best model's corresponding metabolites in the best model
#    mt, nm = findBestModel(optional_nm,rltv_met_vals_mean,j0,j_obs_up,j_obs_low,j_obs_val,mets_list,best_nm)   
#    model_trace[rxn_ids[rxn_idx]]=mt
#    model_coeff_nm[rxn_ids[rxn_idx]]=nm
    nm,md,mt = findBestModel(rltv_met_vals_mean,j0,j_obs_up,j_obs_low,j_obs_val,optional_nm,mets_list,rxn_ids[rxn_idx])
    model_trace[rxn_ids[rxn_idx]]=(md,mt)
    model_coeff_nm[rxn_ids[rxn_idx]]=nm
# Do forest plot of all reaction's inferenced intrinsic turn over numbers
# The forest plot contains point with the mean value for each parameter (corresponding to a metabolite)
# together a credible-interval denoted by a err-bar like style will be plotted 
# The format of the plot likes
#------------------------------------------------------------------------
#                   |                                                   |
#           MetName1|    ----------O--------                            |
#ReactionNo         |                                                   |
#           MetName2|         ---------O------                          |
#                   |                                                   |
#-----------------------------------------------------------------------|
#                   |                                                   |
#ReactNo    MetName1|                         -----O--------            |
#                   |                                                   |
#-----------------------------------------------------------------------|
#......             |                                                   |
#......             |                                                   |
#------------------------------------------------------------------------
#                       0.2     0.4     0.8     1.0     1.2     1.4
#==========================================================================
means,pmin,pmax,vnames=OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()
for react in rxn_ids:
    if not react in model_trace.keys():
        continue
    tc = model_trace[react]
    #if react in ['r_1021','r_0982','r_0091','r_0886']:#these reactions have large error
        #continue
    if not tc:
        continue
    elif not tc[1]:
        continue
    tc = tc[1]    
    pos_p_samples = tc.get_values(varname=tc.varnames[0])
    #p_nums = pos_p_samples.shape[1]
    vnames[react] = model_coeff_nm[react]
    means[react] = pos_p_samples.mean(axis=0)
    pos_p_hpds = pm.stats.hpd(pos_p_samples,0.95)
    if pos_p_samples.ndim ==1:
        pmin[react] = pos_p_hpds[0]
        pmax[react] = pos_p_hpds[1]
    else:
        pmin[react] = pos_p_hpds[:,0]
        pmax[react] = pos_p_hpds[:,1]

# Make a dataframe for plotting
#%matplotlib auto
#reactions = ['r_0534','r_0467','r_0450','r_1054','r_0486','r_0892','r_0893','r_0366','r_0962','r_0889','r_1049','r_1050','r_0958','r_0300','r_0658','r_0831','r_1022','r_0450','r_0451','r_0713']
reactions = ['r_0450','r_1054','r_0486','r_0892','r_0893','r_0366','r_1049','r_1050','r_0958','r_0300','r_0658','r_0450']
plt.rcParams['figure.subplot.left']=0.15# Set left margin in ratio
plt.rcParams['figure.subplot.right']=0.98 # set right margin
plt.rcParams['figure.subplot.bottom']=0.02 # set right margin
plt.rcParams['figure.subplot.top']=0.98 # set right margin

plt.rcParams['axes.labelsize']='medium'
plt.rcParams['ytick.labelsize']='medium'
fig,axs = plt.subplots(len(reactions),figsize=(8,24), sharex=True,gridspec_kw={'hspace':0})

i=0
for react in reactions:
    axs[i].plot(means[react],vnames[react].split('-'),'o')
    axs[i].hlines(vnames[react].split('-'),pmin[react],pmax[react])
    p_nums = len(vnames[react].split('-'))
    axs[i].vlines(0,-0.5,p_nums-0.5,linestyle='dotted',color='r')
    axs[i].set_ylim(-0.5,p_nums-0.5)
    axs[i].set_ylabel(react)
    i += 1
fig.tight_layout()

# Do model prediction ability test
# In this plot, the scatter plot between posterior parameter predicted reaction specific fluxes vs the FBA fluxes
# was drawn for each reactions under all 8 conditions (without that of the reference condition)
# The format of the plot likes
# ---------------|--------------|--------------|
# |           * /|             /|             /|
# |  R2=0.8   */ | R2=0.6     / | R2=0.3     / |
# |           /  |           /  |           /  |
# |         */   |          /   |          /   |
# |         /    |         /    |         /    |
# |        /*    |        /  *  |        /     |
# |       / *    |       /      |       /     *|  j_predicted by eqn. 9
# |     */ *     |      /*  *   |      /       |
# |     /*       |     /        |     /      * |
# |    /*        |    / *       |    /    *    |
# |  */          |  */          |   /  *       |
# |  /*          |  /           | */ *         |
# | /*           | /*           | /            |
# |/             |/  *          |/ *           |
# ----------------------------------------------
#                j_observed by FBA
spec_flux_predict,spec_flux_min,spec_flux_max,spec_flux_observed,ppcs = OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict(),OrderedDict()
rxnid = 0
for react in rxn_ids: #Loop over all reactions
    if not react in model_trace.keys():
        continue # omit non fitted models
    tc = model_trace[react]
    if not tc:
        continue
    elif not tc[1]:
        continue
    md = tc[0]
    tc = tc[1]  # Get the trace  
    
    spec_flux_samples = pm.sample_posterior_predictive(tc,500,md)
    spec_flux_mean = spec_flux_samples['j_obs'].mean(axis=0)
    spec_flux_std = spec_flux_samples['j_obs'].std(axis=0)
    flux_min = spec_flux_mean - spec_flux_std
    flux_max = spec_flux_mean + spec_flux_std
    j_obs_vals = np.delete(specific_flux_vals[rxn_idx],1)
    
    spec_flux_predict[react] = spec_flux_mean
    spec_flux_min[react] = flux_min
    spec_flux_max[react] = flux_max
    spec_flux_observed[react] = j_obs_vals
    ppcs[react] = spec_flux_samples
    
    rxnid += 1

# Make category for each reaction to the corresponding pathway
reac_pathway_dict = dict(model_df[['Rxn_name','Pathway']].values.astype(str))

fig_emp,axs_emp = plt.subplots(4,3,figsize = (12,12))
fig_tca,axs_tca = plt.subplots(2,4,figsize = (16,8))
fig_hmp,axs_hmp = plt.subplots(2,3,figsize = (12,4))
emp_idx,tca_idx,hmp_idx = 0,0,0
r2_values=OrderedDict()
data_df = pd.DataFrame()
for react in rxn_ids:
    if not react in model_trace.keys():
        continue # omit non fitted models
    tc = model_trace[react]
    if not tc:
        continue
    if not react in spec_flux_predict.keys():
        continue
    data_df[react+'_prd'] = spec_flux_predict[react]
    data_df[react+'_obs'] = spec_flux_observed[react]
    
    linx = [spec_flux_observed[react].min(),spec_flux_observed[react].max()]
    liny = [spec_flux_observed[react].min(),spec_flux_observed[react].max()]
    from scipy import stats
    r2_values[react] = stats.pearsonr(spec_flux_observed[react],spec_flux_predict[react])
    if reac_pathway_dict[react] == 'EMP':
        axs_emp[emp_idx//3,np.mod(emp_idx,3)].plot(linx,liny,linestyle='dotted')
        axs_emp[emp_idx//3,np.mod(emp_idx,3)].plot(spec_flux_observed[react],spec_flux_predict[react],'o',label='$R^2$ = {:.3f}\n$p$ = {:.3f}'.format(r2_values[react][0]**2,r2_values[react][1]))
        axs_emp[emp_idx//3,np.mod(emp_idx,3)].vlines(spec_flux_observed[react],spec_flux_min[react],spec_flux_max[react])
        leg_emp = axs_emp[emp_idx//3,np.mod(emp_idx,3)].legend(loc="upper left",handlelength=0, handletextpad=0, fancybox=True)
        for item in leg_emp.legendHandles:
            item.set_visible(False)
        axs_emp[emp_idx//3,np.mod(emp_idx,3)].set_title(react)
        emp_idx += 1
    elif reac_pathway_dict[react] == 'TCA':
        axs_tca[tca_idx//4,np.mod(tca_idx,4)].plot(linx,liny,linestyle='dotted')
        axs_tca[tca_idx//4,np.mod(tca_idx,4)].plot(spec_flux_observed[react],spec_flux_predict[react],'o',label='$R^2$ = {:.3f}\n$p$ = {:.3f}'.format(r2_values[react][0]**2,r2_values[react][1]))
        axs_tca[tca_idx//4,np.mod(tca_idx,4)].vlines(spec_flux_observed[react],spec_flux_min[react],spec_flux_max[react])
        axs_tca[tca_idx//4,np.mod(tca_idx,4)].legend(loc="upper left")
        axs_tca[tca_idx//4,np.mod(tca_idx,4)].set_title(react)
        tca_idx += 1
    elif reac_pathway_dict[react] == 'PPP':
        axs_hmp[hmp_idx//3,np.mod(hmp_idx,3)].plot(linx,liny,linestyle='dotted')
        axs_hmp[hmp_idx//3,np.mod(hmp_idx,3)].plot(spec_flux_observed[react],spec_flux_predict[react],'o',label='$R^2$ = {:.3f}\n$p$ = {:.3f}'.format(r2_values[react][0]**2,r2_values[react][1]))
        axs_hmp[hmp_idx//3,np.mod(hmp_idx,3)].vlines(spec_flux_observed[react],spec_flux_min[react],spec_flux_max[react])
        axs_hmp[hmp_idx//3,np.mod(hmp_idx,3)].legend(loc="upper left")
        axs_hmp[hmp_idx//3,np.mod(hmp_idx,3)].set_title(react)
        hmp_idx += 1
fig_emp.tight_layout()
fig_tca.tight_layout()
fig_hmp.tight_layout()
