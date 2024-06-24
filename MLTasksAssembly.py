# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 21:03:11 2022

@author: tsharma2
"""

from scipy.io import loadmat
import pandas as pd
import numpy as np
import os
from itertools import combinations
import argparse

folderISAAC = "./"
DataPath = 'C:/Users/teesh/OneDrive - Indian Institute of Technology Guwahati/Dataset/EssentialData/'
if os.path.exists(DataPath)!=True:
    DataPath = 'Dataset/EssentialData/'


def tumor_types(Cancer_type):
    Map = {'GBMLGG': ['GBM', 'LGG'],
           'COADREAD': ['COAD', 'READ'],
           'KIPAN': ['KIRC', 'KICH', 'KIRP'],
           'STES': ['ESCA', 'STAD'],
           'PanGI': ['COAD', 'STAD', 'READ', 'ESCA'],
           'PanGyn': ['OV', 'CESC', 'UCS', 'UCEC'],
           'PanSCCs': ['LUSC', 'HNSC', 'ESCA', 'CESC', 'BLCA'],
           'PanPan': ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC',
                           'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG',
                           'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ',
                           'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM']
           }
    if Cancer_type not in Map:
        Map[Cancer_type] = [Cancer_type]

    return Map[Cancer_type]

def get_prognosis_count_RG(dataset,groups,genders):
    
    prog_data = dataset[7]
    prog_data = pd.DataFrame(prog_data)
    pos_source_male = genders[0]+'1'+groups[0]
    neg_source_male = genders[0]+'0'+groups[0]
    pos_target_male = genders[0]+'1'+groups[1]
    neg_target_male = genders[0]+'0'+groups[1]
    pos_source_female = genders[1]+'1'+groups[0]
    neg_source_female = genders[1]+'0'+groups[0]
    pos_target_female = genders[1]+'1'+groups[1]
    neg_target_female = genders[1]+'0'+groups[1]
    
    temp = prog_data.apply(lambda x : True
                      if x[0] == pos_source_male else False, axis = 1)
    count_MALE_SOURCE_P = len(temp[temp == True].index)
    
    temp = prog_data.apply(lambda x : True
                      if x[0] == neg_source_male else False, axis = 1)
    count_MALE_SOURCE_N = len(temp[temp == True].index)
    
    temp = prog_data.apply(lambda x : True
                      if x[0] == pos_target_male else False, axis = 1)
    count_MALE_TARGET_P = len(temp[temp == True].index)
    
    temp = prog_data.apply(lambda x : True
                      if x[0] == neg_target_male else False, axis = 1)
    count_MALE_TARGET_N = len(temp[temp == True].index)
    
    temp = prog_data.apply(lambda x : True
                      if x[0] == pos_source_female else False, axis = 1)
    count_FEMALE_SOURCE_P = len(temp[temp == True].index)
    
    temp = prog_data.apply(lambda x : True
                      if x[0] == neg_source_female else False, axis = 1)
    count_FEMALE_SOURCE_N = len(temp[temp == True].index)
    
    temp = prog_data.apply(lambda x : True
                      if x[0] == pos_target_female else False, axis = 1)
    count_FEMALE_TARGET_P = len(temp[temp == True].index)
    
    temp = prog_data.apply(lambda x : True
                      if x[0] == neg_target_female else False, axis = 1)
    count_FEMALE_TARGET_N = len(temp[temp == True].index)
    
    prog_count = [count_MALE_SOURCE_P, count_MALE_SOURCE_N, count_MALE_TARGET_P, count_MALE_TARGET_N,
                  count_FEMALE_SOURCE_P, count_FEMALE_SOURCE_N, count_FEMALE_TARGET_P, count_FEMALE_TARGET_N]
    
    return prog_count

def get_prognosis_count(dataset,groups):
    
    prog_data = dataset[4]
    prog_data = pd.DataFrame(prog_data)
    pos_source = '1'+groups[0]
    neg_source = '0'+groups[0]
    pos_target = '1'+groups[1]
    neg_target = '0'+groups[1]
    
    temp = prog_data.apply(lambda x : True
                      if x[0] == pos_source else False, axis = 1)
    count_SOURCE_P = len(temp[temp == True].index)
    
    temp = prog_data.apply(lambda x : True
                      if x[0] == neg_source else False, axis = 1)
    count_SOURCE_N = len(temp[temp == True].index)
    
    temp = prog_data.apply(lambda x : True
                      if x[0] == pos_target else False, axis = 1)
    count_TARGET_P = len(temp[temp == True].index)
    
    temp = prog_data.apply(lambda x : True
                      if x[0] == neg_target else False, axis = 1)
    count_TARGET_N = len(temp[temp == True].index)
    
    prog_count = [count_SOURCE_P, count_SOURCE_N, count_TARGET_P, count_TARGET_N]
    
    return prog_count

def race_count(Data):
    
    race_info = pd.DataFrame(Data)
    
    temp = race_info.apply(lambda x : True
                      if x[0] == "WHITE" else False, axis = 1)
    count_WHITE = len(temp[temp == True].index)
    # print('The count of WHITE is: '+str(count_WHITE))
    
    temp = race_info.apply(lambda x : True
                             if x[0] == "BLACK" else False, axis = 1)
    count_BLACK = len(temp[temp == True].index)
    # print('The count of BLACK is: '+str(count_BLACK))
    
    temp = race_info.apply(lambda x : True
                             if x[0] == "ASIAN" else False, axis = 1)
    count_ASIAN = len(temp[temp == True].index)
    # print('The count of ASIAN is: '+str(count_ASIAN))
    
    temp = race_info.apply(lambda x : True
                             if x[0] == "NAT_A" else False, axis = 1)
    count_NAT_A = len(temp[temp == True].index)
    # print('The count of NAT_A is: '+str(count_NAT_A))
    
    temp = race_info.apply(lambda x : True
                             if x[0] == "OTHER" else False, axis = 1)
    count_OTHER = len(temp[temp == True].index)
    # print('The count of OTHER is: '+str(count_OTHER))
    
    #count_race_final = [count_WHITE, count_BLACK, count_ASIAN, count_NAT_A, count_OTHER]
    
    count_race_final = pd.DataFrame({"WHITE":[count_WHITE],
                                    "BLACK":[count_BLACK],
                                    "ASIAN":[count_ASIAN],
                                    "NAT_A":[count_NAT_A],
                                    "OTHER":[count_OTHER]})
    
    return count_race_final

def add_race_info_GA(GAData,groups):
    
    GAData_race = pd.concat(GAData)
    GAData_race = GAData_race[GAData_race['EIGENSTRAT'].isin(['EA', 'AA', 'EAA', 'NA', 'OA'])]
    #Unique_count = GAData_race['EIGENSTRAT'].value_counts()
    # print('The unique race values in GA data are:')
    # print(Unique_count)
    
    GAData_race['race'] = GAData_race['EIGENSTRAT']
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'EA', 'race'] = 'WHITE'
    temp = GAData_race.apply(lambda x : True
                             if x['race'] == "WHITE" else False, axis = 1)
    count_WHITE = len(temp[temp == True].index)
    # print('The count of WHITE is: '+str(count_WHITE))
    
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'AA', 'race'] = 'BLACK'
    temp = GAData_race.apply(lambda x : True
                             if x['race'] == "BLACK" else False, axis = 1)
    count_BLACK = len(temp[temp == True].index)
    # print('The count of BLACK is: '+str(count_BLACK))
    
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'EAA', 'race'] = 'ASIAN'
    temp = GAData_race.apply(lambda x : True
                             if x['race'] == "ASIAN" else False, axis = 1)
    count_ASIAN = len(temp[temp == True].index)
    # print('The count of ASIAN is: '+str(count_ASIAN))
    
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'NA', 'race'] = 'NAT_A'
    temp = GAData_race.apply(lambda x : True
                             if x['race'] == "NAT_A" else False, axis = 1)
    count_NAT_A = len(temp[temp == True].index)
    # print('The count of NAT_A is: '+str(count_NAT_A))
    
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'OA', 'race'] = 'OTHER'
    temp = GAData_race.apply(lambda x : True
                             if x['race'] == "OTHER" else False, axis = 1)
    count_OTHER = len(temp[temp == True].index)
    # print('The count of OTHER is: '+str(count_OTHER))
    
    GAData_race = GAData_race.drop(columns=['EIGENSTRAT'])
    #Unique_count = GAData_race['race'].value_counts()
    # print('The unique race values in GA data after renaming are:')
    # print(Unique_count)
    
    GAData_race = GAData_race[GAData_race['race'].isin(groups)]
    
    return GAData_race, count_WHITE, count_BLACK, count_ASIAN, count_NAT_A, count_OTHER

def add_clinical_outcomes_TCGA(OutcomeData):
    
    OutcomeData.columns = ['E', 'T']
    OutcomeData = OutcomeData[OutcomeData['E'].isin([0, 1])]
    # print('The size of clinical outcome data after filtering values other than 0 and 1 is:')
    # print(np.shape(OutcomeData))
    OutcomeData = OutcomeData.dropna()
    # print('The size of clinical outcome data after filtering rows with missing values is:')
    # print(np.shape(OutcomeData))
    OutcomeData['C'] = 1 - OutcomeData['E']
    OutcomeData.drop(columns=['E'], inplace=True)
    
    return OutcomeData

def keep_patients_race_info(InputData,GAData_race,OutcomeData):
    
    InputData = InputData.join(GAData_race, how='inner')
    InputData = InputData.dropna(axis='columns')
    InputData = InputData.join(OutcomeData, how='inner')
    
    return InputData

def get_n_years(dataset, years):
    
    X, T, C, E, R = dataset['X'], dataset['T'], dataset['C'], dataset['E'], dataset['R']

    df = pd.DataFrame(X)
    df['T'] = T
    df['C'] = C
    df['R'] = R
    df['Y'] = 1

    df = df[~((df['T'] < 365 * years) & (df['C'] == 1))]
    df.loc[df['T'] <= 365 * years, 'Y'] = 0
    df['strat'] = df.apply(lambda row: str(row['Y']) + str(row['R']), axis=1)
    df = df.reset_index(drop=True)

    R = df['R'].values
    Y = df['Y'].values
    y_strat = df['strat'].values
    df = df.drop(columns=['T', 'C', 'R', 'Y', 'strat'])
    X = df.values
    y_sub = R # doese not matter

    return (X, Y.astype('int32'), R, y_sub, y_strat)

def get_n_years_gender(dataset, years):
    
    X, T, C, E, G = dataset['X'], dataset['T'], dataset['C'], dataset['E'], dataset['G']

    df = pd.DataFrame(X)
    df['T'] = T
    df['C'] = C
    df['G'] = G
    df['Y'] = 1

    df = df[~((df['T'] < 365 * years) & (df['C'] == 1))]
    df.loc[df['T'] <= 365 * years, 'Y'] = 0
    df['strat'] = df.apply(lambda row: str(row['Y']) + str(row['G']), axis=1)
    df = df.reset_index(drop=True)

    G = df['G'].values
    Y = df['Y'].values
    y_strat = df['strat'].values
    df = df.drop(columns=['T', 'C', 'G', 'Y', 'strat'])
    X = df.values
    y_sub = G # doese not matter

    return (X, Y.astype('int32'), G, y_sub, y_strat)

def get_n_years_gender_race(dataset, years):
    
    X, T, C, E, R, G = dataset['X'], dataset['T'], dataset['C'], dataset['E'], dataset['R'], dataset['G']

    df = pd.DataFrame(X)
    df['T'] = T
    df['C'] = C
    df['R'] = R
    df['G'] = G
    df['Y'] = 1

    df = df[~((df['T'] < 365 * years) & (df['C'] == 1))]
    df.loc[df['T'] <= 365 * years, 'Y'] = 0
    df['strat'] = df.apply(lambda row: str(row['Y']) + str(row['R']), axis=1)
    df['Gstrat'] = df.apply(lambda row: str(row['Y']) + str(row['G']), axis=1)
    df['GRstrat'] = df.apply(lambda row: str(row['G']) + str(row['Y']) + str(row['R']), axis=1)
    df = df.reset_index(drop=True)

    R = df['R'].values
    G = df['G'].values
    Y = df['Y'].values
    y_strat = df['strat'].values
    Gy_strat = df['Gstrat'].values
    GRy_strat = df['GRstrat'].values
    df = df.drop(columns=['T', 'C', 'R', 'G', 'Y', 'strat', 'Gstrat', 'GRstrat'])
    X = df.values
    y_sub = R # doese not matter

    return (X, Y.astype('int32'), R, y_sub, y_strat, G, Gy_strat, GRy_strat)

def add_race_info_GA_Methylation(MethyAncsData,groups):
    
    GAData_race = pd.concat(MethyAncsData)
    race_groups = ['WHITE',
              'BLACK OR AFRICAN AMERICAN',
              'ASIAN',
              'AMERICAN INDIAN OR ALASKA NATIVE',
              'NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER']
    GAData_race = GAData_race.rename(columns={'race': 'EIGENSTRAT'})
    GAData_race = GAData_race[GAData_race['EIGENSTRAT'].isin(race_groups)]
    #Unique_count = GAData_race['EIGENSTRAT'].value_counts()
    # print('The unique race values in GA data are:')
    # print(Unique_count)
    
    GAData_race['race'] = GAData_race['EIGENSTRAT']
    
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'WHITE', 'race'] = 'WHITE'
    temp = GAData_race.apply(lambda x : True
                             if x['race'] == "WHITE" else False, axis = 1)
    count_WHITE = len(temp[temp == True].index)
    # print('The count of WHITE is: '+str(count_WHITE))
    
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'BLACK OR AFRICAN AMERICAN', 'race'] = 'BLACK'
    temp = GAData_race.apply(lambda x : True
                             if x['race'] == "BLACK" else False, axis = 1)
    count_BLACK = len(temp[temp == True].index)
    # print('The count of BLACK OR AFRICAN AMERICAN is: '+str(count_BLACK))
    
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'ASIAN', 'race'] = 'ASIAN'
    temp = GAData_race.apply(lambda x : True
                             if x['race'] == "ASIAN" else False, axis = 1)
    count_ASIAN = len(temp[temp == True].index)
    # print('The count of ASIAN is: '+str(count_ASIAN))
    
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'AMERICAN INDIAN OR ALASKA NATIVE', 'race'] = 'NAT_A'
    temp = GAData_race.apply(lambda x : True
                             if x['race'] == "NAT_A" else False, axis = 1)
    count_NAT_A = len(temp[temp == True].index)
    # print('The count of AMERICAN INDIAN OR ALASKA NATIVE is: '+str(count_NAT_A))
    
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER', 'race'] = 'OTHER'
    temp = GAData_race.apply(lambda x : True
                             if x['race'] == "OTHER" else False, axis = 1)
    count_OTHER = len(temp[temp == True].index)
    # print('The count of NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER is: '+str(count_OTHER))
    
    GAData_race = GAData_race.drop(columns=['EIGENSTRAT'])
    #Unique_count = GAData_race['race'].value_counts()
    # print('The unique race values in GA data after renaming are:')
    # print(Unique_count)
    
    GAData_race = GAData_race[GAData_race['race'].isin(groups)]
    
    return GAData_race, count_WHITE, count_BLACK, count_ASIAN, count_NAT_A, count_OTHER

def readMethylationDataGenderRace(Cancer_type, Target, genders, groups, DataPath):
    
    MethylationDataPath = DataPath + 'MethylationData/Methylation.mat'
    MethylationData = loadmat(MethylationDataPath)
    
    # extracting input combinations data...
    X, Y, GeneName, SampleName = MethylationData['X'].astype('float32'), MethylationData['CancerType'], MethylationData['FeatureName'][0], MethylationData['SampleName']
    GeneName = [row[0] for row in GeneName]
    SampleName = [row[0][0] for row in SampleName]
    Y = [row[0][0] for row in Y]
    
    MethylationData_X = pd.DataFrame(X, columns=GeneName, index=SampleName)
    MethylationData_Y = pd.DataFrame(Y, index=SampleName, columns=['Disease'])
    MethylationData_Y = MethylationData_Y[MethylationData_Y['Disease'].isin(tumor_types(Cancer_type))]
    MethylationData_in = MethylationData_X.join(MethylationData_Y, how='inner')
    MethylationData_in = MethylationData_in.drop(columns=['Disease'])
    
    index = MethylationData_in.index.values
    index_new = [row[:12] for row in index]
    MethylationData_in.index = index_new
    MethylationData_in = MethylationData_in.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    
    # adding race information...
    MethyAncsDataPath = DataPath + 'MethylationData/MethylationGenetic.xlsx'
    
    # fetching race info from MethylationGenetic.xlsx
    MethyAncsData = [pd.read_excel(MethyAncsDataPath,
                         disease, usecols='A,B',
                         index_col='bcr_patient_barcode',
                         keep_default_na=False)
           for disease in tumor_types(Cancer_type)]
    
    GAData_race = pd.concat(MethyAncsData)
    race_groups = ['WHITE',
              'BLACK OR AFRICAN AMERICAN',
              'ASIAN',
              'AMERICAN INDIAN OR ALASKA NATIVE',
              'NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER']
    GAData_race = GAData_race.rename(columns={'race': 'EIGENSTRAT'})
    GAData_race = GAData_race[GAData_race['EIGENSTRAT'].isin(race_groups)]
    GAData_race['race'] = GAData_race['EIGENSTRAT']
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'WHITE', 'race'] = 'WHITE'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'BLACK OR AFRICAN AMERICAN', 'race'] = 'BLACK'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'ASIAN', 'race'] = 'ASIAN'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'AMERICAN INDIAN OR ALASKA NATIVE', 'race'] = 'NAT_A'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER', 'race'] = 'OTHER'
    GAData_race = GAData_race.drop(columns=['EIGENSTRAT'])
    GAData_race = GAData_race[GAData_race['race'].isin(groups)]
    
    # fetching outcome data from MethylationClinInfo.xlsx
    MethyCIDataPath = DataPath + 'MethylationData/MethylationClinInfo.xlsx'
    if Target=='OS':
        cols = 'A,D,Y,Z'
    elif Target == 'DSS':
        cols = 'A,D,AA,AB'
    elif Target == 'DFI': # this info is very less in methylation data
        cols = 'A,D,AC,AD'
    elif Target == 'PFI':
        cols = 'A,D,AE,AF'
    OutcomeData = pd.read_excel(MethyCIDataPath,
                                usecols=cols,dtype={'OS': np.float64},
                                index_col='bcr_patient_barcode')
    
    # adding clinical outcome endpoints data...
    OutcomeData.columns = ['G', 'E', 'T']
    OutcomeData = OutcomeData[OutcomeData['G'].isin(genders)]
    OutcomeData = OutcomeData[OutcomeData['E'].isin([0, 1])]
    OutcomeData = OutcomeData.dropna()
    OutcomeData['C'] = 1 - OutcomeData['E']
    OutcomeData.drop(columns=['E'], inplace=True)
    
    # Keep patients with race information
    MethylationData_in = MethylationData_in.join(GAData_race, how='inner')
    MethylationData_in = MethylationData_in.dropna(axis='columns')
    MethylationData_in = MethylationData_in.join(OutcomeData, how='inner')
    MethylationData_in = MethylationData_in.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    
    # data
    Data = MethylationData_in
    C = Data['C'].tolist()
    R = Data['race'].tolist()
    G = Data['G'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'race', 'T', 'G'])
    X = Data.values
    X = X.astype('float32')
    InputData = {'X': X,
                  'T': np.asarray(T, dtype=np.float32),
                  'C': np.asarray(C, dtype=np.int32),
                  'E': np.asarray(E, dtype=np.int32),
                  'R': np.asarray(R),
                  'G': np.asarray(G),
                  'Samples': Data.index.values,
                  'FeatureName': list(Data)}
    
    return InputData
    
def readMethylationData(Cancer_type, Target, groups, DataPath):
    
    MethylationDataPath = DataPath + 'MethylationData/Methylation.mat'
    MethylationData = loadmat(MethylationDataPath)
    
    # extracting input combinations data...
    X, Y, GeneName, SampleName = MethylationData['X'].astype('float32'), MethylationData['CancerType'], MethylationData['FeatureName'][0], MethylationData['SampleName']
    GeneName = [row[0] for row in GeneName]
    SampleName = [row[0][0] for row in SampleName]
    Y = [row[0][0] for row in Y]
    
    MethylationData_X = pd.DataFrame(X, columns=GeneName, index=SampleName)
    MethylationData_Y = pd.DataFrame(Y, index=SampleName, columns=['Disease'])
    MethylationData_Y = MethylationData_Y[MethylationData_Y['Disease'].isin(tumor_types(Cancer_type))]
    MethylationData_in = MethylationData_X.join(MethylationData_Y, how='inner')
    MethylationData_in = MethylationData_in.drop(columns=['Disease'])
    # print('The shape of fetched methylation data initially is:')
    # print(MethylationData_in.shape)
    #data_in = MethylationData_in.shape[0]
    
    index = MethylationData_in.index.values
    index_new = [row[:12] for row in index]
    MethylationData_in.index = index_new
    MethylationData_in = MethylationData_in.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    # print('The shape of fetched methylation data after dropping duplicates is:')
    # print(MethylationData_in.shape)
    
    # adding race information...
    MethyAncsDataPath = DataPath + 'MethylationData/MethylationGenetic.xlsx'
    
    # fetching race info from MethylationGenetic.xlsx
    MethyAncsData = [pd.read_excel(MethyAncsDataPath,
                         disease, usecols='A,B',
                         index_col='bcr_patient_barcode',
                         keep_default_na=False)
           for disease in tumor_types(Cancer_type)]
    #GAData_all = pd.concat(MethyAncsData)
    # print('The shape of genetic ancestry data is:')
    # print(GAData_all.shape)
    #GAdata_in = GAData_all.shape[0]
    
    # fetching outcome data from MethylationClinInfo.xlsx
    MethyCIDataPath = DataPath + 'MethylationData/MethylationClinInfo.xlsx'
    
    if Target=='OS':
        cols = 'A,Y,Z'
    elif Target == 'DSS':
        cols = 'A,AA,AB'
    elif Target == 'DFI': # this info is not present in methylation data
        cols = 'A,AC,AD'
    elif Target == 'PFI':
        cols = 'A,AE,AF'
    
    OutcomeData = pd.read_excel(MethyCIDataPath,
                                usecols=cols,dtype={'OS': np.float64},
                                index_col='bcr_patient_barcode')
    
    GAData_race, count_WHITE, count_BLACK, count_ASIAN, count_NAT_A, count_OTHER = add_race_info_GA_Methylation(MethyAncsData,groups)
    # print('The shape of GA data is:')
    # print(GAData_race.shape)
    count_race_GA = [count_WHITE, count_BLACK, count_ASIAN, count_NAT_A, count_OTHER]
    
    # adding clinical outcome endpoints data...
    OutcomeData = add_clinical_outcomes_TCGA(OutcomeData)
    # print('The shape of Outcome data is:')
    # print(OutcomeData.shape)
    
    # Keep patients with race information
    MethylationData_in = MethylationData_in.join(GAData_race, how='inner')
    MethylationData_in = MethylationData_in.dropna(axis='columns')
    MethylationData_in = MethylationData_in.join(OutcomeData, how='inner')
    MethylationData_in = MethylationData_in.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    
    # print('The shape of patients with race data is:')
    # print(MethylationData_in.shape)
    
    # data
    Data = MethylationData_in
    C = Data['C'].tolist()
    R = Data['race'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'race', 'T'])
    X = Data.values
    X = X.astype('float32')
    InputData = {'X': X,
                  'T': np.asarray(T, dtype=np.float32),
                  'C': np.asarray(C, dtype=np.int32),
                  'E': np.asarray(E, dtype=np.int32),
                  'R': np.asarray(R),
                  'Samples': Data.index.values,
                  'FeatureName': list(Data)}
    
    return InputData, count_race_GA

def readMicroRNAData(Cancer_type, Target, groups, DataPath):
    
    # data paths
    MicroRNADataPath = DataPath + 'MicroRNAData/MicroRNA-Expression.mat'
    GADataPath = DataPath + 'Genetic_Ancestry.xlsx'
    OutcomeDataPath = DataPath + 'TCGA-Clinical-Outcome-Endpoints.xlsx'
    #
    MicroRNAData = loadmat(MicroRNADataPath)
    X, Y, GeneName, SampleName = MicroRNAData['X'].astype('float32'), MicroRNAData['CancerType'], MicroRNAData['FeatureName'][0], MicroRNAData['SampleName']
    GeneName = [row[0] for row in GeneName]
    SampleName = [row[0][0] for row in SampleName]
    Y = [row[0][0] for row in Y]
    
    MicroRNAData_X = pd.DataFrame(X, columns=GeneName, index=SampleName)
    MicroRNAData_Y = pd.DataFrame(Y, index=SampleName, columns=['Disease'])
    MicroRNAData_Y = MicroRNAData_Y[MicroRNAData_Y['Disease'].isin(tumor_types(Cancer_type))]
    MicroRNAData_in = MicroRNAData_X.join(MicroRNAData_Y, how='inner')
    MicroRNAData_in = MicroRNAData_in.drop(columns=['Disease'])
    # print('The shape of input data before dropping duplicates is:')
    # print(MicroRNAData_in.shape)
    #data_in = MicroRNAData_in.shape[0]
    
    index = MicroRNAData_in.index.values
    index_new = [row[:12] for row in index]
    MicroRNAData_in.index = index_new
    MicroRNAData_in = MicroRNAData_in.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    # print('The shape of input data after dropping duplicate is:')
    # print(MicroRNAData_in.shape)
    #data_in_drop = MicroRNAData_in.shape[0]
    
    GAData = [pd.read_excel(GADataPath,
                         disease, usecols='A,E',
                         index_col='Patient_ID',
                         keep_default_na=False) 
           for disease in tumor_types(Cancer_type)]
    #GAData_all = pd.concat(GAData)
    # print('The shape of genetic ancestry data is:')
    # print(GAData_all.shape)
    #GAdata_in = GAData_all.shape[0]
    
    # print('fetching Clinical Outcomes data...')
    if Target=='OS':
        cols = 'B,H,I'
    elif Target == 'DSS':
        cols = 'B,J,K'
    elif Target == 'DFI':
        cols = 'B,L,M'
    elif Target == 'PFI':
        cols = 'B,N,O'
    OutcomeData = pd.read_excel(OutcomeDataPath,
                           'TCGA-CDR', usecols=cols,
                           index_col='bcr_patient_barcode')
    
    # adding race information...
    GAData_race, count_WHITE, count_BLACK, count_ASIAN, count_NAT_A, count_OTHER = add_race_info_GA(GAData, groups)
    # print('The shape of GA data is:')
    # print(GAData_race.shape)
    count_race_GA = [count_WHITE, count_BLACK, count_ASIAN, count_NAT_A, count_OTHER]
    
    # adding clinical outcome endpoints data...
    OutcomeData = add_clinical_outcomes_TCGA(OutcomeData)
    # print('The shape of Outcome data is:')
    # print(OutcomeData.shape)
    
    # Keep patients with race information
    MicroRNAData_in = keep_patients_race_info(MicroRNAData_in,GAData_race,OutcomeData)
    # print('The shape of patients with race data is:')
    # print(MicroRNAData_in.shape)
    
    # data
    Data = MicroRNAData_in
    C = Data['C'].tolist()
    R = Data['race'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'race', 'T'])
    X = Data.values
    X = X.astype('float32')
    InputData = {'X': X,
                  'T': np.asarray(T, dtype=np.float32),
                  'C': np.asarray(C, dtype=np.int32),
                  'E': np.asarray(E, dtype=np.int32),
                  'R': np.asarray(R),
                  'Samples': Data.index.values,
                  'FeatureName': list(Data)}
    
    return InputData, count_race_GA

def readmRNAData(Cancer_type, Target, groups, DataPath):
    
    # data paths
    mRNADataPath = DataPath + 'mRNAData/mRNA-Expression.mat'
    GADataPath = DataPath + 'Genetic_Ancestry.xlsx'
    OutcomeDataPath = DataPath + 'TCGA-Clinical-Outcome-Endpoints.xlsx'
    #
    mRNAData = loadmat(mRNADataPath)
    X, Y, GeneName, SampleName = mRNAData['X'].astype('float32'), mRNAData['CancerType'], mRNAData['FeatureName'][0], mRNAData['SampleName']
    GeneName = [row[0] for row in GeneName]
    SampleName = [row[0][0] for row in SampleName]
    Y = [row[0][0] for row in Y]
    
    mRNAData_X = pd.DataFrame(X, columns=GeneName, index=SampleName)
    mRNAData_Y = pd.DataFrame(Y, index=SampleName, columns=['Disease'])
    mRNAData_Y = mRNAData_Y[mRNAData_Y['Disease'].isin(tumor_types(Cancer_type))]
    mRNAData_in = mRNAData_X.join(mRNAData_Y, how='inner')
    mRNAData_in = mRNAData_in.drop(columns=['Disease'])
    # print('The shape of input data before dropping duplicates is:')
    # print(mRNAData_in.shape)
    #data_in = mRNAData_in.shape[0]
    
    index = mRNAData_in.index.values
    index_new = [row[:12] for row in index]
    mRNAData_in.index = index_new
    mRNAData_in = mRNAData_in.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    # print('The shape of input data after dropping duplicate is:')
    # print(mRNAData_in.shape)
    # data_in_drop = mRNAData_in.shape[0]
    
    GAData = [pd.read_excel(GADataPath,
                         disease, usecols='A,E',
                         index_col='Patient_ID',
                         keep_default_na=False) 
           for disease in tumor_types(Cancer_type)]
    #GAData_all = pd.concat(GAData)
    # print('The shape of genetic ancestry data is:')
    # print(GAData_all.shape)
    #GAdata_in = GAData_all.shape[0]
    
    # print('fetching Clinical Outcomes data...')
    if Target=='OS':
        cols = 'B,H,I'
    elif Target == 'DSS':
        cols = 'B,J,K'
    elif Target == 'DFI':
        cols = 'B,L,M'
    elif Target == 'PFI':
        cols = 'B,N,O'
    OutcomeData = pd.read_excel(OutcomeDataPath,
                           'TCGA-CDR', usecols=cols,
                           index_col='bcr_patient_barcode')
    
    # adding race information...
    GAData_race, count_WHITE, count_BLACK, count_ASIAN, count_NAT_A, count_OTHER = add_race_info_GA(GAData, groups)
    # print('The shape of GA data is:')
    # print(GAData_race.shape)
    count_race_GA = [count_WHITE, count_BLACK, count_ASIAN, count_NAT_A, count_OTHER]
    
    # adding clinical outcome endpoints data...
    OutcomeData = add_clinical_outcomes_TCGA(OutcomeData)
    # print('The shape of Outcome data is:')
    # print(OutcomeData.shape)
    
    # Keep patients with race information
    mRNAData_in = keep_patients_race_info(mRNAData_in,GAData_race,OutcomeData)
    # print('The shape of patients with race data is:')
    # print(mRNAData_in.shape)
    
    # data
    Data = mRNAData_in
    C = Data['C'].tolist()
    R = Data['race'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'race', 'T'])
    X = Data.values
    X = X.astype('float32')
    InputData = {'X': X,
                  'T': np.asarray(T, dtype=np.float32),
                  'C': np.asarray(C, dtype=np.int32),
                  'E': np.asarray(E, dtype=np.int32),
                  'R': np.asarray(R),
                  'Samples': Data.index.values,
                  'FeatureName': list(Data)}
    
    return InputData, count_race_GA

def readProteinData(Cancer_type, Target, groups, DataPath):
    
    # data paths
    ProteinDataPath = DataPath + 'ProteinData/Protein-Expression-Data.txt'
    GADataPath = DataPath + 'Genetic_Ancestry.xlsx'
    OutcomeDataPath = DataPath + 'TCGA-Clinical-Outcome-Endpoints.xlsx'
    
    df = pd.read_csv(ProteinDataPath, sep='\t', index_col='SampleID')
    # print('shape of protein data before dropping columns:')
    # print(np.shape(df))
    
    df = df.dropna(axis=1)
    # print('shape of protein data after dropping columns:')
    # print(np.shape(df))
    
    tumorTypes = tumor_types(Cancer_type)
    df = df[df['TumorType'].isin(tumorTypes)]
    # print('shape of protein data for input cancer type:')
    # print(np.shape(df))
    #data_in = np.shape(df)[0]
    
    df = df.drop(columns=['TumorType'])
    index = df.index.values
    index_new = [row[:12] for row in index]
    df.index = index_new
    ProteinData_in = df
    # ProteinData_in = df.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    # print('The shape of input data after dropping duplicate is:')
    # print(ProteinData_in.shape)
    # data_in_drop = ProteinData_in.shape[0]
    
    GAData = [pd.read_excel(GADataPath,
                         disease, usecols='A,E',
                         index_col='Patient_ID',
                         keep_default_na=False) 
           for disease in tumor_types(Cancer_type)]
    #GAData_all = pd.concat(GAData)
    # print('The shape of genetic ancestry data is:')
    # print(GAData_all.shape)
    #GAdata_in = GAData_all.shape[0]
    
    # print('fetching Clinical Outcomes data...')
    if Target=='OS':
        cols = 'B,H,I'
    elif Target == 'DSS':
        cols = 'B,J,K'
    elif Target == 'DFI':
        cols = 'B,L,M'
    elif Target == 'PFI':
        cols = 'B,N,O'
    OutcomeData = pd.read_excel(OutcomeDataPath,
                           'TCGA-CDR', usecols=cols,
                           index_col='bcr_patient_barcode')
    
    # adding race information...
    GAData_race, count_WHITE, count_BLACK, count_ASIAN, count_NAT_A, count_OTHER = add_race_info_GA(GAData, groups)
    # print('The shape of GA data is:')
    # print(GAData_race.shape)
    count_race_GA = [count_WHITE, count_BLACK, count_ASIAN, count_NAT_A, count_OTHER]
    
    # adding clinical outcome endpoints data...
    OutcomeData = add_clinical_outcomes_TCGA(OutcomeData)
    # print('The shape of Outcome data is:')
    # print(OutcomeData.shape)
    
    # Keep patients with race information
    ProteinData_in = keep_patients_race_info(ProteinData_in,GAData_race,OutcomeData)
    # print('The shape of patients with race data is:')
    # print(ProteinData_in.shape)
    
    # data
    Data = ProteinData_in
    C = Data['C'].tolist()
    R = Data['race'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'race', 'T'])
    X = Data.values
    X = X.astype('float32')
    InputData = {'X': X,
                  'T': np.asarray(T, dtype=np.float32),
                  'C': np.asarray(C, dtype=np.int32),
                  'E': np.asarray(E, dtype=np.int32),
                  'R': np.asarray(R),
                  'Samples': Data.index.values,
                  'FeatureName': list(Data)}
    
    return InputData, count_race_GA

def readMutationData(Cancer_type, Target, groups, DataPath):
    
    # data paths
    MutationDataPath = DataPath + 'MutationData/MutationData.xlsx'
    GADataPath = DataPath + 'Genetic_Ancestry.xlsx'
    OutcomeDataPath = DataPath + 'TCGA-Clinical-Outcome-Endpoints.xlsx'
    
    df = pd.read_excel(MutationDataPath,
                         index_col='Patient_ID',
                         keep_default_na=False)
    
    tumorTypes = tumor_types(Cancer_type)
    df = df[df['TumorType'].isin(tumorTypes)]
    #data_in = np.shape(df)[0]
    
    df = df.drop(columns=['TumorType'])
    # index = df.index.values
    # index_new = [row[:12] for row in index]
    # df.index = index_new
    ProteinData_in = df
    
    GAData = [pd.read_excel(GADataPath,
                         disease, usecols='A,E',
                         index_col='Patient_ID',
                         keep_default_na=False) 
           for disease in tumor_types(Cancer_type)]
    #GAData_all = pd.concat(GAData)
    # print('The shape of genetic ancestry data is:')
    # print(GAData_all.shape)
    #GAdata_in = GAData_all.shape[0]
    
    # print('fetching Clinical Outcomes data...')
    if Target=='OS':
        cols = 'B,H,I'
    elif Target == 'DSS':
        cols = 'B,J,K'
    elif Target == 'DFI':
        cols = 'B,L,M'
    elif Target == 'PFI':
        cols = 'B,N,O'
    OutcomeData = pd.read_excel(OutcomeDataPath,
                           'TCGA-CDR', usecols=cols,
                           index_col='bcr_patient_barcode')
    
    # adding race information...
    GAData_race, count_WHITE, count_BLACK, count_ASIAN, count_NAT_A, count_OTHER = add_race_info_GA(GAData, groups)
    # print('The shape of GA data is:')
    # print(GAData_race.shape)
    count_race_GA = [count_WHITE, count_BLACK, count_ASIAN, count_NAT_A, count_OTHER]
    
    # adding clinical outcome endpoints data...
    OutcomeData = add_clinical_outcomes_TCGA(OutcomeData)
    # print('The shape of Outcome data is:')
    # print(OutcomeData.shape)
    
    # Keep patients with race information
    ProteinData_in = keep_patients_race_info(ProteinData_in,GAData_race,OutcomeData)
    # print('The shape of patients with race data is:')
    # print(ProteinData_in.shape)
    
    # data
    Data = ProteinData_in
    C = Data['C'].tolist()
    R = Data['race'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'race', 'T'])
    X = Data.values
    X = X.astype('float32')
    InputData = {'X': X,
                  'T': np.asarray(T, dtype=np.float32),
                  'C': np.asarray(C, dtype=np.int32),
                  'E': np.asarray(E, dtype=np.int32),
                  'R': np.asarray(R),
                  'Samples': Data.index.values,
                  'FeatureName': list(Data)}
    
    return InputData, count_race_GA

def readProteinDataGender(Cancer_type, Target, genders, DataPath):
    
    # data paths
    ProteinDataPath = DataPath + 'ProteinData/Protein-Expression-Data.txt'
    OutcomeDataPath = DataPath + 'TCGA-Clinical-Outcome-Endpoints.xlsx'
    
    df = pd.read_csv(ProteinDataPath, sep='\t', index_col='SampleID')
    df = df.dropna(axis=1)
    tumorTypes = tumor_types(Cancer_type)
    df = df[df['TumorType'].isin(tumorTypes)]
    df = df.drop(columns=['TumorType'])
    index = df.index.values
    index_new = [row[:12] for row in index]
    df.index = index_new
    ProteinData_in = df
    
    if Target=='OS':
        cols = 'B,E,H,I'
    elif Target == 'DSS':
        cols = 'B,E,J,K'
    elif Target == 'DFI':
        cols = 'B,E,L,M'
    elif Target == 'PFI':
        cols = 'B,E,N,O'
    OutcomeData = pd.read_excel(OutcomeDataPath,
                           'TCGA-CDR', usecols=cols,
                           index_col='bcr_patient_barcode')
    OutcomeData_gender = OutcomeData[OutcomeData['gender'].isin(genders)]
    
    OutcomeData_gender.columns = ['G', 'E', 'T']
    OutcomeData_gender = OutcomeData_gender[OutcomeData_gender['E'].isin([0, 1])]
    OutcomeData_gender = OutcomeData_gender.dropna()
    OutcomeData_gender['C'] = 1 - OutcomeData_gender['E']
    OutcomeData_gender.drop(columns=['E'], inplace=True)
    
    ProteinData_in = ProteinData_in.join(OutcomeData_gender, how='inner')
    ProteinData_in = ProteinData_in.dropna(axis='columns')
    
    # data
    Data = ProteinData_in
    C = Data['C'].tolist()
    G = Data['G'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'G', 'T'])
    X = Data.values
    X = X.astype('float32')
    InputData = {'X': X,
                  'T': np.asarray(T, dtype=np.float32),
                  'C': np.asarray(C, dtype=np.int32),
                  'E': np.asarray(E, dtype=np.int32),
                  'G': np.asarray(G),
                  'Samples': Data.index.values,
                  'FeatureName': list(Data)}
    
    return InputData

def readmRNADataGender(Cancer_type, Target, genders, DataPath):
    
    # data paths
    mRNADataPath = DataPath + 'mRNAData/mRNA-Expression.mat'
    OutcomeDataPath = DataPath + 'TCGA-Clinical-Outcome-Endpoints.xlsx'
    
    mRNAData = loadmat(mRNADataPath)
    X, Y, GeneName, SampleName = mRNAData['X'].astype('float32'), mRNAData['CancerType'], mRNAData['FeatureName'][0], mRNAData['SampleName']
    GeneName = [row[0] for row in GeneName]
    SampleName = [row[0][0] for row in SampleName]
    Y = [row[0][0] for row in Y]
    
    mRNAData_X = pd.DataFrame(X, columns=GeneName, index=SampleName)
    mRNAData_Y = pd.DataFrame(Y, index=SampleName, columns=['Disease'])
    mRNAData_Y = mRNAData_Y[mRNAData_Y['Disease'].isin(tumor_types(Cancer_type))]
    mRNAData_in = mRNAData_X.join(mRNAData_Y, how='inner')
    mRNAData_in = mRNAData_in.drop(columns=['Disease'])
    
    index = mRNAData_in.index.values
    index_new = [row[:12] for row in index]
    mRNAData_in.index = index_new
    mRNAData_in = mRNAData_in.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    
    if Target=='OS':
        cols = 'B,E,H,I'
    elif Target == 'DSS':
        cols = 'B,E,J,K'
    elif Target == 'DFI':
        cols = 'B,E,L,M'
    elif Target == 'PFI':
        cols = 'B,E,N,O'
    OutcomeData = pd.read_excel(OutcomeDataPath,
                           'TCGA-CDR', usecols=cols,
                           index_col='bcr_patient_barcode')
    OutcomeData_gender = OutcomeData[OutcomeData['gender'].isin(genders)]
    
    OutcomeData_gender.columns = ['G', 'E', 'T']
    OutcomeData_gender = OutcomeData_gender[OutcomeData_gender['E'].isin([0, 1])]
    OutcomeData_gender = OutcomeData_gender.dropna()
    OutcomeData_gender['C'] = 1 - OutcomeData_gender['E']
    OutcomeData_gender.drop(columns=['E'], inplace=True)
    
    mRNAData_in = mRNAData_in.join(OutcomeData_gender, how='inner')
    mRNAData_in = mRNAData_in.dropna(axis='columns')
    
    # data
    Data = mRNAData_in
    C = Data['C'].tolist()
    G = Data['G'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'G', 'T'])
    X = Data.values
    X = X.astype('float32')
    InputData = {'X': X,
                  'T': np.asarray(T, dtype=np.float32),
                  'C': np.asarray(C, dtype=np.int32),
                  'E': np.asarray(E, dtype=np.int32),
                  'G': np.asarray(G),
                  'Samples': Data.index.values,
                  'FeatureName': list(Data)}
    
    return InputData

def readMicroRNADataGender(Cancer_type, Target, genders, DataPath):
    
    # data paths
    MicroRNADataPath = DataPath + 'MicroRNAData/MicroRNA-Expression.mat'
    OutcomeDataPath = DataPath + 'TCGA-Clinical-Outcome-Endpoints.xlsx'
    #
    MicroRNAData = loadmat(MicroRNADataPath)
    X, Y, GeneName, SampleName = MicroRNAData['X'].astype('float32'), MicroRNAData['CancerType'], MicroRNAData['FeatureName'][0], MicroRNAData['SampleName']
    GeneName = [row[0] for row in GeneName]
    SampleName = [row[0][0] for row in SampleName]
    Y = [row[0][0] for row in Y]
    
    MicroRNAData_X = pd.DataFrame(X, columns=GeneName, index=SampleName)
    MicroRNAData_Y = pd.DataFrame(Y, index=SampleName, columns=['Disease'])
    MicroRNAData_Y = MicroRNAData_Y[MicroRNAData_Y['Disease'].isin(tumor_types(Cancer_type))]
    MicroRNAData_in = MicroRNAData_X.join(MicroRNAData_Y, how='inner')
    MicroRNAData_in = MicroRNAData_in.drop(columns=['Disease'])
    
    index = MicroRNAData_in.index.values
    index_new = [row[:12] for row in index]
    MicroRNAData_in.index = index_new
    MicroRNAData_in = MicroRNAData_in.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    
    if Target=='OS':
        cols = 'B,E,H,I'
    elif Target == 'DSS':
        cols = 'B,E,J,K'
    elif Target == 'DFI':
        cols = 'B,E,L,M'
    elif Target == 'PFI':
        cols = 'B,E,N,O'
    OutcomeData = pd.read_excel(OutcomeDataPath,
                           'TCGA-CDR', usecols=cols,
                           index_col='bcr_patient_barcode')
    OutcomeData_gender = OutcomeData[OutcomeData['gender'].isin(genders)]
    
    OutcomeData_gender.columns = ['G', 'E', 'T']
    OutcomeData_gender = OutcomeData_gender[OutcomeData_gender['E'].isin([0, 1])]
    OutcomeData_gender = OutcomeData_gender.dropna()
    OutcomeData_gender['C'] = 1 - OutcomeData_gender['E']
    OutcomeData_gender.drop(columns=['E'], inplace=True)
    
    MicroRNAData_in = MicroRNAData_in.join(OutcomeData_gender, how='inner')
    MicroRNAData_in = MicroRNAData_in.dropna(axis='columns')
    
    # data
    Data = MicroRNAData_in
    C = Data['C'].tolist()
    G = Data['G'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'G', 'T'])
    X = Data.values
    X = X.astype('float32')
    InputData = {'X': X,
                  'T': np.asarray(T, dtype=np.float32),
                  'C': np.asarray(C, dtype=np.int32),
                  'E': np.asarray(E, dtype=np.int32),
                  'G': np.asarray(G),
                  'Samples': Data.index.values,
                  'FeatureName': list(Data)}
    
    return InputData

def readProteinDataGenderRace(Cancer_type, Target, genders, groups, DataPath):
    
    # data paths
    ProteinDataPath = DataPath + 'ProteinData/Protein-Expression-Data.txt'
    GADataPath = DataPath + 'Genetic_Ancestry.xlsx'
    OutcomeDataPath = DataPath + 'TCGA-Clinical-Outcome-Endpoints.xlsx'
    
    df = pd.read_csv(ProteinDataPath, sep='\t', index_col='SampleID')
    df = df.dropna(axis=1)
    tumorTypes = tumor_types(Cancer_type)
    df = df[df['TumorType'].isin(tumorTypes)]
    df = df.drop(columns=['TumorType'])
    index = df.index.values
    index_new = [row[:12] for row in index]
    df.index = index_new
    ProteinData_in = df
    
    GAData = [pd.read_excel(GADataPath,
                         disease, usecols='A,E',
                         index_col='Patient_ID',
                         keep_default_na=False) 
           for disease in tumor_types(Cancer_type)]
    GAData_race = pd.concat(GAData)
    GAData_race = GAData_race[GAData_race['EIGENSTRAT'].isin(['EA', 'AA', 'EAA', 'NA', 'OA'])]
    GAData_race['race'] = GAData_race['EIGENSTRAT']
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'EA', 'race'] = 'WHITE'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'AA', 'race'] = 'BLACK'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'EAA', 'race'] = 'ASIAN'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'NA', 'race'] = 'NAT_A'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'OA', 'race'] = 'OTHER'
    GAData_race = GAData_race.drop(columns=['EIGENSTRAT'])
    GAData_race = GAData_race[GAData_race['race'].isin(groups)]
    
    # print('fetching Clinical Outcomes data...')
    if Target=='OS':
        cols = 'B,E,H,I'
    elif Target == 'DSS':
        cols = 'B,E,J,K'
    elif Target == 'DFI':
        cols = 'B,E,L,M'
    elif Target == 'PFI':
        cols = 'B,E,N,O'
    OutcomeData = pd.read_excel(OutcomeDataPath,
                           'TCGA-CDR', usecols=cols,
                           index_col='bcr_patient_barcode')
    
    OutcomeData.columns = ['G', 'E', 'T']
    OutcomeData = OutcomeData[OutcomeData['G'].isin(genders)]
    OutcomeData = OutcomeData[OutcomeData['E'].isin([0, 1])]
    OutcomeData = OutcomeData.dropna()
    OutcomeData['C'] = 1 - OutcomeData['E']
    OutcomeData.drop(columns=['E'], inplace=True)
    
    ProteinData_in = ProteinData_in.join(GAData_race, how='inner')
    ProteinData_in = ProteinData_in.dropna(axis='columns')
    ProteinData_in = ProteinData_in.join(OutcomeData, how='inner')
    
    # data
    Data = ProteinData_in
    C = Data['C'].tolist()
    R = Data['race'].tolist()
    G = Data['G'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'race', 'T', 'G'])
    X = Data.values
    X = X.astype('float32')
    InputData = {'X': X,
                  'T': np.asarray(T, dtype=np.float32),
                  'C': np.asarray(C, dtype=np.int32),
                  'E': np.asarray(E, dtype=np.int32),
                  'R': np.asarray(R),
                  'G': np.asarray(G),
                  'Samples': Data.index.values,
                  'FeatureName': list(Data)}
    
    return InputData    

def readmRNADataGenderRace(Cancer_type, Target, genders, groups, DataPath):
    
    # data paths
    mRNADataPath = DataPath + 'mRNAData/mRNA-Expression.mat'
    GADataPath = DataPath + 'Genetic_Ancestry.xlsx'
    OutcomeDataPath = DataPath + 'TCGA-Clinical-Outcome-Endpoints.xlsx'
    #
    mRNAData = loadmat(mRNADataPath)
    X, Y, GeneName, SampleName = mRNAData['X'].astype('float32'), mRNAData['CancerType'], mRNAData['FeatureName'][0], mRNAData['SampleName']
    GeneName = [row[0] for row in GeneName]
    SampleName = [row[0][0] for row in SampleName]
    Y = [row[0][0] for row in Y]
    
    mRNAData_X = pd.DataFrame(X, columns=GeneName, index=SampleName)
    mRNAData_Y = pd.DataFrame(Y, index=SampleName, columns=['Disease'])
    mRNAData_Y = mRNAData_Y[mRNAData_Y['Disease'].isin(tumor_types(Cancer_type))]
    mRNAData_in = mRNAData_X.join(mRNAData_Y, how='inner')
    mRNAData_in = mRNAData_in.drop(columns=['Disease'])
    
    index = mRNAData_in.index.values
    index_new = [row[:12] for row in index]
    mRNAData_in.index = index_new
    mRNAData_in = mRNAData_in.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    
    GAData = [pd.read_excel(GADataPath,
                         disease, usecols='A,E',
                         index_col='Patient_ID',
                         keep_default_na=False) 
           for disease in tumor_types(Cancer_type)]
    GAData_race = pd.concat(GAData)
    GAData_race = GAData_race[GAData_race['EIGENSTRAT'].isin(['EA', 'AA', 'EAA', 'NA', 'OA'])]
    GAData_race['race'] = GAData_race['EIGENSTRAT']
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'EA', 'race'] = 'WHITE'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'AA', 'race'] = 'BLACK'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'EAA', 'race'] = 'ASIAN'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'NA', 'race'] = 'NAT_A'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'OA', 'race'] = 'OTHER'
    GAData_race = GAData_race.drop(columns=['EIGENSTRAT'])
    GAData_race = GAData_race[GAData_race['race'].isin(groups)]
    
    # print('fetching Clinical Outcomes data...')
    if Target=='OS':
        cols = 'B,E,H,I'
    elif Target == 'DSS':
        cols = 'B,E,J,K'
    elif Target == 'DFI':
        cols = 'B,E,L,M'
    elif Target == 'PFI':
        cols = 'B,E,N,O'
    OutcomeData = pd.read_excel(OutcomeDataPath,
                           'TCGA-CDR', usecols=cols,
                           index_col='bcr_patient_barcode')
    
    OutcomeData.columns = ['G', 'E', 'T']
    OutcomeData = OutcomeData[OutcomeData['G'].isin(genders)]
    OutcomeData = OutcomeData[OutcomeData['E'].isin([0, 1])]
    OutcomeData = OutcomeData.dropna()
    OutcomeData['C'] = 1 - OutcomeData['E']
    OutcomeData.drop(columns=['E'], inplace=True)
    
    mRNAData_in = mRNAData_in.join(GAData_race, how='inner')
    mRNAData_in = mRNAData_in.dropna(axis='columns')
    mRNAData_in = mRNAData_in.join(OutcomeData, how='inner')
    
    # data
    Data = mRNAData_in
    C = Data['C'].tolist()
    R = Data['race'].tolist()
    G = Data['G'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'race', 'T', 'G'])
    X = Data.values
    X = X.astype('float32')
    InputData = {'X': X,
                  'T': np.asarray(T, dtype=np.float32),
                  'C': np.asarray(C, dtype=np.int32),
                  'E': np.asarray(E, dtype=np.int32),
                  'R': np.asarray(R),
                  'G': np.asarray(G),
                  'Samples': Data.index.values,
                  'FeatureName': list(Data)}
    
    return InputData    

def readMicroRNADataGenderRace(Cancer_type, Target, genders, groups, DataPath):
    
    # data paths
    MicroRNADataPath = DataPath + 'MicroRNAData/MicroRNA-Expression.mat'
    GADataPath = DataPath + 'Genetic_Ancestry.xlsx'
    OutcomeDataPath = DataPath + 'TCGA-Clinical-Outcome-Endpoints.xlsx'
    #
    MicroRNAData = loadmat(MicroRNADataPath)
    X, Y, GeneName, SampleName = MicroRNAData['X'].astype('float32'), MicroRNAData['CancerType'], MicroRNAData['FeatureName'][0], MicroRNAData['SampleName']
    GeneName = [row[0] for row in GeneName]
    SampleName = [row[0][0] for row in SampleName]
    Y = [row[0][0] for row in Y]
    
    MicroRNAData_X = pd.DataFrame(X, columns=GeneName, index=SampleName)
    MicroRNAData_Y = pd.DataFrame(Y, index=SampleName, columns=['Disease'])
    MicroRNAData_Y = MicroRNAData_Y[MicroRNAData_Y['Disease'].isin(tumor_types(Cancer_type))]
    MicroRNAData_in = MicroRNAData_X.join(MicroRNAData_Y, how='inner')
    MicroRNAData_in = MicroRNAData_in.drop(columns=['Disease'])
    
    index = MicroRNAData_in.index.values
    index_new = [row[:12] for row in index]
    MicroRNAData_in.index = index_new
    MicroRNAData_in = MicroRNAData_in.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    
    GAData = [pd.read_excel(GADataPath,
                         disease, usecols='A,E',
                         index_col='Patient_ID',
                         keep_default_na=False) 
           for disease in tumor_types(Cancer_type)]
    GAData_race = pd.concat(GAData)
    GAData_race = GAData_race[GAData_race['EIGENSTRAT'].isin(['EA', 'AA', 'EAA', 'NA', 'OA'])]
    GAData_race['race'] = GAData_race['EIGENSTRAT']
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'EA', 'race'] = 'WHITE'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'AA', 'race'] = 'BLACK'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'EAA', 'race'] = 'ASIAN'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'NA', 'race'] = 'NAT_A'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'OA', 'race'] = 'OTHER'
    GAData_race = GAData_race.drop(columns=['EIGENSTRAT'])
    GAData_race = GAData_race[GAData_race['race'].isin(groups)]
    
    # print('fetching Clinical Outcomes data...')
    if Target=='OS':
        cols = 'B,E,H,I'
    elif Target == 'DSS':
        cols = 'B,E,J,K'
    elif Target == 'DFI':
        cols = 'B,E,L,M'
    elif Target == 'PFI':
        cols = 'B,E,N,O'
    OutcomeData = pd.read_excel(OutcomeDataPath,
                           'TCGA-CDR', usecols=cols,
                           index_col='bcr_patient_barcode')
    
    OutcomeData.columns = ['G', 'E', 'T']
    OutcomeData = OutcomeData[OutcomeData['G'].isin(genders)]
    OutcomeData = OutcomeData[OutcomeData['E'].isin([0, 1])]
    OutcomeData = OutcomeData.dropna()
    OutcomeData['C'] = 1 - OutcomeData['E']
    OutcomeData.drop(columns=['E'], inplace=True)
    
    MicroRNAData_in = MicroRNAData_in.join(GAData_race, how='inner')
    MicroRNAData_in = MicroRNAData_in.dropna(axis='columns')
    MicroRNAData_in = MicroRNAData_in.join(OutcomeData, how='inner')
    
    # data
    Data = MicroRNAData_in
    C = Data['C'].tolist()
    R = Data['race'].tolist()
    G = Data['G'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'race', 'T', 'G'])
    X = Data.values
    X = X.astype('float32')
    InputData = {'X': X,
                  'T': np.asarray(T, dtype=np.float32),
                  'C': np.asarray(C, dtype=np.int32),
                  'E': np.asarray(E, dtype=np.int32),
                  'R': np.asarray(R),
                  'G': np.asarray(G),
                  'Samples': Data.index.values,
                  'FeatureName': list(Data)}
    
    return InputData    

def run_cv_ancestry(Cancer_type, Feature_type, Target, groups, DataPath):
    
    if Feature_type=='mRNA':
        
        Data, count_race_GA = readmRNAData(Cancer_type, Target, groups, DataPath)
    
    elif Feature_type=='Protein':
        
        Data, count_race_GA = readProteinData(Cancer_type, Target, groups, DataPath)
        
    elif Feature_type=='Methylation':
        
        Data, count_race_GA = readMethylationData(Cancer_type, Target, groups, DataPath)
        
    elif Feature_type=='MicroRNA':
        
        Data, count_race_GA = readMicroRNAData(Cancer_type, Target, groups, DataPath)
        
    elif Feature_type=='Mutation':
        
        Data, count_race_GA = readMutationData(Cancer_type, Target, groups, DataPath)
        
    return Data, count_race_GA

def run_cv_gender(Cancer_type, Feature_type, Target, genders, DataPath):
    
    if Feature_type=='Protein':
        
        Data = readProteinDataGender(Cancer_type, Target, genders, DataPath)
        
    elif Feature_type=='mRNA':
        
        Data = readmRNADataGender(Cancer_type, Target, genders, DataPath)
    
    elif Feature_type=='MicroRNA':
        
        Data = readMicroRNADataGender(Cancer_type, Target, genders, DataPath)
    
    return Data

def run_cv_gender_race(Cancer_type, Feature_type, Target, genders, groups, DataPath):
    
    if Feature_type=='Protein':
        
        Data = readProteinDataGenderRace(Cancer_type, Target, genders, groups, DataPath)
    
    if Feature_type=='mRNA':
        
        Data = readmRNADataGenderRace(Cancer_type, Target, genders, groups, DataPath)
    
    if Feature_type=='MicroRNA':
        
        Data = readMicroRNADataGenderRace(Cancer_type, Target, genders, groups, DataPath)
    
    if Feature_type=='Methylation':
        
        Data = readMethylationDataGenderRace(Cancer_type, Target, genders, groups, DataPath)
    
    return Data

def merge_datasets(datasets):

    data = datasets[0]
    X, T, C, E, R, G, Samples, FeatureName = data['X'], data['T'], data['C'], data['E'], data['R'], data['G'], data['Samples'], data['FeatureName']
    df = pd.DataFrame(X, index=Samples, columns=FeatureName)
    df['T'] = T
    df['C'] = C
    df['E'] = E
    df['R'] = R
    df['G'] = G

    for i in range(1, len(datasets)):
        data1 = datasets[i]
        X1, Samples, FeatureName = data1['X'], data1['Samples'], data1['FeatureName']
        temp = pd.DataFrame(X1, index=Samples, columns=FeatureName)
        df = df.join(temp, how='inner')

    # Packing the data and save it to the disk
    C = df['C'].tolist()
    R = df['R'].tolist()
    G = df['G'].tolist()
    T = df['T'].tolist()
    E = df['E'].tolist()
    df = df.drop(columns=['C', 'R', 'G', 'T', 'E'])
    X = df.values
    X = X.astype('float32')
    data = {'X': X,
            'T': np.asarray(T, dtype=np.float32),
            'C': np.asarray(C, dtype=np.int32),
            'E': np.asarray(E, dtype=np.int32),
            'R': np.asarray(R),
            'G': np.asarray(G),
            'Samples': df.index.values,
            'FeatureName': list(df)}

    return data

def run_cv_gender_race_comb(Cancer_type, Feature_type, Target, genders, groups, DataPath):
    
    datasets = []
    for feature in Feature_type:
        if feature=='Protein':
            Data = readProteinDataGenderRace(Cancer_type, Target, genders, groups, DataPath)
        if feature=='mRNA':
            Data = readmRNADataGenderRace(Cancer_type, Target, genders, groups, DataPath)
        if feature=='MicroRNA':
            Data = readMicroRNADataGenderRace(Cancer_type, Target, genders, groups, DataPath)
        if feature=='Methylation':
            Data = readMethylationDataGenderRace(Cancer_type, Target, genders, groups, DataPath)
        datasets.append(Data)
        
    return merge_datasets(datasets)


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("which_count", type=str, help="gender_race, ancestry, gender")
    parser.add_argument("count_type", type=str, help="single, combination")
    parser.add_argument("num_features", type=int, help="how many features to be combined (2,3,4)")
    parser.add_argument("min_cases", type=int, help="Normally, it is 5")
    parser.add_argument("DDP_group", type=str, help="BLACK, ASIAN, NAT_A")
    args = parser.parse_args()
    print(args)
    
    which_count = args.which_count
    count_type = args.count_type
    num_features = args.num_features
    min_cases = args.min_cases
    DDP_group = args.DDP_group
    
    groups = ("WHITE",DDP_group)
    genders = ("MALE","FEMALE")
    
    pos_col_s_ancest = '1'+groups[0]
    neg_col_s_ancest = '0'+groups[0]
    pos_col_t_ancest = '1'+groups[1]
    neg_col_t_ancest = '0'+groups[1]
    pos_col_s_gender = '1'+genders[0]
    neg_col_s_gender = '0'+genders[0]
    pos_col_t_gender = '1'+genders[1]
    neg_col_t_gender = '0'+genders[1]
    
    Features = ['Protein','mRNA','MicroRNA','Methylation']
    Cancers = ['ACC','BLCA','BRCA','CESC',
            'CHOL','COAD','DLBC','ESCA',
            'GBM','HNSC','KICH','KIRC',
            'KIRP','LAML','LGG','LIHC',
            'LUAD','LUSC','MESO','OV',
            'PAAD','PCPG','PRAD','READ',
            'SARC','SKCM','STAD','TGCT',
            'THCA','THYM','UCEC','UCS','UVM',
            'GBMLGG','COADREAD','KIPAN','STES',
            'PanGI','PanGyn','PanSCCs']
    years_val = [1,2,3,4,5]
    Targets = ['OS','DSS','PFI','DFI']
    
    if count_type=='single':
        count_out = pd.DataFrame()
        count_out_data = pd.DataFrame()
        
        for i_f in range(np.shape(Features)[0]):
            Feature_type = Features[i_f]
            
            for i_t in range(np.shape(Targets)[0]):
                Target = Targets[i_t]
                
                for i_c in range(np.shape(Cancers)[0]):
                    Cancer_type = Cancers[i_c]
                    
                    for i_y in range(np.shape(years_val)[0]):
                        Years = years_val[i_y]
                        
                        print('=================================')
                        print(Cancers[i_c],Features[i_f],Targets[i_t],years_val[i_y])
                        print('=================================')
                        
                        if which_count=='ancestry':
                            dataset, count_race_GA = run_cv_ancestry(Cancer_type, Feature_type, Target, groups, DataPath)
                        elif which_count=='gender':
                            dataset = run_cv_gender(Cancer_type, Feature_type, Target, genders, DataPath)
                        elif which_count=='gender_race':
                            dataset = run_cv_gender_race(Cancer_type, Feature_type, Target, genders, groups, DataPath)
                        
                        if (np.shape(dataset['X'])[0])!=0:
                            
                            if which_count=='ancestry':
                                count_race_final = race_count(dataset['R'])
                                dataset = get_n_years(dataset, Years)
                                X, Y, R, y_sub, y_strat = dataset
                                df = pd.DataFrame(y_strat, columns=['RY'])
                                df['R'] = R
                                df['Y'] = Y
                                print(X.shape)
                                print(df['RY'].value_counts())#race with prognosis counts
                                print(df['R'].value_counts())#race counts
                                print(df['Y'].value_counts())#progonsis counts
                                
                                prog_count = get_prognosis_count(dataset,groups)
                                
                                if (prog_count[0]>=min_cases and 
                                    prog_count[1]>=min_cases and 
                                    prog_count[2]>=min_cases and 
                                    prog_count[3]>=min_cases):
                                    
                                    count_out = pd.DataFrame({"Cancer_type":[Cancer_type],
                                                          "Feature_type":[Feature_type],
                                                          "Target":[Target],
                                                          "Years":Years,
                                                          groups[0]:count_race_final[groups[0]][0],
                                                          groups[1]:count_race_final[groups[1]][0],
                                                          'Y(0)':prog_count[1]+prog_count[3],
                                                          'Y(1)':prog_count[0]+prog_count[2],
                                                          pos_col_s_ancest:prog_count[0],
                                                          neg_col_s_ancest:prog_count[1],
                                                          pos_col_t_ancest:prog_count[2],
                                                          neg_col_t_ancest:prog_count[3]})
                                    
                                    count_out_data = count_out_data.append(count_out,ignore_index=True)
                                    
                            elif which_count=='gender':
                                dataset = get_n_years_gender(dataset, Years)
                                X, Y, G, y_sub, y_strat = dataset
                                df = pd.DataFrame(y_strat, columns=['GY'])
                                df['G'] = G
                                df['Y'] = Y
                                print(X.shape)
                                print(df['GY'].value_counts())#race with prognosis counts
                                print(df['G'].value_counts())#race counts
                                print(df['Y'].value_counts())#progonsis counts
                                
                                prog_count = get_prognosis_count(dataset,genders)
                                
                                if (prog_count[0]>=min_cases and 
                                    prog_count[1]>=min_cases and 
                                    prog_count[2]>=min_cases and 
                                    prog_count[3]>=min_cases):
                                    
                                    count_out = pd.DataFrame({"Cancer_type":[Cancer_type],
                                                          "Feature_type":[Feature_type],
                                                          "Target":[Target],
                                                          "Years":Years,
                                                          genders[0]:prog_count[0]+prog_count[1],
                                                          genders[1]:prog_count[2]+prog_count[3],
                                                          'Y(0)':prog_count[1]+prog_count[3],
                                                          'Y(1)':prog_count[0]+prog_count[2],
                                                          pos_col_s_gender:prog_count[0],
                                                          neg_col_s_gender:prog_count[1],
                                                          pos_col_t_gender:prog_count[2],
                                                          neg_col_t_gender:prog_count[3]})
                                    
                                    count_out_data = count_out_data.append(count_out,ignore_index=True)
                            
                            elif which_count=='gender_race':
                                dataset = get_n_years_gender_race(dataset, Years)
                                X, Y, R, y_sub, y_strat, G, Gy_strat, GRy_strat = dataset
                                df = pd.DataFrame(y_strat, columns=['RY'])
                                df['GRY'] = GRy_strat
                                df['GY'] = Gy_strat
                                df['R'] = R
                                df['G'] = G
                                df['Y'] = Y
                                print(X.shape)
                                print(df['GRY'].value_counts())#gender with prognosis counts
                                print(df['GY'].value_counts())#gender with prognosis counts
                                print(df['G'].value_counts())#gender counts
                                print(df['RY'].value_counts())#race with prognosis counts
                                print(df['R'].value_counts())#race counts
                                print(df['Y'].value_counts())#progonsis counts
                                
                                prog_count = get_prognosis_count_RG(dataset,groups,genders)
                                
                                pos_source_male = genders[0]+'1'+groups[0]
                                neg_source_male = genders[0]+'0'+groups[0]
                                pos_target_male = genders[0]+'1'+groups[1]
                                neg_target_male = genders[0]+'0'+groups[1]
                                pos_source_female = genders[1]+'1'+groups[0]
                                neg_source_female = genders[1]+'0'+groups[0]
                                pos_target_female = genders[1]+'1'+groups[1]
                                neg_target_female = genders[1]+'0'+groups[1]
                                
                                count_out = pd.DataFrame({"Cancer_type":[Cancer_type],
                                                          "Feature_type":[Feature_type],
                                                          "Target":[Target],
                                                          "Years":Years,
                                                          genders[0]:prog_count[0]+prog_count[1]+prog_count[2]+prog_count[3],
                                                          genders[1]:prog_count[4]+prog_count[5]+prog_count[6]+prog_count[7],
                                                          groups[0]:prog_count[0]+prog_count[1]+prog_count[4]+prog_count[5],
                                                          groups[1]:prog_count[2]+prog_count[3]+prog_count[6]+prog_count[7],
                                                          'Y(0)':prog_count[1]+prog_count[3]+prog_count[5]+prog_count[7],
                                                          'Y(1)':prog_count[0]+prog_count[2]+prog_count[4]+prog_count[6],
                                                          pos_col_s_ancest:prog_count[0]+prog_count[4],
                                                          neg_col_s_ancest:prog_count[1]+prog_count[5],
                                                          pos_col_t_ancest:prog_count[2]+prog_count[6],
                                                          neg_col_t_ancest:prog_count[3]+prog_count[7],
                                                          pos_col_s_gender:prog_count[0]+prog_count[2],
                                                          neg_col_s_gender:prog_count[1]+prog_count[3],
                                                          pos_col_t_gender:prog_count[4]+prog_count[6],
                                                          neg_col_t_gender:prog_count[5]+prog_count[7],
                                                          pos_source_male:prog_count[0],
                                                          neg_source_male:prog_count[1],
                                                          pos_target_male:prog_count[2],
                                                          neg_target_male:prog_count[3],
                                                          pos_source_female:prog_count[4],
                                                          neg_source_female:prog_count[5],
                                                          pos_target_female:prog_count[6],
                                                          neg_target_female:prog_count[7]})
                                
                                count_out_data = count_out_data.append(count_out,ignore_index=True)
                                
            filename = folderISAAC+'TCGA-gender_race_count_'+Feature_type+'.xlsx'
            count_out_data.to_excel(filename)
    
    elif count_type=='combination':
        count_out = pd.DataFrame()
        count_out_data = pd.DataFrame()
        comb_feature = combinations(Features,num_features)
        # combination of num_features...
        comb_features = pd.DataFrame()
        Feature_type = pd.DataFrame()
        for cb_2 in list(comb_feature):
            # for f_i in range(0,num_features):
            #     print(cb_2[f_i])
            Feature_type = cb_2
            print('=================================')
            print(Feature_type)
            print('=================================')
            
            for i_t in range(np.shape(Targets)[0]):
                Target = Targets[i_t]
                
                for i_c in range(np.shape(Cancers)[0]):
                    Cancer_type = Cancers[i_c]
                    
                    for i_y in range(np.shape(years_val)[0]):
                        Years = years_val[i_y]
                        
                        print('=================================')
                        print(Cancers[i_c],Targets[i_t],years_val[i_y])
                        print('=================================')
                        
                        dataset = run_cv_gender_race_comb(Cancer_type, Feature_type, Target, genders, groups, DataPath)
                        
                        if (np.shape(dataset['X'])[0])!=0:
                            
                            dataset = get_n_years_gender_race(dataset, Years)
                            X, Y, R, y_sub, y_strat, G, Gy_strat, GRy_strat = dataset
                            df = pd.DataFrame(y_strat, columns=['RY'])
                            df['GRY'] = GRy_strat
                            df['GY'] = Gy_strat
                            df['R'] = R
                            df['G'] = G
                            df['Y'] = Y
                            print(X.shape)
                            print(df['GRY'].value_counts())#gender with prognosis counts
                            print(df['GY'].value_counts())#gender with prognosis counts
                            print(df['G'].value_counts())#gender counts
                            print(df['RY'].value_counts())#race with prognosis counts
                            print(df['R'].value_counts())#race counts
                            print(df['Y'].value_counts())#progonsis counts
                            
                            prog_count = get_prognosis_count_RG(dataset,groups,genders)
                            
                            pos_source_male = genders[0]+'1'+groups[0]
                            neg_source_male = genders[0]+'0'+groups[0]
                            pos_target_male = genders[0]+'1'+groups[1]
                            neg_target_male = genders[0]+'0'+groups[1]
                            pos_source_female = genders[1]+'1'+groups[0]
                            neg_source_female = genders[1]+'0'+groups[0]
                            pos_target_female = genders[1]+'1'+groups[1]
                            neg_target_female = genders[1]+'0'+groups[1]
                            
                            count_out = pd.DataFrame({"Cancer_type":[Cancer_type],
                                                      "Feature_type":[Feature_type],
                                                      "Target":[Target],
                                                      "Years":Years,
                                                      genders[0]:prog_count[0]+prog_count[1]+prog_count[2]+prog_count[3],
                                                      genders[1]:prog_count[4]+prog_count[5]+prog_count[6]+prog_count[7],
                                                      groups[0]:prog_count[0]+prog_count[1]+prog_count[4]+prog_count[5],
                                                      groups[1]:prog_count[2]+prog_count[3]+prog_count[6]+prog_count[7],
                                                      'Y(0)':prog_count[1]+prog_count[3]+prog_count[5]+prog_count[7],
                                                      'Y(1)':prog_count[0]+prog_count[2]+prog_count[4]+prog_count[6],
                                                      pos_col_s_ancest:prog_count[0]+prog_count[4],
                                                      neg_col_s_ancest:prog_count[1]+prog_count[5],
                                                      pos_col_t_ancest:prog_count[2]+prog_count[6],
                                                      neg_col_t_ancest:prog_count[3]+prog_count[7],
                                                      pos_col_s_gender:prog_count[0]+prog_count[2],
                                                      neg_col_s_gender:prog_count[1]+prog_count[3],
                                                      pos_col_t_gender:prog_count[4]+prog_count[6],
                                                      neg_col_t_gender:prog_count[5]+prog_count[7],
                                                      pos_source_male:prog_count[0],
                                                      neg_source_male:prog_count[1],
                                                      pos_target_male:prog_count[2],
                                                      neg_target_male:prog_count[3],
                                                      pos_source_female:prog_count[4],
                                                      neg_source_female:prog_count[5],
                                                      pos_target_female:prog_count[6],
                                                      neg_target_female:prog_count[7]})
                            
                            count_out_data = count_out_data.append(count_out,ignore_index=True)
                        
        filename = folderISAAC+'TCGA-gender_race_count_'+str(num_features)+\
                   '_features_combinations_'+groups[1]+\
                   '_'+Target+'.xlsx'
        count_out_data.to_excel(filename)

if __name__ == '__main__':
    main()



    
    
    
    
    
    