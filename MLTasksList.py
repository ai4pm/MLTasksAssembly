# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 14:56:01 2024

@author: teesh
"""

from scipy.io import loadmat
import pandas as pd
import numpy as np
import os
from itertools import combinations
import argparse
import sys
sys.stdout.flush()

# Define set of input conditions
Features = ['Protein', 'mRNA', 'MicroRNA', 'Methylation']
Cancers = ['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 
            'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 
            'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 
            'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ', 
            'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 
            'UCEC', 'UCS', 'UVM', 'GBMLGG', 'COADREAD', 'KIPAN', 
            'STES', 'PanGI', 'PanGyn', 'PanSCCs']
Endpoints = ['OS', 'DSS', 'PFI', 'DFI']
YearsAll = [1, 2, 3, 4, 5]

# Define global paths
folderISAAC = "./"
DataPath = 'C:/Users/teesh/OneDrive - Indian Institute of Technology Guwahati/Dataset/EssentialData/'
if not os.path.exists(DataPath):
    print("Not working in local system, let's define cluster paths...")
    DataPath = 'Dataset/EssentialData/'
    folderISAAC = "./MLTasksList/"

MethylationDataPath = os.path.join(DataPath, 'MethylationData/Methylation.mat')
MethyAncsDataPath = os.path.join(DataPath, 'MethylationData/MethylationGenetic.xlsx')
MethyCIDataPath = os.path.join(DataPath, 'MethylationData/MethylationClinInfo.xlsx')
MicroRNADataPath = os.path.join(DataPath, 'MicroRNAData/MicroRNA-Expression.mat')
mRNADataPath = os.path.join(DataPath, 'mRNAData/mRNA-Expression.mat')
ProteinDataPath = os.path.join(DataPath, 'ProteinData/Protein-Expression-Data.txt')
GADataPath = os.path.join(DataPath, 'Genetic_Ancestry.xlsx')
OutcomeDataPath = os.path.join(DataPath, 'TCGA-Clinical-Outcome-Endpoints.xlsx')

def tumor_types(Cancer_type):
    Map = {
        'GBMLGG': ['GBM', 'LGG'],
        'COADREAD': ['COAD', 'READ'],
        'KIPAN': ['KIRC', 'KICH', 'KIRP'],
        'STES': ['ESCA', 'STAD'],
        'PanGI': ['COAD', 'STAD', 'READ', 'ESCA'],
        'PanGyn': ['OV', 'CESC', 'UCS', 'UCEC'],
        'PanSCCs': ['LUSC', 'HNSC', 'ESCA', 'CESC', 'BLCA'],
        'PanPan': [
            'ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC',
            'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG',
            'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD', 'PCPG', 'PRAD', 'READ',
            'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS', 'UVM'
        ]
    }
    return Map.get(Cancer_type, [Cancer_type])

def get_prognosis_count(dataset, criteria, groups, genders):
    
    prog_count = []
    
    if criteria == "GenderRace":
        prog_data = pd.DataFrame(dataset[6]) # GRy_strat column
        for gender in genders:
            for outcome in ['1', '0']:
                for group in groups:
                    temp = prog_data.apply(lambda x: x[0] == f"{gender}{outcome}{group}", axis=1)
                    prog_count.append(len(temp[temp]))
                    
    if criteria == "Race":
        prog_data = pd.DataFrame(dataset[3]) # y_strat column
        for outcome in ['1', '0']:
            for group in groups:
                temp = prog_data.apply(lambda x: x[0] == f"{outcome}{group}", axis=1)
                prog_count.append(len(temp[temp]))
                
    if criteria == "Gender":
        prog_data = pd.DataFrame(dataset[5]) # Gy_strat column
        for gender in genders:
            for outcome in ['1', '0']:
                temp = prog_data.apply(lambda x: x[0] == f"{outcome}{gender}", axis=1)
                prog_count.append(len(temp[temp]))
                
    return prog_count

def get_n_years(dataset, years):
    
    X, T, C, E = dataset['X'], dataset['T'], dataset['C'], dataset['E']
    df = pd.DataFrame(X)
    df['T'] = T
    df['C'] = C
    df['Y'] = 1
    R, G = dataset['R'], dataset['G']
    df['R'] = R
    df['G'] = G
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
    
    return (X, Y.astype('int32'), R, y_strat, G, Gy_strat, GRy_strat)

def readMethylationData(Cancer_type, Target, genders, groups, DataPath):
    
    MethylationData = loadmat(MethylationDataPath)
    X, Y, GeneName, SampleName = MethylationData['X'].astype('float32'), MethylationData['CancerType'], MethylationData['FeatureName'][0], MethylationData['SampleName']
    GeneName = [row[0] for row in GeneName]
    SampleName = [row[0][0] for row in SampleName]
    Y = [row[0][0] for row in Y]
    
    MethylationData_X = pd.DataFrame(X, columns=GeneName, index=SampleName)
    MethylationData_Y = pd.DataFrame(Y, index=SampleName, columns=['Disease'])
    MethylationData_Y = MethylationData_Y[MethylationData_Y['Disease'].isin(tumor_types(Cancer_type))]
    MethylationData_in = MethylationData_X.join(MethylationData_Y, how='inner').drop(columns=['Disease'])
    
    MethylationData_in.index = [row[:12] for row in MethylationData_in.index.values]
    MethylationData_in = MethylationData_in.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    
    # if criteria in ["GenderRace", "Race"]:
    MethyAncsData = [
        pd.read_excel(MethyAncsDataPath, disease, usecols='A,B', index_col='bcr_patient_barcode', keep_default_na=False)
        for disease in tumor_types(Cancer_type)
    ]
    GAData_race = pd.concat(MethyAncsData)
    GAData_race = GAData_race.rename(columns={'race': 'EIGENSTRAT'})
    GAData_race['race'] = GAData_race['EIGENSTRAT']
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'WHITE', 'race'] = 'WHITE'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'BLACK OR AFRICAN AMERICAN', 'race'] = 'BLACK'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'ASIAN', 'race'] = 'ASIAN'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'AMERICAN INDIAN OR ALASKA NATIVE', 'race'] = 'NAT_A'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER', 'race'] = 'OTHER'
    GAData_race = GAData_race.drop(columns=['EIGENSTRAT'])
    if groups:
        GAData_race = GAData_race[GAData_race['race'].isin(groups)]

    if Target == 'OS':
        cols = 'A,D,Y,Z'
    elif Target == 'DSS':
        cols = 'A,D,AA,AB'
    elif Target == 'DFI':
        cols = 'A,D,AC,AD'
    elif Target == 'PFI':
        cols = 'A,D,AE,AF'
    OutcomeData = pd.read_excel(MethyCIDataPath, usecols=cols, dtype={'OS': np.float64}, index_col='bcr_patient_barcode')
    
    OutcomeData = OutcomeData[OutcomeData['gender'].isin(genders)]
    OutcomeData.columns = ['G', 'E', 'T']
    OutcomeData = OutcomeData[OutcomeData['E'].isin([0, 1])].dropna()
    OutcomeData['C'] = 1 - OutcomeData['E']
    OutcomeData.drop(columns=['E'], inplace=True)
    MethylationData_in = MethylationData_in.join(GAData_race, how='inner').dropna(axis='columns')
    MethylationData_in = MethylationData_in.join(OutcomeData, how='inner')
    Data = MethylationData_in
    C = Data['C'].tolist()
    R = Data['race'].tolist()
    G = Data['G'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'race', 'T', 'G'])
    X = Data.values.astype('float32')
    
    InputData = {'X': X, 
                 'T': np.asarray(T, dtype=np.float32), 
                 'C': np.asarray(C, dtype=np.int32), 
                 'E': np.asarray(E, dtype=np.int32), 
                 'R': np.asarray(R), 
                 'G': np.asarray(G), 
                 'Samples': Data.index.values, 
                 'FeatureName': list(Data)}
    
    return InputData

def readMicroRNAData(Cancer_type, Target, genders, groups, DataPath):
    
    MicroRNAData = loadmat(MicroRNADataPath)
    X, Y, GeneName, SampleName = MicroRNAData['X'].astype('float32'), MicroRNAData['CancerType'], MicroRNAData['FeatureName'][0], MicroRNAData['SampleName']
    GeneName = [row[0] for row in GeneName]
    SampleName = [row[0][0] for row in SampleName]
    Y = [row[0][0] for row in Y]
    
    MicroRNAData_X = pd.DataFrame(X, columns=GeneName, index=SampleName)
    MicroRNAData_Y = pd.DataFrame(Y, index=SampleName, columns=['Disease'])
    MicroRNAData_Y = MicroRNAData_Y[MicroRNAData_Y['Disease'].isin(tumor_types(Cancer_type))]
    MicroRNAData_in = MicroRNAData_X.join(MicroRNAData_Y, how='inner').drop(columns=['Disease'])
    
    MicroRNAData_in.index = [row[:12] for row in MicroRNAData_in.index.values]
    MicroRNAData_in = MicroRNAData_in.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    
    # if criteria in ["GenderRace", "Race"]:
    GAData = [
        pd.read_excel(GADataPath, disease, usecols='A,E', index_col='Patient_ID', keep_default_na=False)
        for disease in tumor_types(Cancer_type)
    ]
    GAData_race = pd.concat(GAData)
    GAData_race = GAData_race[GAData_race['EIGENSTRAT'].isin(['EA', 'AA', 'EAA', 'NA', 'OA'])]
    GAData_race['race'] = GAData_race['EIGENSTRAT']
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'EA', 'race'] = 'WHITE'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'AA', 'race'] = 'BLACK'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'EAA', 'race'] = 'ASIAN'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'NA', 'race'] = 'NAT_A'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'OA', 'race'] = 'OTHER'
    GAData_race = GAData_race.drop(columns=['EIGENSTRAT'])
    if groups:
        GAData_race = GAData_race[GAData_race['race'].isin(groups)]

    if Target == 'OS':
        cols = 'B,E,H,I'
    elif Target == 'DSS':
        cols = 'B,E,J,K'
    elif Target == 'DFI':
        cols = 'B,E,L,M'
    elif Target == 'PFI':
        cols = 'B,E,N,O'
    OutcomeData = pd.read_excel(OutcomeDataPath, 'TCGA-CDR', usecols=cols, index_col='bcr_patient_barcode')
    
    OutcomeData = OutcomeData[OutcomeData['gender'].isin(genders)]
    OutcomeData.columns = ['G', 'E', 'T']
    OutcomeData = OutcomeData[OutcomeData['E'].isin([0, 1])].dropna()
    OutcomeData['C'] = 1 - OutcomeData['E']
    OutcomeData.drop(columns=['E'], inplace=True)
    MicroRNAData_in = MicroRNAData_in.join(GAData_race, how='inner').dropna(axis='columns')
    MicroRNAData_in = MicroRNAData_in.join(OutcomeData, how='inner')
    Data = MicroRNAData_in
    C = Data['C'].tolist()
    R = Data['race'].tolist()
    G = Data['G'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'race', 'T', 'G'])
    X = Data.values.astype('float32')
    
    InputData = {'X': X, 
                 'T': np.asarray(T, dtype=np.float32), 
                 'C': np.asarray(C, dtype=np.int32), 
                 'E': np.asarray(E, dtype=np.int32), 
                 'R': np.asarray(R), 
                 'G': np.asarray(G), 
                 'Samples': Data.index.values, 
                 'FeatureName': list(Data)}
    
    return InputData

def readmRNAData(Cancer_type, Target, genders, groups, DataPath):
    
    mRNAData = loadmat(mRNADataPath)
    X, Y, GeneName, SampleName = mRNAData['X'].astype('float32'), mRNAData['CancerType'], mRNAData['FeatureName'][0], mRNAData['SampleName']
    GeneName = [row[0] for row in GeneName]
    SampleName = [row[0][0] for row in SampleName]
    Y = [row[0][0] for row in Y]
    
    mRNAData_X = pd.DataFrame(X, columns=GeneName, index=SampleName)
    mRNAData_Y = pd.DataFrame(Y, index=SampleName, columns=['Disease'])
    mRNAData_Y = mRNAData_Y[mRNAData_Y['Disease'].isin(tumor_types(Cancer_type))]
    mRNAData_in = mRNAData_X.join(mRNAData_Y, how='inner').drop(columns=['Disease'])
    
    mRNAData_in.index = [row[:12] for row in mRNAData_in.index.values]
    mRNAData_in = mRNAData_in.reset_index().drop_duplicates(subset='index', keep='first').set_index('index')
    
    # if criteria in ["GenderRace", "Race"]:
    GAData = [
        pd.read_excel(GADataPath, disease, usecols='A,E', index_col='Patient_ID', keep_default_na=False)
        for disease in tumor_types(Cancer_type)
    ]
    GAData_race = pd.concat(GAData)
    GAData_race = GAData_race[GAData_race['EIGENSTRAT'].isin(['EA', 'AA', 'EAA', 'NA', 'OA'])]
    GAData_race['race'] = GAData_race['EIGENSTRAT']
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'EA', 'race'] = 'WHITE'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'AA', 'race'] = 'BLACK'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'EAA', 'race'] = 'ASIAN'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'NA', 'race'] = 'NAT_A'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'OA', 'race'] = 'OTHER'
    GAData_race = GAData_race.drop(columns=['EIGENSTRAT'])
    if groups:
        GAData_race = GAData_race[GAData_race['race'].isin(groups)]

    if Target == 'OS':
        cols = 'B,E,H,I'
    elif Target == 'DSS':
        cols = 'B,E,J,K'
    elif Target == 'DFI':
        cols = 'B,E,L,M'
    elif Target == 'PFI':
        cols = 'B,E,N,O'
    OutcomeData = pd.read_excel(OutcomeDataPath, 'TCGA-CDR', usecols=cols, index_col='bcr_patient_barcode')
    
    OutcomeData = OutcomeData[OutcomeData['gender'].isin(genders)]
    OutcomeData.columns = ['G', 'E', 'T']
    OutcomeData = OutcomeData[OutcomeData['E'].isin([0, 1])].dropna()
    OutcomeData['C'] = 1 - OutcomeData['E']
    OutcomeData.drop(columns=['E'], inplace=True)
    mRNAData_in = mRNAData_in.join(GAData_race, how='inner').dropna(axis='columns')
    mRNAData_in = mRNAData_in.join(OutcomeData, how='inner')
    Data = mRNAData_in
    C = Data['C'].tolist()
    R = Data['race'].tolist()
    G = Data['G'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'race', 'T', 'G'])
    X = Data.values.astype('float32')
    
    InputData = {'X': X, 
                 'T': np.asarray(T, dtype=np.float32), 
                 'C': np.asarray(C, dtype=np.int32), 
                 'E': np.asarray(E, dtype=np.int32), 
                 'R': np.asarray(R), 
                 'G': np.asarray(G), 
                 'Samples': Data.index.values, 
                 'FeatureName': list(Data)}
    
    return InputData

def readProteinData(Cancer_type, Target, groups, genders, DataPath):
    
    df = pd.read_csv(ProteinDataPath, sep='\t', index_col='SampleID').dropna(axis=1)
    tumorTypes = tumor_types(Cancer_type)
    df = df[df['TumorType'].isin(tumorTypes)].drop(columns=['TumorType'])
    df.index = [row[:12] for row in df.index]
    ProteinData_in = df

    # if criteria in ["GenderRace", "Race"]:
    GAData = [
        pd.read_excel(GADataPath, disease, usecols='A,E', index_col='Patient_ID', keep_default_na=False)
        for disease in tumor_types(Cancer_type)
    ]
    GAData_race = pd.concat(GAData)
    GAData_race = GAData_race[GAData_race['EIGENSTRAT'].isin(['EA', 'AA', 'EAA', 'NA', 'OA'])]
    GAData_race['race'] = GAData_race['EIGENSTRAT']
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'EA', 'race'] = 'WHITE'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'AA', 'race'] = 'BLACK'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'EAA', 'race'] = 'ASIAN'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'NA', 'race'] = 'NAT_A'
    GAData_race.loc[GAData_race['EIGENSTRAT'] == 'OA', 'race'] = 'OTHER'
    GAData_race = GAData_race.drop(columns=['EIGENSTRAT'])
    if groups:
        GAData_race = GAData_race[GAData_race['race'].isin(groups)]

    if Target == 'OS':
        cols = 'B,E,H,I'
    elif Target == 'DSS':
        cols = 'B,E,J,K'
    elif Target == 'DFI':
        cols = 'B,E,L,M'
    elif Target == 'PFI':
        cols = 'B,E,N,O'
    OutcomeData = pd.read_excel(OutcomeDataPath, 'TCGA-CDR', usecols=cols, index_col='bcr_patient_barcode')
    
    OutcomeData = OutcomeData[OutcomeData['gender'].isin(genders)]
    OutcomeData.columns = ['G', 'E', 'T']
    OutcomeData = OutcomeData[OutcomeData['E'].isin([0, 1])].dropna()
    OutcomeData['C'] = 1 - OutcomeData['E']
    OutcomeData.drop(columns=['E'], inplace=True)
    ProteinData_in = ProteinData_in.join(GAData_race, how='inner').dropna(axis='columns')
    ProteinData_in = ProteinData_in.join(OutcomeData, how='inner')
    Data = ProteinData_in
    C = Data['C'].tolist()
    R = Data['race'].tolist()
    G = Data['G'].tolist()
    T = Data['T'].tolist()
    E = [1 - c for c in C]
    Data = Data.drop(columns=['C', 'race', 'T', 'G'])
    X = Data.values.astype('float32')
    
    InputData = {'X': X, 
                 'T': np.asarray(T, dtype=np.float32), 
                 'C': np.asarray(C, dtype=np.int32), 
                 'E': np.asarray(E, dtype=np.int32), 
                 'R': np.asarray(R), 
                 'G': np.asarray(G), 
                 'Samples': Data.index.values, 
                 'FeatureName': list(Data)}
    
    return InputData

def run_cv(Cancer_type, Feature_type, Target, groups, genders, DataPath):
    
    if Feature_type == 'mRNA':
        return readmRNAData(Cancer_type, Target, genders=genders, groups=groups, DataPath=DataPath)
    elif Feature_type == 'Protein':
        return readProteinData(Cancer_type, Target, genders=genders, groups=groups, DataPath=DataPath)
    elif Feature_type == 'Methylation':
        return readMethylationData(Cancer_type, Target, genders=genders, groups=groups, DataPath=DataPath)
    elif Feature_type == 'MicroRNA':
        return readMicroRNAData(Cancer_type, Target, genders=genders, groups=groups, DataPath=DataPath)
    else:
        raise ValueError(f"Unsupported Feature_type: {Feature_type}")

def merge_datasets(datasets):
    
    data = datasets[0]
    X, T, C, E, Samples, FeatureName = data['X'], data['T'], data['C'], data['E'], data['Samples'], data['FeatureName']
    
    df = pd.DataFrame(X, index=Samples, columns=FeatureName)
    df['T'] = T
    df['C'] = C
    df['E'] = E
    R, G = data['R'], data['G']
    df['R'] = R
    df['G'] = G

    for i in range(1, len(datasets)):
        data1 = datasets[i]
        X1, Samples, FeatureName = data1['X'], data1['Samples'], data1['FeatureName']
        temp = pd.DataFrame(X1, index=Samples, columns=FeatureName)
        df = df.join(temp, how='inner')

    C = df['C'].tolist()
    T = df['T'].tolist()
    E = df['E'].tolist()
    R = df['R'].tolist()
    G = df['G'].tolist()

    columns_to_drop = ['C', 'T', 'E']
    columns_to_drop += ['R', 'G']
    
    df = df.drop(columns=columns_to_drop)
    X = df.values.astype('float32')
    
    data = {
        'X': X,
        'T': np.asarray(T, dtype=np.float32),
        'C': np.asarray(C, dtype=np.int32),
        'E': np.asarray(E, dtype=np.int32),
        'R': np.asarray(R),
        'G': np.asarray(G),
        'Samples': df.index.values,
        'FeatureName': list(df)
    }
    
    return data

def run_cv_combined(Cancer_type, Feature_type, Target, groups, genders, DataPath):
    
    datasets = []
    for feature in Feature_type:
        if feature == 'Protein':
            datasets.append(readProteinData(Cancer_type, Target, genders=genders, groups=groups, DataPath=DataPath))
        elif feature == 'mRNA':
            datasets.append(readmRNAData(Cancer_type, Target, genders=genders, groups=groups, DataPath=DataPath))
        elif feature == 'MicroRNA':
            datasets.append(readMicroRNAData(Cancer_type, Target, genders=genders, groups=groups, DataPath=DataPath))
        elif feature == 'Methylation':
            datasets.append(readMethylationData(Cancer_type, Target, genders=genders, groups=groups, DataPath=DataPath))

    return merge_datasets(datasets)

def process_dataset(dataset, Years, count_category, groups, genders, 
                    Cancer_type, Feature_type, Target, min_cases=5, 
                    pos_col_s_ancest=None, neg_col_s_ancest=None, 
                    pos_col_t_ancest=None, neg_col_t_ancest=None, 
                    pos_col_s_gender=None, neg_col_s_gender=None, 
                    pos_col_t_gender=None, neg_col_t_gender=None):
    
    # print("Entered into process dataset loop")
    
    count_out = pd.DataFrame()
    
    dataset = get_n_years(dataset, Years)
    X, Y, R, y_strat, G, Gy_strat, GRy_strat = dataset
    df = pd.DataFrame(y_strat, columns=['RY'])
    df['GRY'] = GRy_strat
    df['GY'] = Gy_strat
    df['R'] = R
    df['G'] = G
    df['Y'] = Y
    print(X.shape)
    print(df['GRY'].value_counts())
    print(df['GY'].value_counts())
    print(df['G'].value_counts())
    print(df['RY'].value_counts())
    print(df['R'].value_counts())
    print(df['Y'].value_counts())

    prog_count = get_prognosis_count(dataset, count_category, groups, genders)

    if count_category == 'Race':
        
        if all(pc >= min_cases for pc in prog_count):
            count_out = pd.DataFrame({
                "Cancer Type": [Cancer_type],
                "Omics Feature": [Feature_type],
                "Clinical Outcome Endpoint": [Target],
                "Event Time Threshold (in Years)": Years,
                groups[0]: prog_count[0] + prog_count[2],
                groups[1]: prog_count[1] + prog_count[3],
                'Y(0)': prog_count[1] + prog_count[3],
                'Y(1)': prog_count[0] + prog_count[2],
                pos_col_s_ancest: prog_count[0],
                neg_col_s_ancest: prog_count[1],
                pos_col_t_ancest: prog_count[2],
                neg_col_t_ancest: prog_count[3]
            })

    elif count_category == 'Gender':
        
        if all(pc >= min_cases for pc in prog_count):
            count_out = pd.DataFrame({
                "Cancer Type": [Cancer_type],
                "Omics Feature": [Feature_type],
                "Clinical Outcome Endpoint": [Target],
                "Event Time Threshold (in Years)": Years,
                genders[0]: prog_count[0] + prog_count[1],
                genders[1]: prog_count[2] + prog_count[3],
                'Y(0)': prog_count[1] + prog_count[3],
                'Y(1)': prog_count[0] + prog_count[2],
                pos_col_s_gender: prog_count[0],
                neg_col_s_gender: prog_count[1],
                pos_col_t_gender: prog_count[2],
                neg_col_t_gender: prog_count[3]
            })

    elif count_category == 'GenderRace':
        
        pos_source_male = genders[0] + '1' + groups[0]
        neg_source_male = genders[0] + '0' + groups[0]
        pos_target_male = genders[0] + '1' + groups[1]
        neg_target_male = genders[0] + '0' + groups[1]
        pos_source_female = genders[1] + '1' + groups[0]
        neg_source_female = genders[1] + '0' + groups[0]
        pos_target_female = genders[1] + '1' + groups[1]
        neg_target_female = genders[1] + '0' + groups[1]
        
        if all(pc >= min_cases for pc in prog_count):
            count_out = pd.DataFrame({
                "Cancer Type": [Cancer_type],
                "Omics Feature": [Feature_type],
                "Clinical Outcome Endpoint": [Target],
                "Event Time Threshold (in Years)": Years,
                genders[0]: prog_count[0] + prog_count[1] + prog_count[2] + prog_count[3],
                genders[1]: prog_count[4] + prog_count[5] + prog_count[6] + prog_count[7],
                groups[0]: prog_count[0] + prog_count[1] + prog_count[4] + prog_count[5],
                groups[1]: prog_count[2] + prog_count[3] + prog_count[6] + prog_count[7],
                'Y(0)': prog_count[1] + prog_count[3] + prog_count[5] + prog_count[7],
                'Y(1)': prog_count[0] + prog_count[2] + prog_count[4] + prog_count[6],
                pos_col_s_ancest: prog_count[0] + prog_count[4],
                neg_col_s_ancest: prog_count[1] + prog_count[5],
                pos_col_t_ancest: prog_count[2] + prog_count[6],
                neg_col_t_ancest: prog_count[3] + prog_count[7],
                pos_col_s_gender: prog_count[0] + prog_count[2],
                neg_col_s_gender: prog_count[1] + prog_count[3],
                pos_col_t_gender: prog_count[4] + prog_count[6],
                neg_col_t_gender: prog_count[5] + prog_count[7],
                pos_source_male: prog_count[0],
                neg_source_male: prog_count[1],
                pos_target_male: prog_count[2],
                neg_target_male: prog_count[3],
                pos_source_female: prog_count[4],
                neg_source_female: prog_count[5],
                pos_target_female: prog_count[6],
                neg_target_female: prog_count[7]
            })
    else:
        count_out = None

    return count_out

def main():
    
    parser = argparse.ArgumentParser(description="Process counting arguments.")

    parser.add_argument(
        "--count_category", type=str, choices=["GenderRace", "Race", "Gender"], required=True,
        help="Specify the count category: GenderRace, Race, Gender"
    )
    parser.add_argument(
        "--omicsConfiguration", type=str, choices=["single", "combination"], required=True,
        help="Specify if the data is single or combination of omics feature from the TCGA dataset"
    )
    parser.add_argument(
        "--num_features", type=int, nargs='?', default=None,
        help="Number of features to be combined (2,3,4). Required if omicsConfiguration is combination"
    )
    parser.add_argument(
        "--min_cases", type=int, required=True, default=5,
        help="Minimum number of cases. Normally, it is 5"
    )
    parser.add_argument(
        "--DDP_group", type=str, nargs='?', default=None,
        help="Specify DDP group: BLACK, ASIAN, NAT_A. Required if count_category is GenderRace or Race"
    )
    
    args = parser.parse_args()
    
    print(args, flush=True)

    # Conditional checks
    if args.count_category in ["GenderRace", "Race"] and args.DDP_group is None:
        parser.error("DDP_group is required when count_category is GenderRace or Race.")
    if args.count_category == "Gender" and args.DDP_group is not None:
        parser.error("DDP_group is not required when count_category is Gender.")
    if args.omicsConfiguration == "combination" and args.num_features is None:
        parser.error("num_features is required when omicsConfiguration is combination.")
    if args.omicsConfiguration == "single" and args.num_features is not None:
        parser.error("num_features is not required when omicsConfiguration is single.")

    print(args)

    count_category = args.count_category
    omicsConfiguration = args.omicsConfiguration
    num_features = args.num_features
    min_cases = args.min_cases
    DDP_group = args.DDP_group
    groups = ("WHITE", DDP_group) if DDP_group else ("WHITE","BLACK","ASIAN","NAT_A")
    genders = ("MALE", "FEMALE")
    
    # Debug print to verify correct setup
    print(f"Which count: {count_category}")
    print(f"Count type: {omicsConfiguration}")
    print(f"Number of features: {num_features}")
    print(f"Minimum cases: {min_cases}")
    print(f"DDP group: {DDP_group}")
    print(f"Groups: {groups}")
    print(f"Genders: {genders}")
    
    pos_col_s_ancest = '1' + groups[0]
    neg_col_s_ancest = '0' + groups[0]
    pos_col_t_ancest = '1' + groups[1] if DDP_group else None
    neg_col_t_ancest = '0' + groups[1] if DDP_group else None
    pos_col_s_gender = '1' + genders[0]
    neg_col_s_gender = '0' + genders[0]
    pos_col_t_gender = '1' + genders[1]
    neg_col_t_gender = '0' + genders[1]
    
    if omicsConfiguration == 'single':
        for Feature_type in Features:
            count_out_data = pd.DataFrame()
            for Target in Endpoints:
                for Cancer_type in Cancers:
                    for Years in YearsAll:
                        print(f'=================================\n{Cancer_type} {Feature_type} {Target} {Years}\n=================================')
                        
                        dataset = run_cv(Cancer_type, Feature_type, Target, groups, genders, DataPath)
                        
                        print(dataset['X'].shape[0])
                        if dataset['X'].shape[0] != 0:
                            count_out = process_dataset(dataset, Years, count_category, groups, genders, 
                                                        Cancer_type, Feature_type, Target, min_cases=min_cases, 
                                                        pos_col_s_ancest=pos_col_s_ancest, neg_col_s_ancest=neg_col_s_ancest, 
                                                        pos_col_t_ancest=pos_col_t_ancest, neg_col_t_ancest=neg_col_t_ancest, 
                                                        pos_col_s_gender=pos_col_s_gender, neg_col_s_gender=neg_col_s_gender, 
                                                        pos_col_t_gender=pos_col_t_gender, neg_col_t_gender=neg_col_t_gender)
                            count_out_data = pd.concat([count_out_data, count_out], ignore_index=True)
                            
            filename = os.path.join(folderISAAC, f'TCGA-{count_category}_count_{Feature_type}.xlsx')
            if DDP_group:
                filename = os.path.join(folderISAAC, f'TCGA-{count_category}_count_{Feature_type}_{DDP_group}.xlsx')
            count_out_data.to_excel(filename)
    
    elif omicsConfiguration == 'combination':
        comb_feature = combinations(Features, num_features)
        for Feature_type in comb_feature:
            combName = '_'.join(Feature_type)
            count_out_data = pd.DataFrame()
            for Target in Endpoints:
                for Cancer_type in Cancers:
                    for Years in YearsAll:
                        print(f'=================================\n{Cancer_type} {Target} {Years}\n=================================')
    
                        dataset = run_cv_combined(Cancer_type, Feature_type, Target, groups, genders, DataPath)
    
                        if dataset['X'].shape[0] != 0:
                            count_out = process_dataset(dataset, Years, count_category, groups, genders, 
                                                        Cancer_type, Feature_type, Target, min_cases=min_cases, 
                                                        pos_col_s_ancest=pos_col_s_ancest, neg_col_s_ancest=neg_col_s_ancest, 
                                                        pos_col_t_ancest=pos_col_t_ancest, neg_col_t_ancest=neg_col_t_ancest, 
                                                        pos_col_s_gender=pos_col_s_gender, neg_col_s_gender=neg_col_s_gender, 
                                                        pos_col_t_gender=pos_col_t_gender, neg_col_t_gender=neg_col_t_gender)
                            count_out_data = pd.concat([count_out_data, count_out], ignore_index=True)
            
            filename = os.path.join(folderISAAC, f'TCGA-{count_category}_count_{combName}.xlsx')
            if DDP_group:
                filename = os.path.join(folderISAAC, f'TCGA-{count_category}_count_{combName}_{DDP_group}.xlsx')
            count_out_data.to_excel(filename)

if __name__ == '__main__':
    main()

