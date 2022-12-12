#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import os
import pandas as pd
import re
parser = argparse.ArgumentParser(description="Add the gene names from the plate info files to make the final dataset. The plate info files must be in a folder by themselves and should be .txt files.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--InputFile", help="The CSV output of S_Scores.py")
parser.add_argument("-o", "--Outputpath", help="A path a prefix for the output files E.g ~/Desktop/ChemGAPP would make ChemGAPP_Final_Dataset.txt .")
parser.add_argument("-p", "--PATH", help="The path to the folder containing the plate info files.")
args = vars(parser.parse_args())
inputfile1 = args["InputFile"]
outputfile1 = args["Outputpath"]
PATH1 = args["PATH"]

def Add_gene_names(PATH,inputfile,outputpath):
    def plate_info(file):
        p = pd.read_table(file)
        if len(p.columns) == 5:
            p.columns = ['row','column','strain','info','normalisation_group']
            p = p.set_index(['row','column'])
            p = p.drop(['info','normalisation_group'], axis=1)
            p = p.sort_index()
        if len(p.columns) == 3:
            p.columns = ['row','column','strain']
            p = p.set_index(['row','column'])
            p = p.sort_index()
        return p
    # assign directory
    directory = os.path.expanduser(PATH)
    # iterate over files in
    # that directory
    files = []
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if f.endswith(".txt"):
            if os.path.isfile(f):
                #print(f)
                files.append(f)

    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        return [ atoi(c) for c in re.split(r'(\d+)', text) ]

    files.sort(key=natural_keys)
    plate_DFs = []
    for i in files:
        p = plate_info(i)
        plate_DFs.append(p)
    nm3 = pd.read_csv(inputfile,index_col=[0, 1], header=[0, 1])
    plates = {x[0] for x in nm3.columns}
    nm4 = nm3.copy(deep=False)
    nm4['1','strain'] = 'p'
    columns1 = {x[1] for x in nm4.columns}
    df_with_strains = pd.DataFrame(columns = sorted(columns1))
    for a, n in zip(plate_DFs, sorted(plates)):
        df1 = (nm3.xs((n), axis =1, drop_level=True))
        df2 = pd.merge(df1,a,left_index=True, right_index=True)
        df_with_strains = pd.concat([df_with_strains , df2], ignore_index=True)
    df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
    df_with_strains = df_with_strains.set_index('Gene')

    groupmean = df_with_strains.groupby(level=0).mean()
    df_averaged = pd.DataFrame(groupmean, index=groupmean.index, columns=groupmean.columns)
    df_averaged.index.name = None
    outputfile = os.path.expanduser(outputpath)
    df_with_strains.to_csv(outputfile+"_Final_Dataset.txt", sep='\t')
    df_averaged.to_csv(outputfile+"_Final_Dataset_Averaged.txt", sep='\t')
    return df_with_strains

Add_gene_names(PATH1,inputfile1,outputfile1)