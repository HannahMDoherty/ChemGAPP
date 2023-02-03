#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import os
import pandas as pd
import re

def get_options():

    parser = argparse.ArgumentParser(description="Add the gene names from the plate info files to make the final dataset. The plate info files must be in a folder by themselves and should be .txt files.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--InputFile", help="The CSV output of S_Scores.py")
    parser.add_argument("-o", "--Outputpath", help="A path a prefix for the output files E.g ~/Desktop/ChemGAPP would make ChemGAPP_Final_Dataset.txt .")
    parser.add_argument("-p", "--PATH", help="The path to the folder containing the plate info files.")
    return parser.parse_args()

def main():
    options = get_options()
    PATH = options.PATH
    inputfile = options.InputFile
    outputpath = options.Outputpath
    
    def plate_info(file):
        p = pd.read_table(file)
        #renames the plate info file columns and adds row and column to index before sorting by index.
        if len(p.columns) == 3:
            p.columns = ['row','column','strain']
            p = p.set_index(['row','column'])
            p = p.sort_index()
        else:
            print("Plate information file" + str(file) + "not in correct format.")
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
    
    # looks for digits within the file names so that order of 
    # plate info files is sorted plate1, plate2, plate3 etc. not plate1, plate10, plate2 etc.
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
    #names columns with condition
    df_with_strains = pd.DataFrame(columns = sorted(columns1))
    # splits by each source plate and then assigns gene names from matching plate information file.
    for a, n in zip(plate_DFs, sorted(plates)):
        df1 = (nm3.xs((n), axis =1, drop_level=True))
        df2 = pd.merge(df1,a,left_index=True, right_index=True)
        df_with_strains = pd.concat([df_with_strains , df2], ignore_index=True)
    df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
    df_with_strains = df_with_strains.set_index('Gene')

    #averages scores for rows with the same gene name.
    groupmean = df_with_strains.groupby(level=0).mean()
    df_averaged = pd.DataFrame(groupmean, index=groupmean.index, columns=groupmean.columns)
    df_averaged.index.name = None
    outputfile = os.path.expanduser(outputpath)
    df_with_strains.to_csv(outputfile+"_Final_Dataset.txt", sep='\t')
    df_averaged.to_csv(outputfile+"_Final_Dataset_Averaged.txt", sep='\t')
    return df_with_strains

if __name__ == "__main__":
    main()