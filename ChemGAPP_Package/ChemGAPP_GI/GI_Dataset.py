#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import numpy as np
import pandas as pd
import re
import os
def get_options():
    parser = argparse.ArgumentParser(description="Produces datasets with genetic interaction scores for genetic interaction screens",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--inputfiles", help="Path to IRIS files.")
    parser.add_argument("-p", "--PATH", help="Path for the output files.")
    parser.add_argument("-n", "--nameinfofiles", help="Path to plate information files. Plate info files should be txt files, with the columns: Row, Column, Strain, Replicate, Order, Set.")
    return parser.parse_args()

def GI_dataset():    
    options = get_options()
    indir = options.inputfiles
    PATH = options.PATH
    infodir = options.nameinfofiles
    indir = os.path.expanduser(indir)                     
    m = None
    # cycles through iris files and uses filename to produce column headers.
    for f in os.listdir(indir):
        if f.endswith(".iris"):
            g = pd.read_csv(os.path.join(indir, f),
                    comment='#',
                    index_col=[0, 1],
                    sep='\t')
            if m is None:
                try:
                    m = g['colony size']
                except:
                    m = g['size']
                    m.name = (f.split('_')[0],f.split('.')[0].split('_')[1])
                    m = m.to_frame()
            else:
                try:
                    m1 = g['colony size']
                except:
                    m1 = g['size']
                    m1.name = (f.split('_')[0],
                            f.split('.')[0].split('_')[1])
                    m1 = m1.to_frame()
                #sets them in a dataframe grouped by the secondary gene name and then replicates.
                m = m.join(m1, how='inner')
    
    # renames the plate info file columns and adds row and column to index before sorting by index.
    def plate_info(file):
        p = pd.read_table(file)
        if len(p.columns) != 6:
            print("ERROR: Info file formatting incorrect! Should be 6 columns: Row,Column,Strain,Replicate,Order,Set")
        if len(p.columns) == 6:
            p.columns = ['row','column','strain','Replicate','Order','Set']
            p = p.set_index(['row','column'])
            p = p.sort_index()
        return p


    # looks for digits within the file names so that order of 
    # plate info files is sorted plate1, plate2, plate3 etc. not plate1, plate10, plate2 etc.        
    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        return [ atoi(c) for c in re.split(r'(\d+)', text) ]

    # iterate over plate files in that directory           
    infodir = os.path.expanduser(infodir)
    files = []
    for filename in os.listdir(infodir):
        f = os.path.join(infodir, filename)
        # checking if it is a file
        if f.endswith(".txt"):
            if os.path.isfile(f):
                files.append(f)
    files.sort(key=natural_keys)
    
    #adds the plate files to a list
    plate_DFs = []
    for i in files:
        p = plate_info(i)
        plate_DFs.append(p)
    plates = {x[0:2] for x in m.columns}
    m2 = pd.DataFrame(m.copy(deep=False))
    # splits by secondary gene name 
    columns1 = [x[0] for x in m2.columns]
    df_with_strains = pd.DataFrame(columns = [columns1[0],'strain','Replicate','Order','Set'])
    # adds the gene names, replicate number, set and order to the colony size data based on row and column
    for a, n in zip(plate_DFs, sorted(plates)):
        df1 = (m.xs((n), axis =1, drop_level=True))
        df1 = pd.DataFrame(df1)
        df2 = pd.merge(df1,a,left_index=True, right_index=True)
        df2.columns = [n[0],'strain','Replicate','Order','Set']
        df_with_strains = pd.concat([df_with_strains, df2], ignore_index=True)
    df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
    df2 = df_with_strains.set_index(['Set','Replicate','Order','Gene'])
    set = {x[0] for x in df2.index}
    #splits by set, thus different gene sets have their own output files. 
    for s in sorted(set):
        df = pd.DataFrame()
        df1 = df2.xs((s), axis=0, drop_level=True)
        replicate = {x[0] for x in df1.index}
        for r in sorted(replicate):
            dfb = df1.xs((r), axis=0, drop_level=False)
            dfb = dfb.sort_index(level="Order")
            #calculates the fitness ratios and double expected ratios and adds them to the dataset
            name = (dfb.iloc[1][0]/dfb.iloc[0][0],dfb.iloc[2][0]/dfb.iloc[0][0],(dfb.iloc[1][0]/dfb.iloc[0][0])*(dfb.iloc[2][0]/dfb.iloc[0][0]),dfb.iloc[3][0]/dfb.iloc[0][0])
            columns = list([dfb.iloc[1].name[2],"Secondary Gene","Double Expected","Double Observed"])
            data = []
            zipped = zip(columns, name)
            a_dictionary = dict(zipped)
            data.append(a_dictionary)
            df = df.append(data, True)
        df1.to_csv((PATH+"/"+dfb.iloc[2].name[2]+"_Colony_sizes.csv"))
        df.to_csv((PATH+"/"+dfb.iloc[2].name[2]+"_Interaction_Scores.csv"))

if __name__ == "__main__":
    main()