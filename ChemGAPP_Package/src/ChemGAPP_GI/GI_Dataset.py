#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import numpy as np
import pandas as pd
import re
import os
parser = argparse.ArgumentParser(description="Produces datasets with genetic interaction scores for genetic interaction screens",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--inputfiles", help="Path to IRIS files.")
parser.add_argument("-p", "--PATH", help="Path for the output files.")
parser.add_argument("-n", "--nameinfofiles", help="Path to plate information files. Plate info files should be txt files, with the columns: Row, Column, Strain, Replicate, Order, Set.")
args = vars(parser.parse_args())
indir1 = args["inputfiles"]
path1 = args["PATH"]
infodir1 = args["nameinfofiles"]

def GI_dataset(indir,PATH,infodir):    
    indir = os.path.expanduser(indir)                     
    m = None
    for f in os.listdir(indir):
        if f.endswith(".iris"): 
            #sys.stderr.write(batch + ' ' + f + '\n')
            #print(os.path.join(indir, batch, f))
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
                #sets them in a dataframe grouped by the condition and then columns as the A B C D or E etc
                m = m.join(m1, how='inner')
    def plate_info(file):
        p = pd.read_table(file)
        if len(p.columns) != 6:
            print("ERROR: Info file formatting incorrect! Should be 6 columns: Row,Column,Strain,Replicate,Order,Set")
        if len(p.columns) == 6:
            p.columns = ['row','column','strain','Replicate','Order','Set']
            p = p.set_index(['row','column'])
            p = p.sort_index()
        return p
            
    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        return [ atoi(c) for c in re.split(r'(\d+)', text) ]
                
    #p = plate_info(ident)
    #df2 = pd.merge(m,p,left_index=False, right_index=True)
    #df2 = df2.rename(columns={'strain': 'Gene'})
    #df2 = df2.set_index(['Set','Replicate','Order','Gene'])
    #df2.columns = df2.columns.str.split('(CIP.+)', expand=True).droplevel(2)
    #st.write(df2)
    infodir = os.path.expanduser(infodir)
    files = []
    for filename in os.listdir(infodir):
        f = os.path.join(infodir, filename)
        # checking if it is a file
        if f.endswith(".txt"):
            if os.path.isfile(f):
                #print(f)
                files.append(f)
    files.sort(key=natural_keys)
    plate_DFs = []
    for i in files:
        p = plate_info(i)
        plate_DFs.append(p)
    plates = {x[0:2] for x in m.columns}
    m2 = pd.DataFrame(m.copy(deep=False))
    columns1 = [x[0] for x in m2.columns]
    df_with_strains = pd.DataFrame(columns = [columns1[0],'strain','Replicate','Order','Set'])
    for a, n in zip(plate_DFs, sorted(plates)):
        df1 = (m.xs((n), axis =1, drop_level=True))
        df1 = pd.DataFrame(df1)
        df2 = pd.merge(df1,a,left_index=True, right_index=True)
        #df2 = df2.set_index(['strain','Replicate','Order','Set'])
        df2.columns = [n[0],'strain','Replicate','Order','Set']
        df_with_strains = pd.concat([df_with_strains, df2], ignore_index=True)
    df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
    df2 = df_with_strains.set_index(['Set','Replicate','Order','Gene'])
    set = {x[0] for x in df2.index}
    for s in sorted(set):
        df = pd.DataFrame()
        df1 = df2.xs((s), axis=0, drop_level=True)
        replicate = {x[0] for x in df1.index}
        for r in sorted(replicate):
            dfb = df1.xs((r), axis=0, drop_level=False)
            dfb = dfb.sort_index(level="Order")
            print(dfb.iloc[2].name[2],(dfb.iloc[1][0]/dfb.iloc[0][0],dfb.iloc[2][0]/dfb.iloc[0][0],(dfb.iloc[1][0]/dfb.iloc[0][0])*(dfb.iloc[2][0]/dfb.iloc[0][0]),dfb.iloc[3][0]/dfb.iloc[0][0]))
            name = (dfb.iloc[1][0]/dfb.iloc[0][0],dfb.iloc[2][0]/dfb.iloc[0][0],(dfb.iloc[1][0]/dfb.iloc[0][0])*(dfb.iloc[2][0]/dfb.iloc[0][0]),dfb.iloc[3][0]/dfb.iloc[0][0])
            columns = list([dfb.iloc[1].name[2],"Secondary Gene","Double Expected","Double Observed"])
            data = []
            zipped = zip(columns, name)
            a_dictionary = dict(zipped)
            data.append(a_dictionary)
            df = df.append(data, True)
        df1.to_csv((PATH+"/"+dfb.iloc[2].name[2]+"_Colony_sizes.csv"))
        df.to_csv((PATH+"/"+dfb.iloc[2].name[2]+"_Interaction_Scores.csv"))
GI_dataset(indir1,path1,infodir1)
