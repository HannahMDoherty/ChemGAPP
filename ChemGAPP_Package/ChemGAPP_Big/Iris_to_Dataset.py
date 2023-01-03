#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
parser = argparse.ArgumentParser(description="Takes folder of Iris files and produces Dataset for ChemGAPP",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-p", "--PATH", help="Path to folder which contains IRIS files, IRIS file names should be in the format: CONDITION-concentration-platenumber-batchnumber_replicate.JPG.iris")
parser.add_argument("-o", "--outputfile", help="Name of output file, should be a .txt file")
parser.add_argument("-it", "--IRISphenotype", help="IRIS phenotype from IRIS files you wish to analyse",default='size')
args = vars(parser.parse_args())
PATH1 = args["PATH"]
outputfile1 = args["outputfile"]
phen1 = args["IRISphenotype"]

def Iris_to_dataset(PATH,outputfile,phen):    
    import os
    import sys
    import numpy as np
    import pandas as pd
    m = None
    # cycles through iris files and uses filename to produce column headers.
    indir = os.path.expanduser(PATH)
    for f in sorted(os.listdir(indir)):
        if f.endswith(".iris"): 
            sys.stderr.write(f + '\n')
            print(os.path.join(indir, f))
            g = pd.read_csv(os.path.join(indir, f),
                    comment='#',
                    index_col=[0, 1],
                    sep='\t')
            if m is None:
                try:
                    m = g[phen]
                except:
                    m = g[phen]
                m.name = (int(f.split('-')[2].split('_')[0]),
                          '-'.join(f.split('.')[0].split('-')[:2]).replace(",","."),
                          f.split('.')[0].split('_')[1],"Batch"+ f.split('-')[3].split('_')[0]
                          )
                m = m.to_frame()
            else:
                try:
                    m1 = g[phen]
                except:
                    m1 = g[phen]
                m1.name = (int(f.split('-')[2].split('_')[0]),
                          '-'.join(f.split('.')[0].split('-')[:2]).replace(",","."),
                          f.split('.')[0].split('_')[1],
                          "Batch"+ f.split('-')[3].split('_')[0])
                m1 = m1.to_frame()
            
                m = m.join(m1, how='inner')
    m.to_csv(outputfile)

Iris_to_dataset(PATH1,outputfile1,phen1)