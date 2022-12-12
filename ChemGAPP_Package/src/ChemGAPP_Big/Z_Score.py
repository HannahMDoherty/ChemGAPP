#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
parser = argparse.ArgumentParser(description=" Compares each replicate colony to find outliers within colony size for each plate.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--InputFile", help="The normalised csv file from Check_Normalisation.py")
parser.add_argument("-o", "--OutputFile", help="A CSV file of the dataset where colony sizes are replaced with the colony type values.")
args = vars(parser.parse_args())
inputfile1 = args["InputFile"]
outputfile1 = args["OutputFile"]

def Z_Score(inputfile,outputfile):
    import pandas as pd
    import numpy as np
    import scipy.stats as stats
    def z_score_method(df):
        z = np.abs(stats.zscore(df))
        threshold = 1.4
        outlier = []
        index=0
        for i, v in enumerate(z):
            if np.any(v) == 'nan':
                continue
            else:
                if v > threshold:
                    outlier.append(i)
                else:
                    continue
            return outlier
    m = pd.read_csv(inputfile,
                          index_col=[0, 1],
                          header=[0, 1, 2, 3])
    conditions = {x[0:2] for x in m.columns}
    rounds = 0
    ABSN = pd.DataFrame(index=m.index)
    for c in sorted(conditions):
        rounds = rounds + 1
        print(rounds)
        df1 = m.xs((c), axis =1, drop_level=False)
        ar1=np.array(df1)
        ar2=np.full(np.shape(ar1), "N")
        for j in range(len(ar1)):
            #if np.all(ar1[j]) == 0:
            #    for l, vs in enumerate(ar1[j]):
            #        if vs == 0:
            #            ar2[j][l] = "Z"
            #if np.any(ar1[j]) == 0:
            #    for l, vs in enumerate(ar1[j]):
            #        if vs == 0:
            #            ar2[j][l] = "T"        
            for l, vs in enumerate(ar1[j]):
                if str(vs) == "nan":
                    ar2[j][l] = "X" 
            if np.all(ar1[j]) == 0:
                continue
            else: 
                outlier = z_score_method(ar1[j])
                if outlier != None:
                    if len(outlier) > 0:
                        mean = np.nanmean(ar1[j])
                        if ar1[j][outlier] < mean:
                            ar2[j][outlier] ="S"
                        else:
                            ar2[j][outlier]="B"
        ar_df = pd.DataFrame(ar2, index=m.index)
        ABSN = pd.concat([ABSN, ar_df], axis=1)
    ABSN.columns = (pd.MultiIndex.from_tuples(sorted(m.columns)))
    ABSN.to_csv(outputfile)
    return ABSN

Z_Score(inputfile1,outputfile1)