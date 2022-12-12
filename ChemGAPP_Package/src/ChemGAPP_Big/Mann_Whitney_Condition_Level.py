#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
parser = argparse.ArgumentParser(description="The variance of the replicate means for each condition is calculated and then the average of these variance is calculated for each plate within that conditions, i.e the variance between replicate plate A,B,C,D for plate 1 of condition A, and then the average of plate 1, 2, 3 etc. for condition A."
                                 ,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--InputFile", help="The CSV file with the mean u statistics and p values for each replicate from Mann_Whitney_Plate_Level.py")
parser.add_argument("-o", "--OutputFile", help="A CSV file with the mean variance values for the u statistic and p values of the mann-whitney test for each condition.")
args = vars(parser.parse_args())
inputfile1 = args["InputFile"]
outputfile1 = args["OutputFile"]

def Mann_Whitney_Condition_Level(inputfile,outputfile):
    import pandas as pd
    import numpy as np
    Pmean = pd.read_csv(inputfile)
    Pmean = Pmean.set_index(['Condition','Plate','Batch','Replicate'])
    cond2 = {x[0:3] for x in Pmean.index}
    P_var = pd.DataFrame(columns=['Plate','Condition','Batch','Variance U-Stat','Variance P-Value'])
    for c2 in sorted(cond2):
        df1= Pmean.xs((c2), axis =0, drop_level=False)
        ar1 = np.array(df1)
        pvar_val = np.nanvar(ar1[:,0])
        uvar_val = np.nanvar(ar1[:,1])
        name = (df1.index[0][1],df1.index[0][0],df1.index[0][2],pvar_val,uvar_val)
        columns = list(P_var)
        data = []
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        P_var = P_var.append(data, True)
    P_var = P_var.set_index(['Condition','Batch','Plate'])
    P_var_cond = P_var.groupby(level=[0,1]).mean()
    P_var_cond = P_var_cond.reset_index()
    P_var_cond = P_var_cond.rename(columns={'Variance P-Value':'Mean Variance P-Value','Variance U-Stat':'Mean Variance U-Stat'})
    
    P_var_cond.to_csv(outputfile, index = False)
    return P_var_cond

Mann_Whitney_Condition_Level(inputfile1,outputfile1)