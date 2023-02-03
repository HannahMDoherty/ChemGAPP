#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import pandas as pd
import numpy as np
import scipy.stats as stats

def get_options():
    parser = argparse.ArgumentParser(description="The variance of replicate colony sizes is calculated for each plate and these variance values are averaged for each plate within a condition.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--InputFile", help="The normalised csv file from Check_Normalisation.py")
    parser.add_argument("-o", "--OutputFile", help="A CSV file of the average variances for each condition.")
    return parser.parse_args()


def main():
    options = get_options()
    inputfile = options.InputFile
    outputfile = options.OutputFile
    m = pd.read_csv(inputfile,
                          index_col=[0, 1],
                          header=[0, 1, 2,3])
    m.columns = m.columns.swaplevel(2, 3)
    #makes df with same index as the input normalised dataset file.
    Var_DF = pd.DataFrame(index=m.index)
    conditions = {x[0:3] for x in m.columns}
    rounds = 0 
    #splits into source plate, batch and condition, then compares the variance between the replicates.
    for c in sorted(conditions):
        rounds = rounds + 1
        print(rounds)
        df1 = m.xs((c), axis =1, drop_level=False)
        ar1 = np.array(df1)
        ar2 = np.array([])
        for j in range(0,len(ar1)):
            #The variance of each row is calculated and added to the variance column 
            #if all values are nan then variance is set to nan.
            if np.count_nonzero(~np.isnan(ar1[j])) == 0:
                var = "nan"
            #otherwise calculates variance ingnoring nans
            else:
                var = np.nanvar(ar1[j])
            #appends variance to array of variances
            ar2 = np.append(ar2, var)
        #set array as df
        ar_df = pd.DataFrame(ar2, index=m.index)
        #appends array of variances to variance df
        Var_DF = pd.concat([Var_DF, ar_df], axis=1)
    #sets column names to the source plate, batch and condition.
    Var_DF.columns = (pd.MultiIndex.from_tuples(sorted(conditions)))
    ave_Var_plate = pd.DataFrame(columns=['Condition','Batch','Plate','Average Variance'])
    #calculates the average variance for each condition.
    for f in Var_DF.columns:
        name = (f[1],f[2],f[0],np.nanmean(Var_DF[f].values.astype(float)))
        columns = list(ave_Var_plate)
        data = []
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        ave_Var_plate = ave_Var_plate.append(data, True)
    ave_Var_plate = ave_Var_plate.set_index(['Condition',"Batch",'Plate'])
    cond3 = {x[0:2] for x in ave_Var_plate.index}
    #calculates mean across different source plates and produces df.
    ave_Var_cond = pd.DataFrame(columns=(['Condition','Batch','Average Variance']))
    for cd3 in sorted(cond3):
        dfVC = ave_Var_plate.xs(cd3, axis =0, drop_level=False)
        name = (cd3[0],cd3[1],dfVC['Average Variance'].mean())
        data = []
        columns = list(ave_Var_cond)
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        ave_Var_cond = ave_Var_cond.append(data, True)
    ave_Var_cond.to_csv(outputfile,index=False)
    return ave_Var_cond

if __name__ == "__main__":
    main()