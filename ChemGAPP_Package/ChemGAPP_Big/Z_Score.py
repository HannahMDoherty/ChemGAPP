#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import pandas as pd
import numpy as np
import scipy.stats as stats

def get_options():
    parser = argparse.ArgumentParser(description=" Compares each replicate colony to find outliers within colony size for each plate.",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--InputFile", help="The normalised csv file from Check_Normalisation.py")
    parser.add_argument("-o", "--OutputFile", help="A CSV file of the dataset where colony sizes are replaced with the colony type values.")
    return parser.parse_args()

def main():
    options = get_options()
    inputfile = options.InputFile
    outputfile = options.OutputFile
    def z_score_method(df):
        #performs z-score test and between replicates of same mutant in 
        # same condition plate and takes the absolute value of the z-score value
        z = np.abs(stats.zscore(df))
        threshold = 1.4
        outlier = []
        #goes through each z-score value and position  
        for i, v in enumerate(z):
            #if nan then continue
            if np.any(v) == 'nan':
                continue
            else:
                #if the z-score value is greater than the threshold then 
                # it is an outlier and the postion is appended to 'outlier' list
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
        #makes array that matches the colony size data, however instead of values is all N's for Normal
        ar2=np.full(np.shape(ar1), "N")
        #iterates through row of the columns array.
        for j in range(len(ar1)):
            #iterates through the position and value of row     
            for l, vs in enumerate(ar1[j]):
                #if the value is nan, the matching position in the N array is changed to an X
                if str(vs) == "nan":
                    ar2[j][l] = "X" 
            #if all the replicates in the row are 0 then the row is skipped 
            if np.all(ar1[j]) == 0:
                continue
            else: 
                #performs the z-score method for the replicates and makes list of the outliers
                outlier = z_score_method(ar1[j])
                #if there are outliers then it will find the mean of the replicates and then look 
                # at the positions specified as outliers and see if they are larger than the mean
                #  or smaller and designated an B or S accordingly in the N's array.
                if outlier != None:
                    if len(outlier) > 0:
                        mean = np.nanmean(ar1[j])
                        if ar1[j][outlier] < mean:
                            ar2[j][outlier] ="S"
                        else:
                            ar2[j][outlier]="B"
        #the array of Ns, Bs, Ss and Xs are formatted and saved as a df.
        ar_df = pd.DataFrame(ar2, index=m.index)
        ABSN = pd.concat([ABSN, ar_df], axis=1)
    ABSN.columns = (pd.MultiIndex.from_tuples(sorted(m.columns)))
    ABSN.to_csv(outputfile)
    return ABSN

if __name__ == "__main__":
    main()