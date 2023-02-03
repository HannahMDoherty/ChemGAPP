#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import pandas as pd
import numpy as np
import scipy.stats as stats

def get_options():
    parser = argparse.ArgumentParser(description="Compares the distributions of the colony sizes of replicate plates of the same condition and determines if replicate plates have the same distribution based on the p value of the Mann whitney test.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--InputFile", help="The normalised csv file from Check_Normalisation.py output")
    parser.add_argument("-o", "--OutputFile", help="A CSV file with the u statistics and p values for each comparison.")
    parser.add_argument("-o2", "--OutputFile_Mean", help="A CSV file with the mean u statistics and p values for each replicate")
    return parser.parse_args()


def main():
    options = get_options()
    inputfile = options.InputFile
    outputfile = options.OutputFile
    outputfile2 = options.OutputFile_Mean
    #opens the normalised dataset file
    m = pd.read_csv(inputfile,
                      index_col=[0, 1],
                      header=[0, 1, 2, 3])
    m.columns = m.columns.swaplevel(2,3)
    plates = {x[0:3] for x in m.columns}
    Mann_whit_all = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate 1',
                                      'Replicate 2','U-statistic', 'P-Value'])
    # splits and iterates by plate, condition.
    for p in sorted(plates):
        print(p)
        df3 = m.xs((p), axis =1, drop_level=False)
        print(df3)
        ar1 = np.array(df3)
        for i in list(range(0,len(df3.columns))):
            for j in list(range(0,len(df3.columns))):
                #compares columns that arent the same in pairs
                if i != j:
                    #takes the whole set of values then removes nan values and the corresponding values in the paired set
                    #This allows the mannwhitneyu to work
                    rep1 = ar1[:,i]
                    rep2 = ar1[:,j]
                    rep1 = rep1[~np.isnan(rep1)]
                    rep2 = rep2[~np.isnan(rep2)] 
                    u_statistic, p_value = stats.mannwhitneyu(rep1,rep2)
                    #sets all into a df for output.
                    name = (df3.columns[i][0],df3.columns[i][1],df3.columns[i][2],
                                df3.columns[i][3], df3.columns[j][3],u_statistic, p_value)
                    columns = list(Mann_whit_all)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Mann_whit_all = Mann_whit_all.append(data, True)
                    
    Mann_whit_all = Mann_whit_all.set_index(['Replicate 1','Plate','Condition','Batch'])
    Mann_whit_all = Mann_whit_all.sort_index(0)
    Mann_whit_all.to_csv(outputfile)
    #averages out the p-values to produce p-value mean for each replicate
    Pmean = Mann_whit_all.groupby(level=[0,1,2,3]).mean()
    Pmean = Pmean.reset_index()
    Pmean = Pmean.rename(columns={"Replicate 1":"Replicate","U-statistic":'Mean U-Stat',"P-Value":'Mean P-Value'})
    Pmean.to_csv(outputfile2,index=False)   
    return Pmean

if __name__ == "__main__":
    main()