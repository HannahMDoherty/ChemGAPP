#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
parser = argparse.ArgumentParser(description="Compares the distributions of the colony sizes of replicate plates of the same condition and determines if replicate plates have the same distribution based on the p value of the Mann whitney test.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--InputFile", help="The normalised csv file from Check_Normalisation.py output")
parser.add_argument("-o", "--OutputFile", help="A CSV file with the u statistics and p values for each comparison.")
parser.add_argument("-o2", "--OutputFile_Mean", help="A CSV file with the mean u statistics and p values for each replicate")
args = vars(parser.parse_args())
inputfile1 = args["InputFile"]
outputfile1 = args["OutputFile"]
outputfile21 = args["OutputFile_Mean"]


def Mann_Whitney_Plate_Level(inputfile,outputfile,outputfile2):
    import pandas as pd
    import numpy as np
    import scipy.stats as stats
    m = pd.read_csv(inputfile,
                      index_col=[0, 1],
                      header=[0, 1, 2, 3])
    plates = {x[0:2] for x in m.columns}
    Mann_whit_all = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate 1',
                                      'Replicate 2','U-statistic', 'P-Value'])
    for p in sorted(plates):
        df3 = m.xs((p), axis =1, drop_level=False)
        ar1 = np.array(df3)
        for i in list(range(0,len(df3.columns))):
            for j in list(range(0,len(df3.columns))):
                if i != j:
                    rep1 = ar1[:,i]
                    rep2 = ar1[:,j]
                    rep1 = rep1[~np.isnan(rep1)]
                    rep2 = rep2[~np.isnan(rep2)] 
                    u_statistic, p_value = stats.mannwhitneyu(rep1,rep2)
                    print(u_statistic,p_value)
                    name = (df3.columns[i][0],df3.columns[i][1],df3.columns[i][3],
                                df3.columns[i][2], df3.columns[j][2],u_statistic, p_value)
                    columns = list(Mann_whit_all)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Mann_whit_all = Mann_whit_all.append(data, True)
    Mann_whit_all = Mann_whit_all.set_index(['Replicate 1','Plate','Condition','Batch'])
    Mann_whit_all = Mann_whit_all.sort_index(0)
    Mann_whit_all.to_csv(outputfile)
    Pmean = Mann_whit_all.groupby(level=[0,1,2,3]).mean()
    Pmean = Pmean.reset_index()
    Pmean = Pmean.rename(columns={"Replicate 1":"Replicate","U-statistic":'Mean U-Stat',"P-Value":'Mean P-Value'})
    Pmean.to_csv(outputfile2,index=False)   
    return Pmean

Mann_Whitney_Plate_Level(inputfile1,outputfile1,outputfile21)