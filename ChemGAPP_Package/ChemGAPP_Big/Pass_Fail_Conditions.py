#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import pandas as pd
import numpy as np
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations

def get_options():
    parser = argparse.ArgumentParser(description="The output files of the Mann-Whitney condition level analysis and the condition variance analysis are inputted. The files are tested to see which conditions fail at certain thresholds of variance and Mann-Whitney p value.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-iv", "--InputFile_Variance", help="Output file from Condition_Variance.py")
    parser.add_argument("-imwc", "--InputFile_MWC", help="Output file from Mann_Whitney_Condition_Level.py")
    parser.add_argument("-ov", "--OutputFile_Variance", help="A CSV file showing the conditions and the thresholds at which they pass and fail. Here variances which are greater than the threshold tested fail.")
    parser.add_argument("-omwc", "--OutputFile_MWC", help="A CSV file showing the conditions and the thresholds at which they pass and fail. Here p values which are lower than the threshold tested fail.")
    return parser.parse_args()
    

def main():
    options = get_options()
    inptV = options.InputFile_Variance
    otputV = options.OutputFile_Variance
    inptMWC = options.InputFile_MWC
    otputMWC = options.OutputFile_MWC
    #### Calculating the Mann_Whitney P-Value threshold Pass and Fail values   
    Pmean = pd.read_csv(inptMWC)
    Pmean = Pmean.set_index(['Condition','Batch'])
    thres_ar = [x for x in sorted(Pmean['Mean Variance P-Value']) if x > 0]
    # takes thresholds at intervals across the achieved values
    thres_array = (sorted(thres_ar)[(int(len(thres_ar)*(1/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(4/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(8/24)))],
                   sorted(thres_ar)[int((len(thres_ar)*(12/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(16/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(20/24)))])
    
    mwc_p_f = pd.DataFrame(columns=['Condition','Batch',
                                   'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[0], precision=3)),
                                   'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[1], precision=3)),
                                   'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[2], precision=3)),
                                   'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[3], precision=3)),
                                   'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[4], precision=3)),
                                   'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[5], precision=3))])
    Pmean_array = np.array(Pmean)
    plt2 = {x[0:5] for x in Pmean.index}
    print(thres_array[5],thres_array[4],thres_array[3],thres_array[2],thres_array[1],thres_array[0])
    #Checks mean Variance p-value for each condition greater than thresholds
    #  and marks as F for fail or lower than threshold and marks as P for pass 
    # and creates dataset of Ps and Fs for the various thresholds
    for r,p in zip(range(len(Pmean_array)),sorted(plt2)):
        print(Pmean_array[r][1])
        if Pmean_array[r][1] > thres_array[5]:
            PF5 = 'F' 
        if Pmean_array[r][1] > thres_array[4]:
            PF4 = 'F' 
        if Pmean_array[r][1] > thres_array[3]:
            PF3 = 'F'
        if Pmean_array[r][1] > thres_array[2]:
            PF2 = 'F' 
        if Pmean_array[r][1] > thres_array[1]:
            PF1 = 'F'
        if Pmean_array[r][1] > thres_array[0]:
            PF0 = 'F' 
        if Pmean_array[r][1] <= thres_array[5]:
            PF5 = 'P'
        if Pmean_array[r][1] <= thres_array[4]:
            PF4 = 'P' 
        if Pmean_array[r][1] <= thres_array[3]:
            PF3 = 'P'
        if Pmean_array[r][1] <= thres_array[2]:
            PF2 = 'P' 
        if Pmean_array[r][1] <= thres_array[1]:
            PF1 = 'P'
        if Pmean_array[r][1] <= thres_array[0]:
            PF0 = 'P' 
        name = (p[0],p[1],PF0, PF1,PF2,PF3,PF4,PF5)
        columns = list(mwc_p_f)
        data = []
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        mwc_p_f = mwc_p_f.append(data, True)   
    mwc_p_f = mwc_p_f.sort_values(['Condition','Batch'])
    mwc_p_f = mwc_p_f.reset_index(drop=True)
    mwc_p_f.to_csv(otputMWC, index=False)
    
    ave_Var_cond = pd.read_csv(inptV)
    ave_Var_cond = ave_Var_cond.set_index(['Condition','Batch'])
    thres_ar = [x for x in sorted(ave_Var_cond['Average Variance']) if x > 0]
    # takes thresholds at intervals across the achieved values
    thres_array = (sorted(thres_ar)[(int(len(thres_ar)*(1/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(4/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(8/24)))],
                   sorted(thres_ar)[int((len(thres_ar)*(12/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(16/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(20/24)))])
    
    varc_p_f = pd.DataFrame(columns=['Condition','Batch',
                                   'Average Variance Threshold '+ str(np.format_float_scientific(thres_array[0], precision=3)),
                                   'Average Variance Threshold '+ str(np.format_float_scientific(thres_array[1], precision=3)),
                                   'Average Variance Threshold '+ str(np.format_float_scientific(thres_array[2], precision=3)),
                                   'Average Variance Threshold '+ str(np.format_float_scientific(thres_array[3], precision=3)),
                                   'Average Variance Threshold '+ str(np.format_float_scientific(thres_array[4], precision=3)),
                                   'Average Variance Threshold '+ str(np.format_float_scientific(thres_array[5], precision=3))])
    varc_array = np.array(ave_Var_cond)
    plt2 = {x[0] for x in ave_Var_cond.index}
    #Checks mean Variance for each condition greater than thresholds
    #  and marks as F for fail or lower than threshold and marks as P for pass 
    # and creates dataset of Ps and Fs for the various thresholds
    for r,p in zip(range(len(varc_array)),sorted(ave_Var_cond.index)):
        if varc_array [r] > thres_array[0]:
            PF0 = 'F' 
        if varc_array [r] > thres_array[1]:
            PF1 = 'F' 
        if varc_array [r] > thres_array[2]:
            PF2 = 'F'
        if varc_array [r] > thres_array[3]:
            PF3 = 'F' 
        if varc_array [r] > thres_array[4]:
            PF4 = 'F'
        if varc_array [r] > thres_array[5]:
            PF5 = 'F' 
        if varc_array [r] <= thres_array[0]:
            PF0 = 'P'
        if varc_array [r] <= thres_array[1]:
            PF1 = 'P' 
        if varc_array [r] <= thres_array[2]:
            PF2 = 'P'
        if varc_array [r] <= thres_array[3]:
            PF3 = 'P' 
        if varc_array [r] <= thres_array[4]:
            PF4 = 'P'
        if varc_array [r] <= thres_array[5]:
            PF5 = 'P' 
        name = (p[0],p[1], PF0, PF1,PF2,PF3,PF4,PF5)
        columns = list(varc_p_f)
        data = []
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        varc_p_f = varc_p_f.append(data, True)                              
    varc_p_f = varc_p_f.sort_values(['Condition'])
    varc_p_f = varc_p_f.reset_index(drop=True)
    varc_p_f.to_csv(otputV, index=False)

if __name__ == "__main__":
    main()
# %%
