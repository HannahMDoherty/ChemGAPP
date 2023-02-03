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
    parser = argparse.ArgumentParser(description="The output files of the Mann-Whitney plate level analysis and the Z score analysis are inputted. The files are tested to see which conditions fail at certain thresholds of Normality and Mann-Whitney p value.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-iz", "--InputFile_Z_Score", help="output file from Z_score_count.py")
    parser.add_argument("-imwp", "--InputFile_MWP", help="output file from Mann_Whitney_Plate_Level.py")
    parser.add_argument("-oz", "--OutputFile_Z_Score", help="A CSV file showing the plates and the thresholds at which they pass and fail for the Z-score test. Here normality percentages which are lower than the threshold tested fail.")
    parser.add_argument("-omwp", "--OutputFile_MWP", help="A CSV file showing the plates and the thresholds at which they pass and fail for the Mann-Whitney test. Here p values which are lower than the threshold tested fail.")
    parser.add_argument("-mo", "--Merged_Outputfile", help="A CSV file showing the plates and the thresholds at which they pass and fail for both.")
    return parser.parse_args()


def main():
    options = get_options()
    inptZ = options.InputFile_Z_Score
    otputZ = options.OutputFile_Z_Score
    inptMWP = options.InputFile_MWP
    otputMWP = options.OutputFile_MWP
    otpt_merge = options.Merged_Outputfile
    #### Calculating the Normality threshold Pass and Fail values
    zero_count = pd.read_csv(inptZ)
    zero_count = zero_count.set_index(['Condition','Plate','Replicate','Batch'])
    abnorm_p_f = pd.DataFrame(columns=['Filename','Condition','Plate','Replicate','Batch',
                                   'Normality Threshold 20%','Normality Threshold 30%',
                                   'Normality Threshold 40%','Normality Threshold 50%'
                                   ,'Normality Threshold 60%','Normality Threshold 80%'])

    z_count_array = np.array(zero_count)
    #Checks percentage normality for each condition lower than thresholds
    #  and marks as F for fail or greater than threshold and marks as P for pass 
    # and creates dataset of Ps and Fs for the various thresholds
    plt = {x[0:5] for x in zero_count.index}
    for r,p in zip(range(len(z_count_array)),sorted(plt)):
        print(z_count_array[r][4],z_count_array[r][5])
        if z_count_array[r][4] < 20:
            PF20 = 'F' 
        if z_count_array[r][4] < 30:
            PF30 = 'F'
        if z_count_array[r][4] < 40:
            PF40 = 'F' 
        if z_count_array[r][4] < 50:
            PF50 = 'F'
        if z_count_array[r][4] < 60:
            PF60 = 'F' 
        if z_count_array[r][4] < 80:
            PF80 = 'F'
        if z_count_array[r][4] > 20:
            PF20 = 'P' 
        if z_count_array[r][4] > 30:
            PF30 = 'P'
        if z_count_array[r][4] > 40:
            PF40 = 'P' 
        if z_count_array[r][4] > 50:
            PF50 = 'P'
        if z_count_array[r][4] > 60:
            PF60 = 'P' 
        if z_count_array[r][4] > 80:
            PF80 = 'P'
        name = (p[0]+'-'+str(p[1])+"-"+p[3].replace("Batch","")+'_'+p[2]+'.JPG.iris',p[0]
                ,p[1],p[2],p[3], PF20,PF30,PF40,PF50,PF60,PF80)
        columns = list(abnorm_p_f)
        data = []
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        abnorm_p_f = abnorm_p_f.append(data, True)
    abnorm_p_f.to_csv(otputZ, index=False)
#### Calculating the Mann_Whitney P-Value threshold Pass and Fail values   
    Pmean = pd.read_csv(inptMWP)
    Pmean = Pmean.set_index(['Condition','Plate','Replicate','Batch'])
    thres_ar = [x for x in sorted(Pmean['Mean P-Value']) if x > 0]
    # takes thresholds at intervals across the achieved values
    thres_array = (sorted(thres_ar)[(int(len(thres_ar)*(1/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(4/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(8/24)))],
                   sorted(thres_ar)[int((len(thres_ar)*(12/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(16/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(20/24)))])
    
    mwp_p_f = pd.DataFrame(columns=["Filename",'Condition','Plate','Replicate','Batch',
                                   'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[0], precision=3)),
                                   'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[1], precision=3)),
                                   'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[2], precision=3)),
                                   'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[3], precision=3)),
                                   'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[4], precision=3)),
                                   'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[5], precision=3))])
    Pmean_array = np.array(Pmean)
    plt2 = {x[0:5] for x in Pmean.index}
    #Checks mean p-value for each condition lower than thresholds
    #  and marks as F for fail or greater than threshold and marks as P for pass 
    # and creates dataset of Ps and Fs for the various thresholds
    for r,p in zip(range(len(Pmean_array)),sorted(plt2)):
        if Pmean_array[r][1] < thres_array[5]:
            PF5 = 'F' 
        if Pmean_array[r][1] < thres_array[4]:
            PF4 = 'F' 
        if Pmean_array[r][1] < thres_array[3]:
            PF3 = 'F'
        if Pmean_array[r][1] < thres_array[2]:
            PF2 = 'F' 
        if Pmean_array[r][1] < thres_array[1]:
            PF1 = 'F'
        if Pmean_array[r][1] < thres_array[0]:
            PF0 = 'F' 
        if Pmean_array[r][1] > thres_array[5]:
            PF5 = 'P'
        if Pmean_array[r][1] > thres_array[4]:
            PF4 = 'P' 
        if Pmean_array[r][1] > thres_array[3]:
            PF3 = 'P'
        if Pmean_array[r][1] > thres_array[2]:
            PF2 = 'P' 
        if Pmean_array[r][1] > thres_array[1]:
            PF1 = 'P'
        if Pmean_array[r][1] > thres_array[0]:
            PF0 = 'P' 
        name = ((p[0]+'-'+str(p[1])+"-"+p[3].replace("Batch","")+'_'+p[2]+'.JPG.iris'),p[0]
                ,p[1],p[2],p[3], PF0, PF1,PF2,PF3,PF4,PF5)
        columns = list(mwp_p_f)
        data = []
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        mwp_p_f = mwp_p_f.append(data, True)                              
    mwp_p_f = mwp_p_f.sort_values(['Condition','Plate'])
    mwp_p_f = mwp_p_f.reset_index(drop=True)
    mwp_p_f.to_csv(otputMWP, index=False)
    
    #Merges the two datasets together and shows percentage fails across both tests
    p_f_merge = pd.merge(mwp_p_f,abnorm_p_f)
    p_f_merge['Percentage Fails'] = 'N'
    for index, row in p_f_merge.iterrows():
        P = 0 
        F = 0
        for cols in row:
            if cols == 'P':
                P = P +1
            if cols == 'F':
                F = F +1
        p_f_merge.iloc[index]['Percentage Fails'] = (F/(F+P))*100
    p_f_merge_ind = p_f_merge.set_index(['Condition','Plate'])
    plt4 = {x[0:2] for x in p_f_merge_ind.index}
    p_f_merge_reps = pd.DataFrame(columns = ['Filename', 'Condition', 'Plate', 'Replicate', 'Batch',
           'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[0], precision=3)),
           'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[1], precision=3)),
           'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[2], precision=3)),
           'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[3], precision=3)),
           'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[4], precision=3)),
           'Mann_Whitney Threshold '+ str(np.format_float_scientific(thres_array[5], precision=3)),
           'Normality Threshold 20%', 'Normality Threshold 30%',
           'Normality Threshold 40%', 'Normality Threshold 50%',
           'Normality Threshold 60%', 'Normality Threshold 80%',
           'Percentage Fails','Replicates within Condition/Plate'])
    for p4 in sorted(plt4):
        dfpfm = p_f_merge_ind.xs((p4), axis =0, drop_level=False)
        dfpfm['Replicates within Condition/Plate'] = "N"
        replist = dfpfm['Replicate'].values
        dfpfm = dfpfm.reset_index()
        first_column = dfpfm.pop('Filename')
        dfpfm.insert(0, 'Filename', first_column)
        for index, row in dfpfm.iterrows():
            dfpfm.loc[:,('Replicates within Condition/Plate')] = str(replist)
        p_f_merge_reps = p_f_merge_reps.append(dfpfm)
    p_f_merge_reps.to_csv(otpt_merge, index=False)
    return p_f_merge_reps

if __name__ == "__main__":
    main()