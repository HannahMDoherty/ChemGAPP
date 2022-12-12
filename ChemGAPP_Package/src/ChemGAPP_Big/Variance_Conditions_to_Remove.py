#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
parser = argparse.ArgumentParser(description="Outputs a list of conditions which were removed at a certain chosen threshold for the variance test. Also outputs a new dataset to go back into the process of normalisation and scoring, but with detrimental plates removed.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--InputFile", help="output from Condition_Variance.py")
parser.add_argument("-o", "--OutputFile", help="A CSV file with the name of the plates that were removed and their file names.")
parser.add_argument("-od", "--Original_Dataset", help="The original .csv dataset used in the first stage or the output of MW_plates_to_remove.py or Z_plates_to_remove.py to remove more plates")
parser.add_argument("-or", "--Output_removed", help="A .csv dataset with detrimental plates removed.")
parser.add_argument("-t", "--Threshold",type=int, help="A chosen threshold, usually based off of the bar chart produced by Bar_plot_Condition.py.")
args = vars(parser.parse_args())
input_V1 = args["InputFile"]
outpt1 = args["OutputFile"]
input_original1 = args["Original_Dataset"]
outpt_removed1 = args["Output_removed"]
threshold1 = args["Threshold"]

def Variance_conditions_to_remove(input_V,threshold,outpt,input_original,outpt_removed):
    import pandas as pd
    import numpy as np
    input_DF_2 = pd.read_csv(input_V)
    m = pd.read_csv(input_original,index_col=[0, 1],header=[0, 1, 2, 3])
    m.columns = m.columns.swaplevel(0,1)
    m.columns = m.columns.swaplevel(1, 3)
    Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                    'Average Variance'])
    lst = []
    for row in input_DF_2.iterrows():
        if row[1][2] > threshold:
            abc=(str(row[1][0]),row[1][1])
            lst.append(abc)
            name = (row[1][0],row[1][1],row[1][2])
            columns = list(Plts_rm)
            data = []
            zipped = zip(columns, name)
            a_dictionary = dict(zipped)
            data.append(a_dictionary)
            Plts_rm = Plts_rm.append(data, True)
    Plts_rm.to_csv(outpt)
    lst2 =  {x[0:2] for x in m.columns}
    lst3 = []
    for i in lst:
        if i in lst2:
            lst3.append(i)
    n = m.drop(columns=lst3, axis=1)
    n.columns = n.columns.swaplevel(3,1)
    n.columns = n.columns.swaplevel(1,0)
    n.to_csv(outpt_removed)
    return n
    
Variance_conditions_to_remove(input_V1,threshold1,outpt1,input_original1,outpt_removed1)
