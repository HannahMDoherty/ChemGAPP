#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
parser = argparse.ArgumentParser(description="Outputs a list of plates which were removed at a certain chosen threshold for the Z-score test. Also outputs a new dataset to go back into the process of normalisation and scoring, but with detrimental plates removed.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--InputFile", help="output from Z_score_count.py")
parser.add_argument("-o", "--OutputFile", help="A CSV file with the name of the plates that were removed and their file names.")
parser.add_argument("-od", "--Original_Dataset", help="The original .csv dataset used in the first stage or the output of MW_plates_to_remove.py to remove more plates")
parser.add_argument("-or", "--Output_removed", help="A .csv dataset with detrimental plates removed.")
parser.add_argument("-t", "--Threshold",type=float, help="A chosen threshold, usually based off of the bar chart produced by Bar_plot_Plate.py.")
args = vars(parser.parse_args())
input_z1 = args["InputFile"]
outpt1 = args["OutputFile"]
input_original1 = args["Original_Dataset"]
outpt_removed1 = args["Output_removed"]
threshold1 = args["Threshold"]

def Z_plates_to_remove(input_z,threshold,outpt,input_original,outpt_removed):
    import pandas as pd
    import numpy as np
    input_DF_2 = pd.read_csv(input_z)
    m = pd.read_csv(input_original,index_col=[0, 1],header=[0, 1, 2, 3])
    Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                    'Normality Percentage','File Name'])
    lst = []
    for row in input_DF_2.iterrows():
        if row[1][8] < threshold:
            abc=(str(row[1][0]),row[1][1],row[1][2],row[1][3])
            lst.append(abc)
            name = (row[1][0],row[1][1],row[1][3],row[1][2],
                    row[1][8],(row[1][1].replace(".",",").replace(" ","-")+'-'+str(row[1][0])+"-"+str(row[1][3].replace("Batch",""))+'_'+row[1][2]+'.JPG.iris'))
            columns = list(Plts_rm)
            data = []
            zipped = zip(columns, name)
            a_dictionary = dict(zipped)
            data.append(a_dictionary)
            Plts_rm = Plts_rm.append(data, True)
    Plts_rm.to_csv(outpt)
    lst2 =  m.columns.tolist()
    lst3 = []
    for i in lst:
        if i in lst2:
            lst3.append(i)
    n = m.drop(columns=lst3, axis=1)
    n.to_csv(outpt_removed)
    return Plts_rm
Z_plates_to_remove(input_z1,threshold1,outpt1,input_original1,outpt_removed1)
