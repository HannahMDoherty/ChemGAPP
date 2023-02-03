#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import pandas as pd

def get_options():
    parser = argparse.ArgumentParser(description="Outputs a list of conditions which were removed at a certain chosen threshold for the Mann Whitney Condition Level test. Also outputs a new dataset to go back into the process of normalisation and scoring, but with detrimental plates removed.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--InputFile", help="output from Mann_Whitney_Condition_Level.py")
    parser.add_argument("-o", "--OutputFile", help="A CSV file with the name of the plates that were removed and their file names.")
    parser.add_argument("-od", "--Original_Dataset", help="The original .csv dataset used in the first stage or the output of MW_Plates_to_Remove.py or Z_Plates_to_Remove.py or Variance_Conditions_to_Remove.py to remove more plates")
    parser.add_argument("-or", "--Output_removed", help="A .csv dataset with detrimental plates removed.")
    parser.add_argument("-t", "--Threshold",type=float, help="A chosen threshold, usually based off of the bar chart produced by Bar_plot_Condition.py.")
    return parser.parse_args()
    
def main():
    options = get_options()
    input_V = options.InputFile
    outpt = options.OutputFile
    input_original = options.Original_Dataset
    outpt_removed = options.Output_removed
    threshold = options.Threshold
    #reads the output from Mann_Whitney_Condition_Level.py
    input_DF_2 = pd.read_csv(input_V)
    #reads the original dataset before normalisation
    m = pd.read_csv(input_original,index_col=[0, 1],header=[0, 1, 2, 3])
    #swaps levels so later we can remove by just condition and batch
    m.columns = m.columns.swaplevel(0,1)
    m.columns = m.columns.swaplevel(1, 3)
    Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                    'Mean Variance P-Value','Mean Variance U-Stat'])
    #goes down row by row and checks if mean variance p-value is larger than the chosen threshold,
    # if so it is appended to a row of Plts_rm dataframe
    lst = []
    for row in input_DF_2.iterrows():
        if row[1][3] > threshold:
            #appends the condition and batch name to a list
            abc=(str(row[1][0]),row[1][1])
            lst.append(abc)
            name = (row[1][0],row[1][1],row[1][3],row[1][2])
            columns = list(Plts_rm)
            data = []
            zipped = zip(columns, name)
            a_dictionary = dict(zipped)
            data.append(a_dictionary)
            Plts_rm = Plts_rm.append(data, True)
    Plts_rm.to_csv(outpt)
    #checks if condition and batch in the original 
    # dataset supplied as you can supply a dataset that has had things removed already from another test.
    lst2 =  {x[0:2] for x in m.columns}
    lst3 = []
    for i in lst:
        if i in lst2:
            lst3.append(i)
    #drops columns which match the condition and batch.
    n = m.drop(columns=lst3, axis=1)
    n.columns = n.columns.swaplevel(3,1)
    n.columns = n.columns.swaplevel(1,0)
    n.to_csv(outpt_removed)
    
if __name__ == "__main__":
    main()