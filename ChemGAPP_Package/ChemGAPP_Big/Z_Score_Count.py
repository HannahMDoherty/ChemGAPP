#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import pandas as pd
import numpy as np

def get_options():
    parser = argparse.ArgumentParser(description="Counts the number of each colony type within each plate and the percentage of each colony type.",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--InputFile", help="The CSV output file from Z_Score.py")
    parser.add_argument("-o", "--OutputFile", help="A CSV file with the counts and percentages of the each colony type.")
    return parser.parse_args()
    

def main():
    options = get_options()
    inputfile = options.InputFile
    outputfile = options.OutputFile
    #reads the csv file of the datasets of Ns,Bs,Ss and Xs from Z_score.py
    ABSN = pd.read_csv(inputfile,
                      index_col=[0, 1],
                      header=[0, 1, 2, 3])
    print("File Uploaded")
    #create empty dataframe for data to append into.
    zero_count = pd.DataFrame(columns=['Plate','Condition','Replicate','Batch','Normal','Bigger','Smaller',
                                   'NaN','% Normal','% Bigger','% Smaller','% NaN'])
    replicate = {x[0:4] for x in ABSN.columns}
    rounds = 0
    #iterates through every column of the dataset and counts the N, B, S 
    # and X values then calculates the percentage of the plate they represent
    #these are then appended into the empty dataset for each plate condition replicate and batch. 
    for r in sorted(replicate):
        rounds = rounds + 1
        print(rounds,r)
        df1 = ABSN.xs((r), axis =1, drop_level=False)
        ar1 = np.array(df1)
        count_S = int(np.count_nonzero(ar1 == "S", axis=0))
        count_B = int(np.count_nonzero(ar1 == "B", axis=0))
        count_N = int(np.count_nonzero(ar1 == "N", axis=0))
        count_nan = int(np.count_nonzero(ar1 == "X", axis=0))
        total = (count_S+count_N+count_B+count_nan)
        B_perc = (count_B/total)*100
        S_perc = (count_S/total)*100
        N_perc = (count_N/total)*100
        Nan_perc = (count_nan/total)*100
        name = (df1.columns[0][0],df1.columns[0][1],
                    df1.columns[0][2],df1.columns[0][3],
                    count_N, count_B, count_S, count_nan,N_perc,B_perc,
                    S_perc,Nan_perc)
        columns = list(zero_count)
        data = []
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        zero_count = zero_count.append(data, True)
        zero_count.to_csv(outputfile, index=False)
    return zero_count

if __name__ == "__main__":
    main()