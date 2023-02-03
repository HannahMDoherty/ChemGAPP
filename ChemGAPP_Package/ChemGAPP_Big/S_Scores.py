#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import pandas as pd
import numpy as np

def get_options():
    parser = argparse.ArgumentParser(description="Computes the S-scores from the normalised dataset. ",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--InputFile", help="The normalised csv file from Check_Normalisation.py")
    parser.add_argument("-o", "--OutputFile", help="A CSV file of the dataset as S-scores")
    parser.add_argument("-s", "--Scale", help="Scale scores to and inter-quartile range of 1.35. Options: True, False", default=True)
    return parser.parse_args()
    
    
    
def main():
    options = get_options()
    inputfile = options.InputFile
    outputfile = options.OutputFile
    scale = options.Scale
    nm = pd.read_csv(inputfile,
                      index_col=[0, 1],
                      header=[0, 1, 2, 3])
    nm = nm.reindex(sorted(nm.columns), axis=1)
    #replaces inf values with nan
    nm = nm.replace(np.inf, np.nan)
    nm_array = np.array(nm)
    nmlen = len(nm)

    plates = {x[0:2] for x in nm.columns}
    repnumber = np.array([])
    #splits by condition and source plate and counts number of replicates then finds median number of replicates
    for p in sorted(plates):
        df1 = (nm.xs((p), axis =1, drop_level=False))
        ar1 = np.array(df1)
        length = ar1.shape[1]
        repnumber = np.append(repnumber,length)
    ncont = np.nanmedian(repnumber) 
    print("ncont calculated")  
    

    plates2 = {x[0] for x in nm.columns}
    varcontfull = []
    #splits dataset by source plate and iterates through
    for p in sorted(plates2):
        varcontdf = np.zeros((nmlen, 1))
        varcontdf.shape = (nmlen,1)
        #produces df with all conditions for a specific source plate.
        df1 = (nm.xs((p), axis =1, drop_level=False))
        #swaps levels so batch included in splitting for attributes
        df1 = df1.swaplevel(0,1,axis=1)
        df1 = df1.swaplevel(1,3,axis=1)
        #splits by condition and batch
        conds = {x[0:2] for x in df1.columns}
        for p2 in sorted(conds):
            df2 = (df1.xs((p2), axis =1, drop_level=False))
            ar1 = np.array(df2)
            #calculates variance for each row as an array
            varvallist = np.nanvar(ar1,axis=1)
            varvallist = np.array(varvallist)
            varvallist.shape = (nmlen,1)
            #produces df with variance for every mutant in every condition plate
            varcontdf = np.concatenate((varcontdf, varvallist),axis=1)
        varcontdf = pd.DataFrame(varcontdf)
        varcontdf = varcontdf.iloc[: , 1:] 
        ar2 = np.array(varcontdf)
        #calculates the mean variance for each row to make a list of vcont values.
        varcontlist = np.nanmedian(varcontdf,axis=1)
        varcontfull.append(varcontlist)
    print("vcont calculated")
        
    ucontfull = []
    #splits by plate so only same mutants in each row
    for p in sorted(plates2):
        ucontlist = []
        #calculates the median colony size for all rows across all conditions to make list of ucont values
        df1 = (nm.xs((p), axis =1, drop_level=False))
        ar1 = np.array(df1)
        uvals = np.nanmedian(ar1,axis=1)
        ucontlist = np.append(ucontlist,uvals)
        ucontfull.append(ucontlist)
    print("ucont calculated")

    
    
    ################################################
    
    
    
    #calculates median relative error for all mutants in all conditions
    medrelerror = (np.nanstd(nm_array,dtype=np.float64)/(np.nanmean(nm_array,dtype=np.float64)))
    rounds = 0
    final_s_score = pd.DataFrame(index=nm.index)
    #splits by condition and source plate
    for p in sorted(plates):
        #produce empty array for S-scores
        s_score_array = np.array([])
        rounds = rounds + 1
        print(rounds, p)
        df1 = (nm.xs((p), axis =1, drop_level=False))
        ar1 = np.array(df1)
        #produce minimum for varexp, mean standard deviation squared.
        varmin = np.nanmean(np.nanstd(ar1,axis=1))
        varmin2 = varmin*varmin
        # runs through the rows of the array zipped to the ucont and varcont list values associated to that row.
        for r,k,vc in zip(range(len(ar1)), ucontfull[int(p[0])-1],varcontfull[int(p[0])-1]):
            #calculate var exp by getting the std squared of the replicate colonysizes
            varexpstd = np.nanstd(ar1[r])
            varexp = varexpstd*varexpstd
            #if varexp lower than the minimium bound then varexp becomes the minimum bound
            if varexp < varmin2:
                varexp = varmin2
            #nexp is the number of replicates for the plate and condition. thus takes the shape value as this matches this.
            nexp = (ar1.shape[1])
            ncont = np.nanmedian(repnumber)
            # calculates uexp by finding average of the replciates colony sizes for mutant of interest
            uexp = np.nanmean(ar1[r])
            ucont = k
            #minimum bound for varcont is ucont multiplied by the median relative error
            minboundvarcont = ucont*medrelerror 
            #thus if vc larger then varcont is vc else the minimum bound is taken
            if vc > minboundvarcont:
                varcont = vc
            else:
                varcont = minboundvarcont
            #calculates svar and s-scores using the variables for each mutant of interest in the condition of interest
            svar = (varexp*(nexp-1)+varcont*(ncont-1))/(nexp+ncont-2)
            S_score = (uexp-ucont)/(np.sqrt(svar/nexp+svar/ncont))
            #appends s-score to an array to build up each row of each column to concatenate into a full dataset
            s_score_array = np.append(s_score_array,S_score)  
        ar_df = pd.DataFrame(s_score_array, index=nm.index)
        final_s_score = pd.concat([final_s_score, ar_df], axis=1)
    #adds columns back in for just the conditions and source plate numbers
    final_s_score.columns = (pd.MultiIndex.from_tuples(sorted(nm.columns.droplevel([2,3]).unique())))
    final_s_score = final_s_score.replace([np.inf, -np.inf], np.nan)
    if scale == True:
        mlen = len(final_s_score)
        ardf = np.zeros((mlen, 1))
        ardf.shape = (mlen,1)
        #scales everything such that the interquartile range of each column is equal to 1.35
        for c in sorted(final_s_score.columns):
            df1 = pd.DataFrame(final_s_score.xs(c,axis=1,drop_level=False))
            ar1 = np.array(df1)
            ar2 = np.array(df1)
            mean = np.nanmean(ar1)
            iqr = (np.nanpercentile(ar1,[75])-np.nanpercentile(ar1,[25]))
            for i in range(len(ar1)):
                ar2[i] = (ar1[i]*(1.35/iqr))
            ar2.shape = (mlen,1)
            ardf = np.concatenate((ardf,ar2), axis=1)
        ardf = pd.DataFrame(ardf, index=final_s_score.index)
        ardf = ardf.iloc[: , 1:]
        ardf.columns = (pd.MultiIndex.from_tuples(sorted(final_s_score.columns)))
        ardf = ardf[sorted(ardf)]
        ardf.to_csv(outputfile)
    if scale == "False":
        ardf = final_s_score
        ardf.to_csv(outputfile)
    return ardf

if __name__ == "__main__":
    main()