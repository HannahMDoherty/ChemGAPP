#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
parser = argparse.ArgumentParser(description="Computes the S-scores from the normalised dataset. ",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--InputFile", help="The normalised csv file from Check_Normalisation.py")
parser.add_argument("-o", "--OutputFile", help="A CSV file of the dataset as S-scores")
args = vars(parser.parse_args())
inputfile1 = args["InputFile"]
outputfile1 = args["OutputFile"]

def S_Scores(inputfile,outputfile):
    import pandas as pd
    import numpy as np
    nm = pd.read_csv(inputfile,
                      index_col=[0, 1],
                      header=[0, 1, 2, 3])
    nm = nm.reindex(sorted(nm.columns), axis=1)
    nm = nm.replace(np.inf, np.nan)
    nm_array = np.array(nm)
    nmlen = len(nm)

    plates = {x[0:2] for x in nm.columns}
    repnumber = np.array([])
    for p in sorted(plates):
        df1 = (nm.xs((p), axis =1, drop_level=False))
        ar1 = np.array(df1)
        length = ar1.shape[1]
        repnumber = np.append(repnumber,length)
    ncont = np.nanmedian(repnumber) 
    print("ncont calculated")  
    

    plates2 = {x[0] for x in nm.columns}
    varcontfull = []
    for p in sorted(plates2):
        varcontdf = np.zeros((nmlen, 1))
        varcontdf.shape = (nmlen,1)
        df1 = (nm.xs((p), axis =1, drop_level=False))
        df1 = df1.swaplevel(0,1,axis=1)
        df1 = df1.swaplevel(1,3,axis=1)
        conds = {x[0:2] for x in df1.columns}
        for p2 in sorted(conds):
            df2 = (df1.xs((p2), axis =1, drop_level=False))
            ar1 = np.array(df2)
            varvallist = np.nanvar(ar1,axis=1)
            varvallist = np.array(varvallist)
            varvallist.shape = (nmlen,1)
            varcontdf = np.concatenate((varcontdf, varvallist),axis=1)
        varcontdf = pd.DataFrame(varcontdf)
        varcontdf = varcontdf.iloc[: , 1:] 
        ar2 = np.array(varcontdf)
        varcontlist = np.nanmedian(varcontdf,axis=1)
        varcontfull.append(varcontlist)
    print("vcont calculated")
        
    ucontfull = []
    for p in sorted(plates2):
        ucontlist = []
        df1 = (nm.xs((p), axis =1, drop_level=False))
        #print(df1)
        ar1 = np.array(df1)
        uvals = np.nanmedian(ar1,axis=1)
        ucontlist = np.append(ucontlist,uvals)
        ucontfull.append(ucontlist)
    print("ucont calculated")

    
    
    ################################################
    
    
    nm3 = pd.DataFrame(columns=sorted(plates), index=nm.index)
    medrelerror = (np.nanstd(nm_array,dtype=np.float64)/(np.nanmean(nm_array,dtype=np.float64)))
    rounds = 0
    final_s_score = pd.DataFrame(index=nm.index)
    for p in sorted(plates):
        s_score_array = np.array([])
        rounds = rounds + 1
        print(rounds, p)
        df1 = (nm.xs((p), axis =1, drop_level=False))
        ar1 = np.array(df1)
        varmin = np.nanmean(np.nanstd(ar1,axis=1))
        varmin2 = varmin*varmin
        for r,k,vc in zip(range(len(ar1)), ucontfull[int(p[0])-1],varcontfull[int(p[0])-1]):
            varexpstd = np.nanstd(ar1[r])
            #varexp = ((np.nanmean(ar1[r])-np.nanmin(ar1[r]))*(np.nanmax(ar1[r])-np.nanmean(ar1[r])))
            varexp = varexpstd*varexpstd
            if varexp < varmin2:
                varexp = varmin2
            nexp = (ar1.shape[1])
            ncont = np.nanmedian(repnumber)
            uexp = np.nanmean(ar1[r])
            ucont = k
            minboundvarcont = ucont*medrelerror 
            if vc > minboundvarcont:
                varcont = vc
            else:
                varcont = minboundvarcont
            svar = (varexp*(nexp-1)+varcont*(ncont-1))/(nexp+ncont-2)
            S_score = (uexp-ucont)/(np.sqrt(svar/nexp+svar/ncont))
            s_score_array = np.append(s_score_array,S_score)  
        ar_df = pd.DataFrame(s_score_array, index=nm.index)
        final_s_score = pd.concat([final_s_score, ar_df], axis=1)
    final_s_score.columns = (pd.MultiIndex.from_tuples(sorted(nm.columns.droplevel([2,3]).unique())))
    final_s_score = final_s_score.replace([np.inf, -np.inf], np.nan)
    mlen = len(final_s_score)
    ardf = np.zeros((mlen, 1))
    ardf.shape = (mlen,1)
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
    ardf
    ardf.to_csv(outputfile)
    return ardf

S_Scores(inputfile1,outputfile1)