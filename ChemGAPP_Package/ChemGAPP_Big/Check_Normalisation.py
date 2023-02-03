#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
from pickle import FALSE
import pandas as pd
import numpy as np
import scipy
import scipy.stats
from scipy.stats import wilcoxon
from scipy.stats import ranksums

def get_options():
    parser = argparse.ArgumentParser(description="Checks each plate individually to see if outer-edge normalisation is required due to plate effects.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--InputFile", help="The CSV file from Iris_to_dataset.py of the colony dataset with conditions across the top and row/column coordinates downwards")
    parser.add_argument("-o", "--OutputFile", help="CSV file of the normalised colony sizes.")
    parser.add_argument("-m", "--max_colony_size", help="Maximum colony size allowed, any colony larger than this will be set to this maximum",default=False,type=int)
    return parser.parse_args()

def main():
    options = get_options()
    inputfile = options.InputFile
    outputfile = options.OutputFile
    maximum_size = options.max_colony_size
    m = pd.read_csv(inputfile,index_col=[0, 1],header=[0, 1, 2, 3])
    m = m.apply(pd.to_numeric)
    m = m[sorted(m)]
    m_array = np.array(m)
    #calculates the size of the inputted plates
    mlen = len(m)
    if mlen == 1536:
        rlen = 32
        clen = 48
    if mlen == 384:
        rlen = 16
        clen = 24
    if mlen == 96:
        rlen = 8
        clen = 12
    n = pd.read_csv(inputfile,index_col=[0, 1],header=[0, 1, 2, 3])
    n = n.apply(pd.to_numeric)
    n = n[sorted(n)]
    m_ind = m.reset_index()
    # produces dataframe and array excluding the outer two rows and columns for each plate for plate middle mean (PMM) calculation
    mask = (m_ind['row'] != 1) & (m_ind['row'] != rlen) & (m_ind['column'] != 1) & (m_ind['column'] != clen) & (m_ind['row'] != 2) & (m_ind['row'] != (rlen-1)) & (m_ind['column'] != 2) & (m_ind['column'] != (clen-1))
    df_pmm = m_ind[mask]
    df_pmm = df_pmm.set_index(['row','column'])
    # produces dataframe and array of the outer two rows and columns for each plate
    df_outer = m_ind[-mask]
    df_outer = df_outer.set_index(['row','column'])
    pmm_array = np.array(df_pmm)
    #calculates median colony size for the entire dataset
    m_medain = np.nanmedian(m_array)
    
    #Change all false zeros to NaN in the dataset
    conditions = {x[0:2] for x in m.columns}
    n2 = np.zeros((mlen, 0))
    n2.shape = (mlen,0)
    for c in sorted(conditions):
        df1 = n.xs((c), axis =1, drop_level=False)
        ar1=np.array(df1)
        ar1 = ar1.astype(np.float64)
        for j in range(len(ar1)):
            if np.all(ar1[j]) == 0:
                if not np.any(ar1[j]) == 0:
                    for l, vs in enumerate(ar1[j]):
                        if vs == 0:
                            ar1[j][l] = "NaN"
        n2 = np.concatenate((n2,ar1),axis=1)
    n2 = pd.DataFrame(n2, index=m.index)
    n2.columns = (pd.MultiIndex.from_tuples(sorted(n.columns)))
    n2 = n2[sorted(n2)]
    n_array = np.array(n2)
    
    ardf = np.zeros((mlen, 1))
    ardf.shape = (mlen,1)
    rounds = 0
    #runs through each plate individually matching the plates for the outer and inner dataframes
    for c1,c2,i in zip(sorted(df_pmm.columns),sorted(df_outer.columns), range(len(m.columns))):
        rounds = rounds + 1
        print(rounds)
        ar1 = pmm_array[:,i]
        ar2 = n_array[:,i]
        # finds the colony sizes within the 40th and 60th percentiles for PMM calculation
        pmm_40 = np.nanpercentile(ar1, 40)
        pmm_60 = np.nanpercentile(ar1, 60)
        mask2= (ar1 >= pmm_40) & (ar1 <= pmm_60) 
        pmm_perc_values = ar1[mask2]
        # finds mean of these 40th-60th percentile values = PMM
        PMM = pmm_perc_values.mean()
        #compares the distributions of outer colonies and inner colonies to check if first step of normalisation is required
        df1 = df_pmm.xs((c1), axis =1, drop_level=False)
        arA=np.array(df1)
        arA= np.array(arA).flatten()
        df2 = df_outer.xs((c2), axis =1, drop_level=False)
        arB=np.array(df2)
        arB= np.array(arB).flatten()
        w, p = ranksums(arA, arB)
        # if siginificantly different performs the first step and second step. 
        if p < 0.05:
            print("Different Dist")
            for ind, j in zip(m.index,range(len(ar2))):
                # for each colony within the outer two edges this calculates the median of the row and column in which they are located
                if ind[0] == 1 or ind[0] == rlen or ind[1] == 1 or ind[1] == clen or ind[0] == 2 or ind[0] == (rlen-1) or ind[1] == 2 or ind[1] == (clen-1):  
                    if ind[0] == 1:
                        p_median_list_1 = [list(m_array[c:c+1,i]) for c, ind in zip(range(len(m_array)),m.index) if ind[0] == 1]
                        p_median = np.nanmedian(p_median_list_1)
                    elif ind[0] == rlen:
                        p_median_list_32 = [list(m_array[c:c+1,i]) for c, ind in zip(range(len(m_array)),m.index) if ind[0] == rlen]
                        p_median = np.nanmedian(p_median_list_32)
                    elif ind[1] == 1:
                        p_median_list_1_1 = [list(m_array[c:c+1,i]) for c, ind in zip(range(len(m_array)),m.index) if ind[1] == 1]
                        p_median = np.nanmedian(p_median_list_1_1)
                    elif ind[1] == clen:
                        p_median_list_48 = [list(m_array[c:c+1,i]) for c, ind in zip(range(len(m_array)),m.index) if ind[1] == clen]
                        p_median = np.nanmedian(p_median_list_48)
                    elif ind[0] == 2:
                        p_median_list_2 = [list(m_array[c:c+1,i]) for c, ind in zip(range(len(m_array)),m.index) if ind[0] == 2]
                        p_median = np.nanmedian(p_median_list_2)
                    elif ind[0] == (rlen-1):
                        p_median_list_31 = [list(m_array[c:c+1,i]) for c, ind in zip(range(len(m_array)),m.index) if ind[0] == (rlen-1)]
                        p_median = np.nanmedian(p_median_list_31)
                    elif ind[1] == 2:
                        p_median_list_2_2 = [list(m_array[c:c+1,i]) for c, ind in zip(range(len(m_array)),m.index) if ind[1] == 2]
                        p_median = np.nanmedian(p_median_list_2_2)
                    elif ind[1] == (clen-1):
                        p_median_list_47 = [list(m_array[c:c+1,i]) for c, ind in zip(range(len(m_array)),m.index) if ind[1] == (clen-1)]
                        p_median = np.nanmedian(p_median_list_47)
                    # colonys within the outer two edges are then scaled such that the median of the row or column are equal to the PMM
                    # It then scales such that the the PMM is equal the median colony size of the entire datatset.
                    ar2[j] = (((ar2[j])*(PMM/p_median))*(m_medain/PMM))

                    if maximum_size != False:
                        if ar2[j] > maximum_size:
                            ar2[j] = maximum_size
                # If colonies are within the centre of the plate this scales such that the the PMM is equal
                # the median colony size of the entire datatset.
                elif ind[0] != 1 or ind[0] != rlen or ind[1] != 1 or ind[1] != clen or ind[0] != 2 or ind[0] != (rlen-1) or ind[1] != 2 or ind[1] != (clen-1):
                    ar2[j] = ((ar2[j])*(m_medain/PMM))
                    if maximum_size != False:
                        if ar2[j] > maximum_size:
                            ar2[j] = maximum_size
            ar2.shape = (mlen,1)
            ardf = np.concatenate((ardf,ar2), axis=1)
        # if outer edge and inner colonies are not signficantly different 
        # then just scales entire plate such that the PMM is equal to the median colony size of entire dataset
        else:
            print(c1,"Same Dist")
            for ind, j in zip(m.index,range(len(ar2))):
                ar2[j] = ((ar2[j])*(m_medain/PMM))
                if maximum_size != False:
                    if ar2[j] > maximum_size:
                        ar2[j] = maximum_size
            ar2.shape = (mlen,1)
            ardf = np.concatenate((ardf,ar2), axis=1)
    ardf = pd.DataFrame(ardf, index=m.index)
    ardf = ardf.iloc[: , 1:]
    ardf.columns = (pd.MultiIndex.from_tuples(sorted(m.columns)))
    ardf = ardf[sorted(ardf)] 
    ardf.to_csv(outputfile)
    return ardf
    
if __name__ == "__main__":
    main()