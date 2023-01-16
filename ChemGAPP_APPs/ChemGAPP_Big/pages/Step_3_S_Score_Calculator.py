import streamlit as st
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import wilcoxon
from scipy.stats import ranksums
import os
import re
import seaborn as sns
import matplotlib.pyplot as plt
st.set_page_config(layout="wide")
st.title('ChemGAPP Big: S-Score Calculator')

def s_scores_calc(ipt,scale):
    nm = ipt
    nm = nm.reindex(sorted(nm.columns), axis=1)
    #replaces inf values with nan
    nm = nm.replace(np.inf, np.nan)
    nm_array = np.array(nm)
    nmlen = len(nm)
    plates = {x[0:2] for x in nm.columns}
    repnumber = np.array([])
    #splits by condition and source plate and counts number 
    # of replicates then finds median number of replicates
    for p in sorted(plates):
        df1 = (nm.xs((p), axis =1, drop_level=False))
        ar1 = np.array(df1)
        length = ar1.shape[1]
        repnumber = np.append(repnumber,length)
    ncont = np.nanmedian(repnumber)
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
        conds = {x[0:2] for x in df1.columns}
        #splits by condition and batch
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
        print(varcontdf, len(varcontdf))
        ar2 = np.array(varcontdf)
        #calculates the mean variance for each row to make a list of vcont values.
        varcontlist = np.nanmedian(varcontdf,axis=1)
        varcontfull.append(varcontlist)

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
    if scale == "Yes":
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
    if scale == "No":
        ardf = final_s_score
    return ardf
#renames the plate info file columns and adds row and column to index before sorting by index.
def plate_info(file):
            p = pd.read_table(file)
            if len(p.columns) == 3:
                p.columns = ['row','column','strain']
                p = p.set_index(['row','column'])
                p = p.sort_index()
            else:
                print("Plate information file" + str(file) + "not in correct format.")
            return p
# looks for digits within the file names so that order of 
# plate info files is sorted plate1, plate2, plate3 etc. not plate1, plate10, plate2 etc.
def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text.name) ]

my_expandera = st.expander(label="Instructions:", expanded=True)
my_expandera.markdown(""" 
1- Upload plate information files.

-----

2- Select if you want to score the original dataset, the curated dataset or both.
""")

plate_info_files = st.sidebar.file_uploader("Upload Plate Information Files", accept_multiple_files=True)
genre = st.sidebar.radio(
     "Which dataset would you like to score?",
     ('None','Original', 'Curated', 'Both'))
scaling = st.sidebar.radio(
     "Would you like to scale the inter-quartile range to 1.35?",
     ('Yes','No'))
complete = st.sidebar.button(label="Begin!")
if complete:
    if genre == 'Original':
        my_expander2 = st.expander(label="Initial Dataset:")
        element4 = my_expander2.write(st.session_state.normalised_dataset)
        elementss = st.info("Calculating S-Scores...")
        nm = st.session_state.normalised_dataset
        ardf = s_scores_calc(nm,scaling)
        ardf.to_csv(st.session_state.outputfile+"_S_scores.csv")
        plate_info_files.sort(key=natural_keys)
        plate_DFs = []
        for i in plate_info_files:
            p = plate_info(i)
            plate_DFs.append(p)
        nm3 = pd.read_csv(st.session_state.outputfile+"_S_scores.csv",index_col=[0, 1], header=[0, 1])
        plates = {x[0] for x in nm3.columns}
        nm4 = nm3.copy(deep=False)
        nm4['1','strain'] = 'p'
        columns1 = {x[1] for x in nm4.columns}
        #names columns with condition
        df_with_strains = pd.DataFrame(columns = sorted(columns1))
        # splits by each source plate and then assigns gene names from matching plate information file.
        for a, n in zip(plate_DFs, sorted(plates)):
            df1 = (nm3.xs((n), axis =1, drop_level=True))
            df2 = pd.merge(df1,a,left_index=True, right_index=True)
            df_with_strains = pd.concat([df_with_strains , df2], ignore_index=True)
        df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
        df_with_strains = df_with_strains.set_index('Gene')
        df_with_strains = df_with_strains.astype(float)
        df_with_strains.to_csv(st.session_state.outputfile+"_Final_Dataset.txt", sep="\t")
        #averages scores for rows with the same gene name.
        groupmean = df_with_strains.groupby(level=0).mean()
        df_averaged = pd.DataFrame(groupmean, index=groupmean.index, columns=groupmean.columns)
        df_averaged.index.name = None
        df_averaged.to_csv(st.session_state.outputfile+"_Final_Dataset_Averaged.txt", sep="\t")

        if "final_DF" not in st.session_state:
            st.session_state.final_DF = df_with_strains
        elif "final_DF" in st.session_state:
            st.session_state.final_DF = df_with_strains
        
        if "final_ave_DF" not in st.session_state:
            st.session_state.final_ave_DF = df_averaged
        elif "final_ave_DF" in st.session_state:
            st.session_state.final_ave_DF = df_averaged

        my_expander3 = st.expander(label="S-Scores:")    
        element5 = my_expander3.write(df_averaged)
        elementss.empty()

    if genre == 'Curated':
        my_expander2 = st.expander(label="Curated Dataset:")
        element4 = my_expander2.write(st.session_state.normalised_curated_dataset)
        elementss = st.info("Calculating S-Scores...")
        nm = st.session_state.normalised_curated_dataset
        ardf = s_scores_calc(nm,scaling)
        ardf.to_csv(st.session_state.outputfile+"_Curated_S_scores.csv")
        plate_info_files.sort(key=natural_keys)
        plate_DFs = []
        for i in plate_info_files:
            p = plate_info(i)
            plate_DFs.append(p)
        nm3 = pd.read_csv(st.session_state.outputfile+"_Curated_S_scores.csv",index_col=[0, 1], header=[0, 1])
        plates = {x[0] for x in nm3.columns}
        nm4 = nm3.copy(deep=False)
        nm4['1','strain'] = 'p'
        columns1 = {x[1] for x in nm4.columns}
        #names columns with condition
        df_with_strains = pd.DataFrame(columns = sorted(columns1))
        # splits by each source plate and then assigns gene names from matching plate information file.
        for a, n in zip(plate_DFs, sorted(plates)):
            df1 = (nm3.xs((n), axis =1, drop_level=True))
            df2 = pd.merge(df1,a,left_index=True, right_index=True)
            df_with_strains = pd.concat([df_with_strains , df2], ignore_index=True)
        df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
        df_with_strains = df_with_strains.set_index('Gene')
        df_with_strains = df_with_strains.astype(float)
        df_with_strains.to_csv(st.session_state.outputfile+"_Final_Curated_Dataset.txt", sep="\t")
        #averages scores for rows with the same gene name.
        groupmean = df_with_strains.groupby(level=0).mean()
        df_averaged = pd.DataFrame(groupmean, index=groupmean.index, columns=groupmean.columns)
        df_averaged.index.name = None
        df_averaged.to_csv(st.session_state.outputfile+"_Final_Curated_Dataset_Averaged.txt", sep="\t")
        my_expander3 = st.expander(label="Curated S-Scores:")    
        element5 = my_expander3.write(df_averaged)
        elementss.empty()

        if "curated_DF" not in st.session_state:
            st.session_state.curated_DF = df_with_strains
        elif "curated_DF" in st.session_state:
            st.session_state.curated_DF = df_with_strains
        
        if "curated_ave_DF" not in st.session_state:
            st.session_state.curated_ave_DF = df_averaged
        elif "curated_ave_DF" in st.session_state:
            st.session_state.curated_ave_DF = df_averaged

    if genre == 'Both':
        my_expander2 = st.expander(label="Initial Dataset:")
        element4 = my_expander2.write(st.session_state.normalised_dataset)
        my_expander3 = st.expander(label="Curated Dataset:")    
        element5 = my_expander3.write(st.session_state.normalised_curated_dataset)
        elementss = st.info("Calculating S-Scores...")
        nm = st.session_state.normalised_dataset
        ardf = s_scores_calc(nm,scaling)
        ardf.to_csv(st.session_state.outputfile+"_S_scores.csv")
    
        plate_info_files.sort(key=natural_keys)
        plate_DFs = []
        for i in plate_info_files:
            p = plate_info(i)
            plate_DFs.append(p)
        nm3 = pd.read_csv(st.session_state.outputfile+"_S_scores.csv",index_col=[0, 1], header=[0, 1])
        plates = {x[0] for x in nm3.columns}
        nm4 = nm3.copy(deep=False)
        nm4['1','strain'] = 'p'
        columns1 = {x[1] for x in nm4.columns}
        #names columns with condition
        df_with_strains = pd.DataFrame(columns = sorted(columns1))
        # splits by each source plate and then assigns gene names from matching plate information file.
        for a, n in zip(plate_DFs, sorted(plates)):
            df1 = (nm3.xs((n), axis =1, drop_level=True))
            df2 = pd.merge(df1,a,left_index=True, right_index=True)
            df_with_strains = pd.concat([df_with_strains , df2], ignore_index=True)
        df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
        df_with_strains = df_with_strains.set_index('Gene')
        df_with_strains = df_with_strains.astype(float)
        #averages scores for rows with the same gene name.
        groupmean = df_with_strains.groupby(level=0).mean()
        df_averaged = pd.DataFrame(groupmean, index=groupmean.index, columns=groupmean.columns)
        df_averaged.index.name = None
        df_with_strains.to_csv(st.session_state.outputfile+"_Final_Dataset.txt", sep="\t")
        df_averaged.to_csv(st.session_state.outputfile+"_Final_Dataset_Averaged.txt", sep="\t")   
        
        if "final_DF" not in st.session_state:
            st.session_state.final_DF = df_with_strains
        elif "final_DF" in st.session_state:
            st.session_state.final_DF = df_with_strains
        
        if "final_ave_DF" not in st.session_state:
            st.session_state.final_ave_DF = df_averaged
        elif "final_ave_DF" in st.session_state:
            st.session_state.final_ave_DF = df_averaged
            
        my_expander3 = st.expander(label="S-Scores:")    
        element5 = my_expander3.write(df_averaged)
        elementss.empty()
        elementss = st.info("Calculating S-Scores...")
        
        nm = st.session_state.normalised_curated_dataset
        ardf = s_scores_calc(nm,scaling)
        ardf.to_csv(st.session_state.outputfile+"_Curated_S_scores.csv")
        nm3 = pd.read_csv(st.session_state.outputfile+"_Curated_S_scores.csv",index_col=[0, 1], header=[0, 1])
        plates = {x[0] for x in nm3.columns}
        nm4 = nm3.copy(deep=False)
        nm4['1','strain'] = 'p'
        columns1 = {x[1] for x in nm4.columns}
        #names columns with condition
        df_with_strains = pd.DataFrame(columns = sorted(columns1))
        # splits by each source plate and then assigns gene names from matching plate information file.
        for a, n in zip(plate_DFs, sorted(plates)):
            df1 = (nm3.xs((n), axis =1, drop_level=True))
            df2 = pd.merge(df1,a,left_index=True, right_index=True)
            df_with_strains = pd.concat([df_with_strains , df2], ignore_index=True)
        df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
        df_with_strains = df_with_strains.set_index('Gene')
        df_with_strains = df_with_strains.astype(float)
        #averages scores for rows with the same gene name.
        groupmean = df_with_strains.groupby(level=0).mean()
        df_averaged = pd.DataFrame(groupmean, index=groupmean.index, columns=groupmean.columns)
        df_averaged.index.name = None
        df_with_strains.to_csv(st.session_state.outputfile+"_Final_Curated_Dataset.txt", sep="\t")
        df_averaged.to_csv(st.session_state.outputfile+"_Final_Curated_Dataset_Averaged.txt", sep="\t")

        if "curated_DF" not in st.session_state:
            st.session_state.curated_DF = df_with_strains
        elif "curated_DF" in st.session_state:
            st.session_state.curated_DF = df_with_strains
        
        if "curated_ave_DF" not in st.session_state:
            st.session_state.curated_ave_DF = df_averaged
        elif "curated_ave_DF" in st.session_state:
            st.session_state.curated_ave_DF = df_averaged

        my_expander4 = st.expander(label="Curated S-Scores:")    
        element6 = my_expander4.write(df_averaged)
        elementss.empty()
        st.success("Please Continue to Step 4!")