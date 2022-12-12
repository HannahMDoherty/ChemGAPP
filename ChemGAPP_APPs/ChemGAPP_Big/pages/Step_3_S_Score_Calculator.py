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
#st.markdown("# Main page ")
#st.sidebar.markdown("# Main page")
st.title('ChemGAPP Big: S-Score Calculator')

def s_scores_calc(ipt):
    nm = ipt
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
    #print("ncont calculated")
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
        print(varcontdf, len(varcontdf))
        ar2 = np.array(varcontdf)
        varcontlist = np.nanmedian(varcontdf,axis=1)
        print(varcontlist)
        varcontfull.append(varcontlist)
    #print("vcont calculated")
    ucontfull = []
    for p in sorted(plates2):
        ucontlist = []
        df1 = (nm.xs((p), axis =1, drop_level=False))
        ar1 = np.array(df1)
        uvals = np.nanmedian(ar1,axis=1)
        ucontlist = np.append(ucontlist,uvals)
        ucontfull.append(ucontlist)
    #print("ucont calculated")
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
    return ardf
def plate_info(file):
            p = pd.read_table(file)
            if len(p.columns) == 5:
                p.columns = ['row','column','strain','info','normalisation_group']
                p = p.set_index(['row','column'])
                p = p.drop(['info','normalisation_group'], axis=1)
                p = p.sort_index()
            if len(p.columns) == 3:
                p.columns = ['row','column','strain']
                p = p.set_index(['row','column'])
                p = p.sort_index()
            return p
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

complete = st.sidebar.button(label="Begin!")
if complete:
    if genre == 'Original':
        my_expander2 = st.expander(label="Initial Dataset:")
        element4 = my_expander2.write(st.session_state.normalised_dataset)
        elementss = st.info("Calculating S-Scores...")
        nm = st.session_state.normalised_dataset
        ardf = s_scores_calc(nm)
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
        df_with_strains = pd.DataFrame(columns = sorted(columns1))
        for a, n in zip(plate_DFs, sorted(plates)):
            df1 = (nm3.xs((n), axis =1, drop_level=True))
            df2 = pd.merge(df1,a,left_index=True, right_index=True)
            df_with_strains = pd.concat([df_with_strains , df2], ignore_index=True)
        df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
        df_with_strains = df_with_strains.set_index('Gene')
        df_with_strains = df_with_strains.astype(float)
        df_with_strains.to_csv(st.session_state.outputfile+"_Final_Dataset.txt", sep="\t")
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
        ardf = s_scores_calc(nm)
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
        df_with_strains = pd.DataFrame(columns = sorted(columns1))
        for a, n in zip(plate_DFs, sorted(plates)):
            df1 = (nm3.xs((n), axis =1, drop_level=True))
            df2 = pd.merge(df1,a,left_index=True, right_index=True)
            df_with_strains = pd.concat([df_with_strains , df2], ignore_index=True)
        df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
        df_with_strains = df_with_strains.set_index('Gene')
        df_with_strains = df_with_strains.astype(float)
        df_with_strains.to_csv(st.session_state.outputfile+"_Final_Curated_Dataset.txt", sep="\t")
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
        ardf = s_scores_calc(nm)
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
        df_with_strains = pd.DataFrame(columns = sorted(columns1))
        for a, n in zip(plate_DFs, sorted(plates)):
            df1 = (nm3.xs((n), axis =1, drop_level=True))
            df2 = pd.merge(df1,a,left_index=True, right_index=True)
            df_with_strains = pd.concat([df_with_strains , df2], ignore_index=True)
        df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
        df_with_strains = df_with_strains.set_index('Gene')
        df_with_strains = df_with_strains.astype(float)
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
        ardf = s_scores_calc(nm)
        ardf.to_csv(st.session_state.outputfile+"_Curated_S_scores.csv")
        #plate_DFs = []
        #for i in plate_info_files:
        #    p = plate_info(i)
        #    plate_DFs.append(p)
        nm3 = pd.read_csv(st.session_state.outputfile+"_Curated_S_scores.csv",index_col=[0, 1], header=[0, 1])
        plates = {x[0] for x in nm3.columns}
        nm4 = nm3.copy(deep=False)
        nm4['1','strain'] = 'p'
        columns1 = {x[1] for x in nm4.columns}
        df_with_strains = pd.DataFrame(columns = sorted(columns1))
        for a, n in zip(plate_DFs, sorted(plates)):
            df1 = (nm3.xs((n), axis =1, drop_level=True))
            df2 = pd.merge(df1,a,left_index=True, right_index=True)
            df_with_strains = pd.concat([df_with_strains , df2], ignore_index=True)
        df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
        df_with_strains = df_with_strains.set_index('Gene')
        df_with_strains = df_with_strains.astype(float)
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