import streamlit as st
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import wilcoxon
from scipy.stats import ranksums
import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from io import BytesIO, StringIO
from pathlib import Path


st.set_page_config(layout="wide")
st.title('ChemGAPP Big: Normalisation')
uploaded_files = st.sidebar.file_uploader("Upload multiple IRIS files", accept_multiple_files=True)
my_expandera = st.expander(label="Instructions:", expanded=True)
my_expandera.markdown(""" 
1- Upload your IRIS files.

Ensure IRIS file names are in the format:

`CONDITION-concentration-platenumber-batchnumber_replicate.JPG.iris`

E.g. `AMPICILLIN-50 mM-6-1_B.JPG.iris`

Where concentrations have decimals, use a comma instead of a period:

E.g. `AMPICILLIN-0,5 mM-6-1_B.JPG.iris`

Where a concentration is not relevant, put two dashes between condition and plate number:

E.g. `LB--1-2_A.JPG.iris`

If only one source plate and/or only one batch was produced, assign 1 for these:

E.g. `AMPICILLIN-0,5 mM-1-1_B.JPG.iris`

`platenumber` refers to the source plate number, i.e which mutants are on the experiment plate. This will match the plate information file number in later steps.


-----

2- Enter a path to the folder you would like to save the output files to. Ensure you include a prefix which will be added to the start of all output file names e.g: ~/Projects/Test

-----

3- Input the IRIS phenotype to analyse. The spelling should exactly match the IRIS file spelling and capitalisations. Default = size. 
""")

if "outputfile" not in st.session_state:
        st.session_state.outputfile = []
def new_outfile_changed():
        if st.session_state.new_outfile:
            st.session_state.outputfile = st.session_state.new_outfile
def countPyFiles(path):
    count = 0
    indir = os.path.expanduser(path)
    for x in os.listdir(indir):
        if x.endswith(".iris"):
            count = count + 1
    return count
def getAllFiles(path):
    indir = os.path.expanduser(path)
    count = 0
    for root, dirs, files in os.walk(indir):
        #if dirs.startswith('Batch'):
        for name in dirs:
            count = count + countPyFiles(os.path.join(root,name))
    return count
m = None

otptfile = st.sidebar.text_input(label="Output Filename Prefix (Including Path):", on_change=new_outfile_changed, key="new_outfile")
phen = st.sidebar.text_input(label="IRIS Phenotype to analyse:")
complete = st.sidebar.button(label="Begin!")
if complete:
    count = len(uploaded_files)
    percent_complete = 0 
    jump = 1/count
    my_bar = st.sidebar.progress(0)
    element2 = st.sidebar.info("Files are being proccessed please wait!")
    m = None
    # cycles through iris files and uses filename to produce column headers.
    for f in uploaded_files:
        g = pd.read_csv(f,
                comment='#',
                index_col=[0, 1],
                sep='\t')
        if m is None:
            m = g[phen]
            m.name = (int(f.name.split('-')[2].split('_')[0]),
                      '-'.join(f.name.split('.')[0].split('-')[:2]).replace(",",".").replace("-"," "),
                      f.name.split('.')[0].split('_')[1],"Batch"+ f.name.split('-')[3].split('_')[0])
            m = m.to_frame()
        else:
            m1 = g[phen]
            m1.name = (int(f.name.split('-')[2].split('_')[0]),
                      '-'.join(f.name.split('.')[0].split('-')[:2]).replace(",",".").replace("-"," "),
                      f.name.split('.')[0].split('_')[1],"Batch"+ f.name.split('-')[3].split('_')[0])
            m1 = m1.to_frame()
            m = m.join(m1, how='inner')
        percent_complete = percent_complete + jump
        if percent_complete > 1:
            percent_complete = round(percent_complete, 2)
        my_bar.progress(percent_complete)
    #m = m.replace({0:np.nan}) 
    element2.empty()
    element3 = st.sidebar.success("Files Processed!")
    my_expander2 = st.expander(label="Initial Dataset:")
    if "initial_dataset" not in st.session_state:
        st.session_state.initial_dataset = m
    elif "initial_dataset" in st.session_state:
        st.session_state.initial_dataset = m
    element4 = my_expander2.write(st.session_state.initial_dataset)
    m.to_csv(otptfile+".csv")
    my_bar.empty()
    #element.empty()
    #my_expander.empty()
    element3.empty()
    elementn = st.info("Normalising Data...")
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
    n = m
    n = n.apply(pd.to_numeric)
    n = n[sorted(n)]
    n_array = np.array(n)
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
    df_outer_norm = pd.DataFrame(index=m.index)
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
                # If colonies are within the centre of the plate this scales such that the the PMM is equal
                # the median colony size of the entire datatset.    
                elif ind[0] != 1 or ind[0] != rlen or ind[1] != 1 or ind[1] != clen or ind[0] != 2 or ind[0] != (rlen-1) or ind[1] != 2 or ind[1] != (clen-1):
                    ar2[j] = ((ar2[j])*(m_medain/PMM))
            ar2.shape = (mlen,1)
            ardf = np.concatenate((ardf,ar2), axis=1)
        # if outer edge and inner colonies are not signficantly different 
        # then just scales entire plate such that the PMM is equal to the median colony size of entire dataset
        else:
            for ind, j in zip(m.index,range(len(ar2))):
                ar2[j] = ((ar2[j])*(m_medain/PMM))
            ar2.shape = (mlen,1)
            ardf = np.concatenate((ardf,ar2), axis=1)
    ardf = pd.DataFrame(ardf, index=m.index)
    ardf = ardf.iloc[: , 1:]
    ardf.columns = (pd.MultiIndex.from_tuples(sorted(m.columns)))
    ardf = ardf[sorted(ardf)]  
    elementn.empty()   
    my_expander3 = st.expander(label="Normalised Dataset:")
    if "normalised_dataset" not in st.session_state:
        st.session_state.normalised_dataset = ardf
    elif "normalised_dataset" in st.session_state:
        st.session_state.normalised_dataset = ardf
    my_expander3.write(st.session_state.normalised_dataset)
    ardf.to_csv(otptfile+"_Normalised.csv")
    elementz = st.info("Z Score Analysis Running...")
    def z_score_method(df):
        #performs z-score test and between replicates of same mutant in 
        # same condition plate and takes the absolute value of the z-score value
        z = np.abs(stats.zscore(df))
        threshold = 1.4
        outlier = []
        #goes through each z-score value and position 
        for i, v in enumerate(z):
            #if nan then continue
            if np.any(v) == 'nan':
                continue
            else:
                #if the z-score value is greater than the threshold then 
                # it is an outlier and the postion is appended to 'outlier' list
                if v > threshold:
                    outlier.append(i)
                else:
                    continue
            return outlier
    conditions = {x[0:2] for x in ardf.columns}
    replicate = {x[0:4] for x in ardf.columns}
    Plate_level = pd.DataFrame(index=ardf.index)
    ABSN = pd.DataFrame(index=ardf.index)
    for c in sorted(conditions):
        df1 = ardf.xs((c), axis =1, drop_level=False)
        ar1=np.array(df1)
        #makes array that matches the colony size data, however instead of values is all N's for Normal
        ar2=np.full(np.shape(ar1), "N")
        #iterates through row of the columns array.
        for j in range(len(ar1)):     
            #iterates through the position and value of row 
            for l, vs in enumerate(ar1[j]):
                #if the value is nan, the matching position in the N array is changed to an X
                if str(vs) == "nan":
                    ar2[j][l] = "X" 
            #if all the replicates in the row are 0 then the row is skipped 
            if np.all(ar1[j]) == 0:
                continue
            else:
                #performs the z-score method for the replicates and makes list of the outliers 
                outlier = z_score_method(ar1[j])
                #if there are outliers then it will find the mean of the replicates and then look 
                # at the positions specified as outliers and see if they are larger than the mean
                #  or smaller and designated an B or S accordingly in the N's array.
                if outlier != None:
                    if len(outlier) > 0:
                        mean = np.nanmean(ar1[j])
                        if ar1[j][outlier] < mean:
                            ar2[j][outlier] ="S"
                        else:
                            ar2[j][outlier]="B"
        #the array of Ns, Bs, Ss and Xs are formatted and saved as a df.
        ar_df = pd.DataFrame(ar2, index=ardf.index)
        ABSN = pd.concat([ABSN, ar_df], axis=1)
    ABSN.columns = (pd.MultiIndex.from_tuples(sorted(ardf.columns)))
    ABSN.to_csv(otptfile+"_Z_score.csv")
    zero_count = pd.DataFrame(columns=['Plate','Condition','Replicate','Batch','Normal','Bigger','Smaller',
                                'NaN','% Normal','% Bigger','% Smaller','% NaN'])
    replicate = {x[0:4] for x in ABSN.columns}
    #iterates through every column of the dataset and counts the N, B, S 
    # and X values then calculates the percentage of the plate they represent
    #these are then appended into the empty dataset for each plate condition replicate and batch. 
    for r in sorted(replicate):
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
        elementz.empty()
    my_expander4 = st.expander("Z Score Colony Type Counts:")
    my_expander4.write(zero_count)
    zero_count.to_csv(otptfile+"_Z_score_type_counts.csv", index=False)
    if "Z_Count" not in st.session_state:
        st.session_state.Z_Count = zero_count
    elif "Z_Count" in st.session_state:
        st.session_state.Z_Count = zero_count
    elementmwp = st.info("Mann Whitney Plate Level Analysis Running...")
    ardf_swapped = st.session_state.normalised_dataset
    ardf_swapped.columns = ardf_swapped.columns.swaplevel(2,3)
    plates = {x[0:3] for x in ardf_swapped.columns}
    Mann_whit_all = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate 1',
                                      'Replicate 2','U-statistic', 'P-Value'])
    # splits and iterates by plate, condition.
    for p in sorted(plates):
        df3 = ardf_swapped.xs((p), axis =1, drop_level=False)
        ar1 = np.array(df3)
        for i in list(range(0,len(df3.columns))):
            for j in list(range(0,len(df3.columns))):
                #compares columns that arent the same in pairs
                if i != j:
                    #takes the whole set of values then removes nan values and the corresponding values in the paired set
                    #This allows the mannwhitneyu to work
                    rep1 = ar1[:,i]
                    rep2 = ar1[:,j]
                    rep1 = rep1[~np.isnan(rep1)]
                    rep2 = rep2[~np.isnan(rep2)]  
                    u_statistic, p_value = stats.mannwhitneyu(rep1,rep2)
                    #sets all into a df for output.
                    name = (df3.columns[i][0],df3.columns[i][1],df3.columns[i][2],
                                df3.columns[i][3], df3.columns[j][3],u_statistic, p_value)
                    columns = list(Mann_whit_all)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Mann_whit_all = Mann_whit_all.append(data, True)
    Mann_whit_all = Mann_whit_all.set_index(['Replicate 1','Plate','Condition','Batch'])
    Mann_whit_all = Mann_whit_all.sort_index(0)
    #averages out the p-values to produce p-value mean for each replicate
    Pmean = Mann_whit_all.groupby(level=[0,1,2,3]).mean()
    Pmean = Pmean.reset_index()
    Pmean = Pmean.rename(columns={"Replicate 1":"Replicate","U-statistic":'Mean U-Stat',"P-Value":'Mean P-Value'})
    elementmwp.empty()
    my_expander5 = st.expander("Mann Whitney Mean p-values:")
    my_expander5.write(Pmean)
    Pmean.to_csv(otptfile+"_Mann_Whitney_Plate_level.csv",index=False)
    if "MW_Plate_Results" not in st.session_state:
        st.session_state.MW_Plate_Results = Pmean
    elif "MW_Plate_Results" in st.session_state:
        st.session_state.MW_Plate_Results = Pmean   
# condition level
    elementmwc = st.info("Mann Whitney Condition Level Analysis Running....")
    Pmean = Pmean.set_index(['Condition','Plate','Batch','Replicate'])
    cond2 = {x[0:3] for x in Pmean.index}
    P_var = pd.DataFrame(columns=['Plate','Condition','Batch','Variance U-Stat','Variance P-Value'])
    #iterates through dataset splitting by plate, condition and batch number.
    for c2 in sorted(cond2):
        df1= Pmean.xs((c2), axis =0, drop_level=False)
        ar1 = np.array(df1)
        #calculates the variance between the replicates p-value means.
        pvar_val = np.nanvar(ar1[:,0])
        #calculates the variance between the replicates u-value means.
        uvar_val = np.nanvar(ar1[:,1])
        #sets information into a dataframe
        name = (df1.index[0][1],df1.index[0][0],df1.index[0][2],pvar_val,uvar_val)
        columns = list(P_var)
        data = []
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        P_var = P_var.append(data, True)
    P_var = P_var.set_index(['Condition','Batch','Plate'])
    P_var_cond = P_var.groupby(level=[0,1]).mean()
    P_var_cond = P_var_cond.reset_index()
    P_var_cond = P_var_cond.rename(columns={'Variance P-Value':'Mean Variance P-Value','Variance U-Stat':'Mean Variance U-Stat'})
    P_var_cond.to_csv(otptfile+"_MW_Condtion_Level.csv", index = False)
    elementmwc.empty()
    my_expander6 = st.expander("Mann Whitney Mean p-value Variance:")
    my_expander6.write(P_var_cond)
    if "MW_cond_Results" not in st.session_state:
        st.session_state.MW_cond_Results = P_var_cond
    elif "MW_cond_Results" in st.session_state:
        st.session_state.MW_cond_Results = P_var_cond
    elementv = st.info("Condition Level Variance Analysis Running...")
    #makes df with same index as the input normalised dataset file.
    Var_DF = pd.DataFrame(index=ardf.index)
    conditions = {x[0:3] for x in ardf.columns}

    #splits into source plate, batch and condition, then compares the variance between the replicates.
    for c in sorted(conditions):

        df1 = ardf.xs((c), axis =1, drop_level=False)
        ar1 = np.array(df1)
        ar2 = np.array([])
        for j in range(0,len(ar1)):
            #The variance of each row is calculated and added to the variance column 
            #if all values are nan then variance is set to nan.
            if np.count_nonzero(~np.isnan(ar1[j])) == 0:
                var = "nan"
            #otherwise calculates variance ingnoring nans
            else:
                var = np.nanvar(ar1[j])
            #appends variance to array of variances
            ar2 = np.append(ar2, var)
        ar_df = pd.DataFrame(ar2, index=ardf.index)
        Var_DF = pd.concat([Var_DF, ar_df], axis=1)
    #sets column names to the source plate, batch and condition.
    Var_DF.columns = (pd.MultiIndex.from_tuples(sorted(conditions)))
    ave_Var_plate = pd.DataFrame(columns=['Condition','Batch','Plate','Average Variance'])
    #calculates the average variance for each condition.
    for f in Var_DF.columns:
        name = (f[1],f[2],f[0],np.nanmean(Var_DF[f].values.astype(float)))
        columns = list(ave_Var_plate)
        data = []
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        ave_Var_plate = ave_Var_plate.append(data, True)
    ave_Var_plate = ave_Var_plate.set_index(['Condition',"Batch",'Plate'])
    cond3 = {x[0:2] for x in ave_Var_plate.index}
    ave_Var_cond = pd.DataFrame(columns=(['Condition','Batch','Average Variance']))
    #calculates mean across different source plates and produces df.
    for cd3 in sorted(cond3):
        dfVC = ave_Var_plate.xs(cd3, axis =0, drop_level=False)
        name = (cd3[0],cd3[1],dfVC['Average Variance'].mean())
        data = []
        columns = list(ave_Var_cond)
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        ave_Var_cond = ave_Var_cond.append(data, True)
    ave_Var_cond.to_csv(otptfile+"_Average_Variance.csv",index=False)
    if "Var_cond_Results" not in st.session_state:
        st.session_state.Var_cond_Results = ave_Var_cond
    elif "Var_cond_Results" in st.session_state:
        st.session_state.Var_cond_Results = ave_Var_cond
    elementv.empty()
    my_expander7 = st.expander("Condition Variance:")
    my_expander7.write(ave_Var_cond)
    ###### PASS FAIL 
    zero_count = zero_count.set_index(['Condition','Plate','Replicate','Batch'])
    abnorm_p_f = pd.DataFrame(columns=['Filename','Condition','Plate','Replicate','Batch',
                                   'Normality Threshold 20%','Normality Threshold 30%',
                                   'Normality Threshold 40%','Normality Threshold 50%'
                                   ,'Normality Threshold 60%','Normality Threshold 80%'])
    z_count_array = np.array(zero_count)
    #Checks percentage normality for each condition less than thresholds
    #  and marks as F for fail or greater than threshold and marks as P for pass 
    # and creates dataset of Ps and Fs for the various thresholds
    plt = {x[0:5] for x in zero_count.index}
    for r,p in zip(range(len(z_count_array)),sorted(plt)):
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
        if z_count_array[r][4] >= 20:
            PF20 = 'P' 
        if z_count_array[r][4] >= 30:
            PF30 = 'P'
        if z_count_array[r][4] >= 40:
            PF40 = 'P' 
        if z_count_array[r][4] >= 50:
            PF50 = 'P'
        if z_count_array[r][4] >= 60:
            PF60 = 'P' 
        if z_count_array[r][4] >= 80:
            PF80 = 'P'
        name = ((str(p[0]).replace(" ","-").replace(".",",")+'-'+str(p[1])+'-'+str(p[3]).replace("Batch","")+'_'+str(p[2])+'.JPG.iris'),p[0]
                ,p[1],p[2],p[3], PF20,PF30,PF40,PF50,PF60,PF80)
        columns = list(abnorm_p_f)
        data = []
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        abnorm_p_f = abnorm_p_f.append(data, True)
    abnorm_p_f.to_csv(otptfile+"_Z_score_Pass_Fail.csv", index=False)
#### Calculating the Mann_Whitney P-Value threshold Pass and Fail values
    Pmean = pd.read_csv(otptfile+"_Mann_Whitney_Plate_level.csv")
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
                                   'MW Threshold '+ str(np.format_float_scientific(thres_array[0], precision=3)),
                                   'MW Threshold '+ str(np.format_float_scientific(thres_array[1], precision=3)),
                                   'MW Threshold '+ str(np.format_float_scientific(thres_array[2], precision=3)),
                                   'MW Threshold '+ str(np.format_float_scientific(thres_array[3], precision=3)),
                                   'MW Threshold '+ str(np.format_float_scientific(thres_array[4], precision=3)),
                                   'MW Threshold '+ str(np.format_float_scientific(thres_array[5], precision=3))])
    Pmean_array = np.array(Pmean)
    plt2 = {x[0:5] for x in Pmean.index}
    #Checks mean p-value for each condition less than thresholds
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
        if Pmean_array[r][1] >= thres_array[5]:
            PF5 = 'P'
        if Pmean_array[r][1] >= thres_array[4]:
            PF4 = 'P' 
        if Pmean_array[r][1] >= thres_array[3]:
            PF3 = 'P'
        if Pmean_array[r][1] >= thres_array[2]:
            PF2 = 'P' 
        if Pmean_array[r][1] >= thres_array[1]:
            PF1 = 'P'
        if Pmean_array[r][1] >= thres_array[0]:
            PF0 = 'P' 
        name = ((str(p[0]).replace(" ","-").replace(".",",")+'-'+str(p[1])+'-'+str(p[3]).replace("Batch","")+'_'+str(p[2])+'.JPG.iris'),p[0]
                ,p[1],p[2],p[3], PF0, PF1,PF2,PF3,PF4,PF5)
        columns = list(mwp_p_f)
        data = []
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        mwp_p_f = mwp_p_f.append(data, True)    
    mwp_p_f = mwp_p_f.sort_values(['Condition','Plate'])
    mwp_p_f = mwp_p_f.reset_index(drop=True)
    mwp_p_f.to_csv(otptfile+"_MW_Pass_Fail.csv", index=False)
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
           'MW Threshold '+ str(np.format_float_scientific(thres_array[0], precision=3)),
           'MW Threshold '+ str(np.format_float_scientific(thres_array[1], precision=3)),
           'MW Threshold '+ str(np.format_float_scientific(thres_array[2], precision=3)),
           'MW Threshold '+ str(np.format_float_scientific(thres_array[3], precision=3)),
           'MW Threshold '+ str(np.format_float_scientific(thres_array[4], precision=3)),
           'MW Threshold '+ str(np.format_float_scientific(thres_array[5], precision=3)),
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
    p_f_merge_reps.to_csv(otptfile+"_Pass_Fail.csv", index=False)
    ### producing the bar plot for plate level tests
    Pass_Fail = p_f_merge_reps
    Pass_Fail['Plate'] = pd.to_numeric(Pass_Fail['Plate'])
    Pass_Fail = Pass_Fail.set_index(['Condition','Plate'])
    plt5 = {x[0:2] for x in Pass_Fail.index}
    bar = pd.DataFrame(columns = ['Condition','Plate','Batch','Replicates', 'Threshold','Percentage'])
    for p5 in sorted(plt5):
        dfpf14 = Pass_Fail.xs((p5), axis=0, drop_level=False)
        #splits the data by condition and source plate number so number of F's 
        # can be counted for each threshold across replicates to produce percentage
        # fails by dividing by overall count and multiplying by 100. This then allows
        # the plotting of bar charts with percentage of plates within a condition failing at certain thresholds
        for cols in dfpf14.columns:
            if cols[0] == 'M' or cols[0] == 'N':
                failed = ((len(dfpf14[dfpf14[cols] == 'F'])/dfpf14[cols].count())*100)
                name=(dfpf14.index[0][0],dfpf14.index[0][1],dfpf14['Batch'][0],dfpf14['Replicates within Condition/Plate'][0],cols,failed)
                columns = list(bar)
                data = []
                zipped = zip(columns, name)
                a_dictionary = dict(zipped)
                data.append(a_dictionary)
                bar = bar.append(data, True) 
        bar = bar.round(2)
    bar_list = bar['Percentage'].astype(float).values
    b3 = {x for x in bar_list}
    b3 = list(sorted(b3))
    bar['Percentage'] = pd.Categorical(bar['Percentage'], categories=b3,ordered=True)
    sns.set_theme(style="whitegrid")
    #stacks the percentages so you can see the how much each percentage contributes at different thresholds.
    b = sns.displot(bar, x='Threshold', hue='Percentage',
                    multiple='stack', legend = True, palette=sns.color_palette("Spectral_r",len(b3)))
    b.set_axis_labels("", "Count of Condition With % Plates Lost")
    b.set_xticklabels(rotation=90)
    st.write("Count of Conditions with Different Percentages of Plates Lost at Different Thresholds")
    col1, col2,col3 = st.columns((2,1,1))
    col1.write("Plate Level")
    if "plate_level_bar" not in st.session_state:
        st.session_state.plate_level_bar = b
    elif "plate_level_bar" in st.session_state:
        st.session_state.plate_level_bar = b
    plate_level = col1.pyplot(st.session_state.plate_level_bar)
    Pmean = pd.read_csv(otptfile+"_MW_Condtion_Level.csv")
    Pmean = Pmean.set_index(['Condition','Batch'])
    thres_ar = [x for x in sorted(Pmean['Mean Variance P-Value']) if x > 0]
    # takes thresholds at intervals across the achieved values
    thres_array = (sorted(thres_ar)[(int(len(thres_ar)*(1/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(4/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(8/24)))],
                   sorted(thres_ar)[int((len(thres_ar)*(12/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(16/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(20/24)))])
    mwc_p_f = pd.DataFrame(columns=['Condition','Batch',
                                   'Mann Whitney Threshold '+ str(np.format_float_scientific(thres_array[0], precision=3)),
                                   'Mann Whitney Threshold '+ str(np.format_float_scientific(thres_array[1], precision=3)),
                                   'Mann Whitney Threshold '+ str(np.format_float_scientific(thres_array[2], precision=3)),
                                   'Mann Whitney Threshold '+ str(np.format_float_scientific(thres_array[3], precision=3)),
                                   'Mann Whitney Threshold '+ str(np.format_float_scientific(thres_array[4], precision=3)),
                                   'Mann Whitney Threshold '+ str(np.format_float_scientific(thres_array[5], precision=3))])
    Pmean_array = np.array(Pmean)
    plt2 = {x[0:5] for x in Pmean.index}
    #Checks mean Variance p-value for each condition greater than thresholds
    #  and marks as F for fail or lower than threshold and marks as P for pass 
    # and creates dataset of Ps and Fs for the various thresholds
    for r,p in zip(range(len(Pmean_array)),sorted(plt2)):
        if Pmean_array[r][1] > thres_array[5]:
            PF5 = 'F'
        if Pmean_array[r][1] > thres_array[4]:
            PF4 = 'F'
        if Pmean_array[r][1] > thres_array[3]:
            PF3 = 'F'
        if Pmean_array[r][1] > thres_array[2]:
            PF2 = 'F'
        if Pmean_array[r][1] > thres_array[1]:
            PF1 = 'F'
        if Pmean_array[r][1] > thres_array[0]:
            PF0 = 'F'
        if Pmean_array[r][1] <= thres_array[5]:
            PF5 = 'P'
        if Pmean_array[r][1] <= thres_array[4]:
            PF4 = 'P'
        if Pmean_array[r][1] <= thres_array[3]:
            PF3 = 'P'
        if Pmean_array[r][1] <= thres_array[2]:
            PF2 = 'P'
        if Pmean_array[r][1] <= thres_array[1]:
            PF1 = 'P'
        if Pmean_array[r][1] <= thres_array[0]:
            PF0 = 'P' 
        name = (p[0],p[1],PF0, PF1,PF2,PF3,PF4,PF5)
        columns = list(mwc_p_f)
        data = []
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        mwc_p_f = mwc_p_f.append(data, True)  
    mwc_p_f = mwc_p_f.sort_values(['Condition','Batch'])
    mwc_p_f = mwc_p_f.reset_index(drop=True)
    mwc_p_f.to_csv(otptfile+"_MW_Condition_PASS_FAIL.csv", index=False)
    ave_Var_cond = pd.read_csv(otptfile+"_Average_Variance.csv")
    ave_Var_cond = ave_Var_cond.set_index(['Condition','Batch'])
    thres_ar = [x for x in sorted(ave_Var_cond['Average Variance']) if x > 0]
     # takes thresholds at intervals across the achieved values
    thres_array = (sorted(thres_ar)[(int(len(thres_ar)*(1/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(4/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(8/24)))],
                   sorted(thres_ar)[int((len(thres_ar)*(12/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(16/24)))],
                   sorted(thres_ar)[(int(len(thres_ar)*(20/24)))])
                
    varc_p_f = pd.DataFrame(columns=['Condition','Batch',
                                   'Average Variance Threshold '+ str(np.format_float_scientific(thres_array[0], precision=3)),
                                   'Average Variance Threshold '+ str(np.format_float_scientific(thres_array[1], precision=3)),
                                   'Average Variance Threshold '+ str(np.format_float_scientific(thres_array[2], precision=3)),
                                   'Average Variance Threshold '+ str(np.format_float_scientific(thres_array[3], precision=3)),
                                   'Average Variance Threshold '+ str(np.format_float_scientific(thres_array[4], precision=3)),
                                   'Average Variance Threshold '+ str(np.format_float_scientific(thres_array[5], precision=3))])
    varc_array = np.array(ave_Var_cond)
    plt2 = {x[0] for x in ave_Var_cond.index}
    #Checks mean Variance for each condition greater than thresholds
    #  and marks as F for fail or lower than threshold and marks as P for pass 
    # and creates dataset of Ps and Fs for the various thresholds
    for r,p in zip(range(len(varc_array)),sorted(ave_Var_cond.index)):
        if varc_array [r] > thres_array[0]:
            PF0 = 'F' 
        if varc_array [r] > thres_array[1]:
            PF1 = 'F' 
        if varc_array [r] > thres_array[2]:
            PF2 = 'F'
        if varc_array [r] > thres_array[3]:
            PF3 = 'F' 
        if varc_array [r] > thres_array[4]:
            PF4 = 'F'
        if varc_array [r] > thres_array[5]:
            PF5 = 'F' 
        if varc_array [r] <= thres_array[0]:
            PF0 = 'P'
        if varc_array [r] <= thres_array[1]:
            PF1 = 'P' 
        if varc_array [r] <= thres_array[2]:
            PF2 = 'P'
        if varc_array [r] <= thres_array[3]:
            PF3 = 'P' 
        if varc_array [r] <= thres_array[4]:
            PF4 = 'P'
        if varc_array [r] <= thres_array[5]:
            PF5 = 'P' 
        name = (p[0],p[1], PF0, PF1,PF2,PF3,PF4,PF5)
        columns = list(varc_p_f)
        data = []
        zipped = zip(columns, name)
        a_dictionary = dict(zipped)
        data.append(a_dictionary)
        varc_p_f = varc_p_f.append(data, True)                              
    varc_p_f = varc_p_f.sort_values(['Condition'])
    varc_p_f = varc_p_f.reset_index(drop=True)
    varc_p_f.to_csv(otptfile+"_Var_Condition_PASS_FAIL.csv", index=False)
    Pass_Fail = mwc_p_f
    Pass_Fail = Pass_Fail.set_index(['Condition','Batch'])
    bar = pd.DataFrame(columns = ['Condition','Batch', 'Threshold','Percentage'])
    #iterates over rows by the condition and batch
    for p5 in sorted(Pass_Fail.index):
        dfpf14 = Pass_Fail.xs((p5), axis=0, drop_level=False)
        dfpf14 = pd.DataFrame(dfpf14)
        #transposes the data so number of F's can be counted to produce percentage fails by dividing by overall count and multiplying by 100.
        dfpf14 = pd.DataFrame.transpose(dfpf14)
        for cols in dfpf14.columns:
            if cols[0] == 'M' or cols[0] == 'N' or cols[0] == 'A':
                failed = ((len(dfpf14[dfpf14[cols] == 'F'])/dfpf14[cols].count())*100)
                name=(dfpf14.index[0][0],dfpf14.index[0][1],cols,failed)
                columns = list(bar)
                data = []
                zipped = zip(columns, name)
                a_dictionary = dict(zipped)
                data.append(a_dictionary)
                bar = bar.append(data, True) 
        bar = bar.round(2)
    bar_list = bar['Percentage'].astype(float).values
    b3 = {x for x in bar_list}
    b3 = list(sorted(b3))
    bar['Percentage'] = pd.Categorical(bar['Percentage'], categories=b3,ordered=True)
    sns.set_theme(style="whitegrid")
    b = sns.displot(bar, x='Threshold', hue='Percentage',
                    multiple='stack', legend = False, palette=("white","orange"))
    b.set_axis_labels("", "Count of Conditions Lost")
    b.set_xticklabels(rotation=90)
    col2.write("Condition Level: Mann Whitney")
    if "mann_cond_bar" not in st.session_state:
        st.session_state.mann_cond_bar = b
    elif "mann_cond_bar" in st.session_state:
        st.session_state.mann_cond_bar = b
    mwc_level = col2.pyplot(st.session_state.mann_cond_bar)
    Pass_Fail = varc_p_f
    Pass_Fail = Pass_Fail.set_index(['Condition','Batch'])
    bar = pd.DataFrame(columns = ['Condition','Batch', 'Threshold','Percentage'])
    #iterates over rows by the condition and batch
    for p5 in sorted(Pass_Fail.index):
        dfpf14 = Pass_Fail.xs((p5), axis=0, drop_level=False)
        dfpf14 = pd.DataFrame(dfpf14)
        #transposes the data so number of F's can be counted to produce percentage fails by dividing by overall count and multiplying by 100.
        dfpf14 = pd.DataFrame.transpose(dfpf14)
        for cols in dfpf14.columns:
            if cols[0] == 'M' or cols[0] == 'N' or cols[0] == 'A':
                failed = ((len(dfpf14[dfpf14[cols] == 'F'])/dfpf14[cols].count())*100)
                name=(dfpf14.index[0][0],dfpf14.index[0][1],cols,failed)
                columns = list(bar)
                data = []
                zipped = zip(columns, name)
                a_dictionary = dict(zipped)
                data.append(a_dictionary)
                bar = bar.append(data, True) 
        bar = bar.round(2)
    bar_list = bar['Percentage'].astype(float).values
    b3 = {x for x in bar_list}
    b3 = list(sorted(b3))
    bar['Percentage'] = pd.Categorical(bar['Percentage'], categories=b3,ordered=True)
    sns.set_theme(style="whitegrid")
    b = sns.displot(bar, x='Threshold', hue='Percentage',
                    multiple='stack', legend = False, palette=("white","orange"))
    b.set_axis_labels("", "Count of Conditions Lost")
    b.set_xticklabels(rotation=90)
    col3.write("Condition Level: Variance")
    if "var_cond_bar" not in st.session_state:
        st.session_state.var_cond_bar = b
    elif "var_cond_bar" in st.session_state:
        st.session_state.var_cond_bar = b
    mwc_level = col3.pyplot(st.session_state.var_cond_bar)
    ardf.columns = ardf.columns.swaplevel(3, 2)
    st.success("Please Continue to Step 2!")