import streamlit as st
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import wilcoxon
from scipy.stats import ranksums
import os
import seaborn as sns
import matplotlib.pyplot as plt
import mpl_scatter_density
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats.stats import pearsonr
from itertools import combinations

st.set_page_config(layout="wide")
#st.markdown("# Main page ")
#st.sidebar.markdown("# Main page")
st.title('ChemGAPP Big: Threshold Selector')

my_expanderb = st.expander(label="Instructions:", expanded=True)
my_expanderb.markdown("""
1- Type the threshold values into the corresponding boxes based on the bar plots below. 

-----

2- Then select which statistical tests you would like to use to remove detrimental data. *Any combination can be selected.* 

-----

3- Press `Begin Quality Control Tests!`

-----

4- Check the replicate reproducibilty 'r' value for improvement after curation. Aim to maximise this r value with your curation.

"""
)

if "mwpt" not in st.session_state:
        st.session_state.mwpt = []
def new_mwpt_changed():
        if st.session_state.new_mwpt:
            st.session_state.mwpt = float(st.session_state.new_mwpt)
if "mwct" not in st.session_state:
        st.session_state.mwct = []
def new_mwct_changed():
        if st.session_state.new_mwct:
            st.session_state.mwct = float(st.session_state.new_mwct)
if "vt" not in st.session_state:
        st.session_state.vt = []
def new_vt_changed():
        if st.session_state.new_vt:
            st.session_state.vt = float(st.session_state.new_vt)
if "np" not in st.session_state:
        st.session_state.np = []
def new_np_changed():
        if st.session_state.new_np:
            st.session_state.np = float(st.session_state.new_np)
if st.session_state.outputfile:
    my_expander2 = st.expander(label="Initial Dataset:")
    element4 = my_expander2.write(st.session_state.initial_dataset)
    my_expander4 = st.expander(label="Bar Plots:")
    my_expander4.write("Count of Conditions with Different Percentages of Plates Lost at Different Thresholds")
    col1, col2,col3 = my_expander4.columns((2,1,1))
    col1.write("Plate Level")
    plate_level = col1.pyplot(st.session_state.plate_level_bar)
    col2.write("Condition Level: Mann Whitney")
    mwc_level = col2.pyplot(st.session_state.mann_cond_bar)
    col3.write("Condition Level: Variance")
    mwc_level = col3.pyplot(st.session_state.var_cond_bar)
    inp = st.sidebar.text_input('Percentage Normality Threshold:',on_change=new_np_changed, key="new_np")
    imwpt = st.sidebar.text_input('Mann Whitney Plate Level Threshold:',on_change=new_mwpt_changed, key="new_mwpt")
    imwct = st.sidebar.text_input('Mann Whitney Condition Level Threshold:',on_change=new_mwct_changed, key="new_mwct")
    ivt = st.sidebar.text_input('Average Variance Threshold:',on_change=new_vt_changed, key="new_vt")

    options = st.sidebar.multiselect(
        'Please select all statistical tests would you like to use to remove data:',
        ['Percentage Normality', 'Mann Whitney Plate Level', 'Mann Whitney Condition Level', 'Average Variance'])
    complete = st.sidebar.button(label="Begin Plate Removal!")
    if complete:
        if len(options) == 0:
            st.warning("No tests selected! If you do not wish to remove data, continue to Step 3 and select 'Original'")
        if set(options) == {'Percentage Normality'}:
            info1 = st.info("Running Plate Removal...")
            input_DF_2 = st.session_state.Z_Count
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Normality Percentage','File Name'])
            lst = []
            #goes down row by row and checks if percentage normality is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][8] < st.session_state.np:
                    #appends the condition, source plate, replicate and batch name to a list
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
            Plts_rm.to_csv(st.session_state.outputfile+"_Normality_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            
            thresholdset = ("Normality Threshold = "+ str(st.session_state.np))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="Curated Dataset:")
            element4 = my_expander2.write(n)
        elif set(options) == {'Mann Whitney Plate Level'}:
            info1 = st.info("Running Plate Removal...")
            input_DF_2 = st.session_state.MW_Plate_Results
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Mean P-Value','Mean U-Stat','File Name'])
            lst = []
            #goes down row by row and checks if mean p-value is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][5] < st.session_state.mwpt:
                    #appends the condition, source plate, replicate and batch name to a list
                    abc=(str(row[1][1]),row[1][2],row[1][0],row[1][3])
                    lst.append(abc)
                    name = (row[1][1],row[1][2],row[1][3],row[1][0],
                        row[1][5],row[1][4],(row[1][2].replace(".",",").replace(" ","-")+'-'+str(row[1][1])+"-"+str(row[1][3].replace("Batch",""))+'_'+row[1][0]+'.JPG.iris'))
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_MWP_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            
            thresholdset = ("MW Plate Threshold = "+ str(st.session_state.mwpt))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="New Dataset:")
            element4 = my_expander2.write(n)
        elif set(options) == {'Average Variance'}:
            info1 = st.info("Running Plate Removal...")
            input_DF_2 = st.session_state.Var_cond_Results
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Average Variance'])
            lst = []
            #goes down row by row and checks if mean variance is larger than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][2] > st.session_state.vt:
                    #appends the condition and batch name to a list                    
                    abc=(str(row[1][0]),row[1][1])
                    lst.append(abc)
                    name = (row[1][0],row[1][1],row[1][2])
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_Var_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            thresholdset = ("Average Variance Threshold = "+ str(st.session_state.vt))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="New Dataset:")
            element4 = my_expander2.write(n)
        elif set(options) == {'Mann Whitney Condition Level'}:
            info1 = st.info("Running Plate Removal...")
            input_DF_2 = st.session_state.MW_cond_Results
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Mean Variance P-Value','Mean Variance U-Stat'])
            lst = []
            #goes down row by row and checks if mean variance p-value is larger than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][3] > st.session_state.mwct:
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
            Plts_rm.to_csv(st.session_state.outputfile+"_MWC_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            thresholdset = ("MW Condition Threshold = "+ str(st.session_state.mwct))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="New Dataset:")
            element4 = my_expander2.write(n)
        elif set(options) == {'Mann Whitney Plate Level','Percentage Normality'}:
            info1 = st.info("Running Plate Removal...")
            input_DF_2 = st.session_state.Z_Count
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Normality Percentage','File Name'])
            lst = []
            #goes down row by row and checks if percentage normality is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][8] <  st.session_state.np:
                    #appends the condition, source plate, replicate and batch name to a list
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
            Plts_rm.to_csv(st.session_state.outputfile+"_Normality_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            input_DF_2 = st.session_state.MW_Plate_Results
            m = n
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Mean P-Value','Mean U-Stat','File Name'])
            lst = []
            #goes down row by row and checks if mean p-value is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][5] < st.session_state.mwpt:
                    #appends the condition, source plate, replicate and batch name to a list
                    abc=(str(row[1][1]),row[1][2],row[1][0],row[1][3])
                    lst.append(abc)
                    name = (row[1][1],row[1][2],row[1][3],row[1][0],
                            row[1][5],row[1][4],(row[1][2].replace(".",",").replace(" ","-")+'-'+str(row[1][1])+"-"+str(row[1][3].replace("Batch",""))+'_'+row[1][0]+'.JPG.iris'))
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_MWP_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            thresholdset = ("Normality Threshold = "+ str(st.session_state.np)+", MW Plate Threshold = "+ str(st.session_state.mwpt))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="New Dataset:")
            element4 = my_expander2.write(n)
        
        elif set(options) == {'Mann Whitney Condition Level','Mann Whitney Plate Level','Percentage Normality'}:
            info1 = st.info("Running Plate Removal...")
            input_DF_2 = st.session_state.Z_Count
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Normality Percentage','File Name'])
            lst = []
            #goes down row by row and checks if percentage normality is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe 
            for row in input_DF_2.iterrows():
                if row[1][8] <  st.session_state.np:
                    #appends the condition, source plate, replicate and batch name to a list
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
            Plts_rm.to_csv(st.session_state.outputfile+"_Normality_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            input_DF_2 = st.session_state.MW_Plate_Results
            m = n
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Mean P-Value','Mean U-Stat','File Name'])
            lst = []
            #goes down row by row and checks if mean p-value is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][5] < st.session_state.mwpt:
                    #appends the condition, source plate, replicate and batch name to a list
                    abc=(str(row[1][1]),row[1][2],row[1][0],row[1][3])
                    lst.append(abc)
                    name = (row[1][1],row[1][2],row[1][3],row[1][0],
                            row[1][5],row[1][4],(row[1][2].replace(".",",").replace(" ","-")+'-'+str(row[1][1])+"-"+str(row[1][3].replace("Batch",""))+'_'+row[1][0]+'.JPG.iris'))
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_MWP_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            input_DF_2 = st.session_state.MW_cond_Results
            m = n
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Mean Variance P-Value','Mean Variance U-Stat'])
            lst = []
            #goes down row by row and checks if mean variance p-value is larger than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][3] > st.session_state.mwct:
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
            Plts_rm.to_csv(st.session_state.outputfile+"_MWC_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            thresholdset = ("Normality Threshold = "+ str(st.session_state.np)+", MW Plate Threshold = "+ str(st.session_state.mwpt)+", MW Condition Threshold = "+ str(st.session_state.mwct))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="New Dataset:")
            element4 = my_expander2.write(n)
        
        elif set(options) == {'Average Variance','Mann Whitney Condition Level','Mann Whitney Plate Level','Percentage Normality'}:
            info1 = st.info("Running Plate Removal...")
            input_DF_2 = st.session_state.Z_Count
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Normality Percentage','File Name'])
            lst = []
            #goes down row by row and checks if percentage normality is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe 
            for row in input_DF_2.iterrows():
                if row[1][8] <  st.session_state.np:
                    #appends the condition, source plate, replicate and batch name to a list
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
            Plts_rm.to_csv(st.session_state.outputfile+"_Normality_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            input_DF_2 = st.session_state.MW_Plate_Results
            m = n
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Mean P-Value','Mean U-Stat','File Name'])
            lst = []
            #goes down row by row and checks if mean p-value is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][5] < st.session_state.mwpt:
                    #appends the condition, source plate, replicate and batch name to a list
                    abc=(str(row[1][1]),row[1][2],row[1][0],row[1][3])
                    lst.append(abc)
                    name = (row[1][1],row[1][2],row[1][3],row[1][0],
                            row[1][5],row[1][4],(row[1][2].replace(".",",").replace(" ","-")+'-'+str(row[1][1])+"-"+str(row[1][3].replace("Batch",""))+'_'+row[1][0]+'.JPG.iris'))
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_MWP_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            input_DF_2 = st.session_state.MW_cond_Results
            m = n
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Mean Variance P-Value','Mean Variance U-Stat'])
            lst = []
            #goes down row by row and checks if mean variance p-value is larger than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][3] > st.session_state.mwct:
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
            Plts_rm.to_csv(st.session_state.outputfile+"_MWC_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            input_DF_2 = st.session_state.Var_cond_Results
            m = n
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Average Variance'])
            lst = []
            for row in input_DF_2.iterrows():
                if row[1][2] > st.session_state.vt:
                    #appends the condition and batch name to a list
                    abc=(str(row[1][0]),row[1][1])
                    lst.append(abc)
                    name = (row[1][0],row[1][1],row[1][2])
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_Var_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            thresholdset = ("Normality Threshold = "+ str(st.session_state.np)+", MW Plate Threshold = "+ str(st.session_state.mwpt)+", MW Condition Threshold = "+ str(st.session_state.mwct)+ ", Average Variance Threshold = "+ str(st.session_state.vt))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="New Dataset:")
            element4 = my_expander2.write(n)

        elif set(options) == {'Mann Whitney Condition Level','Percentage Normality'}:
            info1 = st.info("Running Plate Removal...")
            input_DF_2 = st.session_state.Z_Count
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Normality Percentage','File Name'])
            lst = []
            #goes down row by row and checks if percentage normality is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe 
            for row in input_DF_2.iterrows():
                if row[1][8] <  st.session_state.np:
                    #appends the condition, source plate, replicate and batch name to a list
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
            Plts_rm.to_csv(st.session_state.outputfile+"_Normality_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            input_DF_2 = st.session_state.MW_cond_Results
            m = n
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Mean Variance P-Value','Mean Variance U-Stat'])
            lst = []
            #goes down row by row and checks if mean variance p-value is larger than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][3] > st.session_state.mwct:
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
            Plts_rm.to_csv(st.session_state.outputfile+"_MWC_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            thresholdset = ("Normality Threshold = "+ str(st.session_state.np)+", MW Condition Threshold = "+ str(st.session_state.mwct))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="New Dataset:")
            element4 = my_expander2.write(n)

        elif set(options) == {'Mann Whitney Condition Level','Mann Whitney Plate Level'}:
            info1 = st.info("Running Plate Removal...")
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            input_DF_2 = st.session_state.MW_Plate_Results
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Mean P-Value','Mean U-Stat','File Name'])
            lst = []
            #goes down row by row and checks if mean p-value is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][5] < st.session_state.mwpt:
                    #appends the condition, source plate, replicate and batch name to a list
                    abc=(str(row[1][1]),row[1][2],row[1][0],row[1][3])
                    lst.append(abc)
                    name = (row[1][1],row[1][2],row[1][3],row[1][0],
                            row[1][5],row[1][4],(row[1][2].replace(".",",").replace(" ","-")+'-'+str(row[1][1])+"-"+str(row[1][3].replace("Batch",""))+'_'+row[1][0]+'.JPG.iris'))
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_MWP_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            input_DF_2 = st.session_state.MW_cond_Results
            m = n
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Mean Variance P-Value','Mean Variance U-Stat'])
            lst = []
            #goes down row by row and checks if mean variance p-value is larger than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][3] > st.session_state.mwct:
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
            Plts_rm.to_csv(st.session_state.outputfile+"_MWC_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            thresholdset = ("MW Plate Threshold = "+ str(st.session_state.mwpt)+", MW Condition Threshold = "+ str(st.session_state.mwct))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="New Dataset:")
            element4 = my_expander2.write(n)

        elif set(options) == {'Average Variance','Percentage Normality'}:
            info1 = st.info("Running Plate Removal...")
            input_DF_2 = st.session_state.Z_Count
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Normality Percentage','File Name'])
            lst = []
            #goes down row by row and checks if percentage normality is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe 
            for row in input_DF_2.iterrows():
                if row[1][8] <  st.session_state.np:
                    #appends the condition, source plate, replicate and batch name to a list
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
            Plts_rm.to_csv(st.session_state.outputfile+"_Normality_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            input_DF_2 = st.session_state.Var_cond_Results
            m = n
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Average Variance'])
            lst = []
            for row in input_DF_2.iterrows():
                if row[1][2] > st.session_state.vt:
                    #appends the condition and batch name to a list
                    abc=(str(row[1][0]),row[1][1])
                    lst.append(abc)
                    name = (row[1][0],row[1][1],row[1][2])
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_Var_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            thresholdset = ("Normality Threshold = "+ str(st.session_state.np)+ ", Average Variance Threshold = "+ str(st.session_state.vt))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="New Dataset:")
            element4 = my_expander2.write(n)
        
        elif set(options) == {'Average Variance','Mann Whitney Plate Level'}:
            info1 = st.info("Running Plate Removal...")
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            input_DF_2 = st.session_state.MW_Plate_Results
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Mean P-Value','Mean U-Stat','File Name'])
            lst = []
            #goes down row by row and checks if mean p-value is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][5] < st.session_state.mwpt:
                    #appends the condition, source plate, replicate and batch name to a list
                    abc=(str(row[1][1]),row[1][2],row[1][0],row[1][3])
                    lst.append(abc)
                    name = (row[1][1],row[1][2],row[1][3],row[1][0],
                            row[1][5],row[1][4],(row[1][2].replace(".",",").replace(" ","-")+'-'+str(row[1][1])+"-"+str(row[1][3].replace("Batch",""))+'_'+row[1][0]+'.JPG.iris'))
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_MWP_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            input_DF_2 = st.session_state.Var_cond_Results
            m = n
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Average Variance'])
            lst = []
            for row in input_DF_2.iterrows():
                if row[1][2] > st.session_state.vt:
                    #appends the condition and batch name to a list
                    abc=(str(row[1][0]),row[1][1])
                    lst.append(abc)
                    name = (row[1][0],row[1][1],row[1][2])
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_Var_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            thresholdset = ("MW Plate Threshold = "+ str(st.session_state.mwpt)+", Average Variance Threshold = "+ str(st.session_state.vt))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="New Dataset:")
            element4 = my_expander2.write(n)
        
        elif set(options) == {'Average Variance','Mann Whitney Condition Level'}:
            info1 = st.info("Running Plate Removal...")
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            input_DF_2 = st.session_state.MW_cond_Results
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Mean Variance P-Value','Mean Variance U-Stat'])
            lst = []
            #goes down row by row and checks if mean variance p-value is larger than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][3] > st.session_state.mwct:
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
            Plts_rm.to_csv(st.session_state.outputfile+"_MWC_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            input_DF_2 = st.session_state.Var_cond_Results
            m = n
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Average Variance'])
            lst = []
            for row in input_DF_2.iterrows():
                if row[1][2] > st.session_state.vt:
                    #appends the condition and batch name to a list
                    abc=(str(row[1][0]),row[1][1])
                    lst.append(abc)
                    name = (row[1][0],row[1][1],row[1][2])
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_Var_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            thresholdset = ("MW Condition Threshold = "+ str(st.session_state.mwct)+ ", Average Variance Threshold = "+ str(st.session_state.vt))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="New Dataset:")
            element4 = my_expander2.write(n)

        elif set(options) == {'Average Variance','Mann Whitney Plate Level','Percentage Normality'}:
            info1 = st.info("Running Plate Removal...")
            input_DF_2 = st.session_state.Z_Count
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Normality Percentage','File Name'])
            lst = []
            #goes down row by row and checks if percentage normality is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe 
            for row in input_DF_2.iterrows():
                if row[1][8] <  st.session_state.np:
                    #appends the condition, source plate, replicate and batch name to a list
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
            Plts_rm.to_csv(st.session_state.outputfile+"_Normality_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            input_DF_2 = st.session_state.MW_Plate_Results
            m = n
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Mean P-Value','Mean U-Stat','File Name'])
            lst = []
            #goes down row by row and checks if mean p-value is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][5] < st.session_state.mwpt:
                    #appends the condition, source plate, replicate and batch name to a list
                    abc=(str(row[1][1]),row[1][2],row[1][0],row[1][3])
                    lst.append(abc)
                    name = (row[1][1],row[1][2],row[1][3],row[1][0],
                            row[1][5],row[1][4],(row[1][2].replace(".",",").replace(" ","-")+'-'+str(row[1][1])+"-"+str(row[1][3].replace("Batch",""))+'_'+row[1][0]+'.JPG.iris'))
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_MWP_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            input_DF_2 = st.session_state.Var_cond_Results
            m = n
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Average Variance'])
            lst = []
            for row in input_DF_2.iterrows():
                if row[1][2] > st.session_state.vt:
                    #appends the condition and batch name to a list
                    abc=(str(row[1][0]),row[1][1])
                    lst.append(abc)
                    name = (row[1][0],row[1][1],row[1][2])
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_Var_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            thresholdset = ("Normality Threshold = "+ str(st.session_state.np)+", MW Plate Threshold = "+ str(st.session_state.mwpt)+ ", Average Variance Threshold = "+ str(st.session_state.vt))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="New Dataset:")
            element4 = my_expander2.write(n)

        elif set(options) == {'Average Variance','Mann Whitney Condition Level','Mann Whitney Plate Level'}:
            info1 = st.info("Running Plate Removal...")
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            input_DF_2 = st.session_state.MW_Plate_Results
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Mean P-Value','Mean U-Stat','File Name'])
            lst = []
            #goes down row by row and checks if mean p-value is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][5] < st.session_state.mwpt:
                    #appends the condition, source plate, replicate and batch name to a list
                    abc=(str(row[1][1]),row[1][2],row[1][0],row[1][3])
                    lst.append(abc)
                    name = (row[1][1],row[1][2],row[1][3],row[1][0],
                            row[1][5],row[1][4],(row[1][2].replace(".",",").replace(" ","-")+'-'+str(row[1][1])+"-"+str(row[1][3].replace("Batch",""))+'_'+row[1][0]+'.JPG.iris'))
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_MWP_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            input_DF_2 = st.session_state.MW_cond_Results
            m = n
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Mean Variance P-Value','Mean Variance U-Stat'])
            lst = []
            #goes down row by row and checks if mean variance p-value is larger than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][3] > st.session_state.mwct:
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
            Plts_rm.to_csv(st.session_state.outputfile+"_MWC_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            input_DF_2 = st.session_state.Var_cond_Results
            m = n
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Average Variance'])
            lst = []
            for row in input_DF_2.iterrows():
                if row[1][2] > st.session_state.vt:
                    #appends the condition and batch name to a list
                    abc=(str(row[1][0]),row[1][1])
                    lst.append(abc)
                    name = (row[1][0],row[1][1],row[1][2])
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_Var_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            thresholdset = ("MW Plate Threshold = "+ str(st.session_state.mwpt)+", MW Condition Threshold = "+ str(st.session_state.mwct)+ ", Average Variance Threshold = "+ str(st.session_state.vt))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="New Dataset:")
            element4 = my_expander2.write(n)
        
        elif set(options) == {'Average Variance','Mann Whitney Condition Level','Percentage Normality'}:
            info1 = st.info("Running Plate Removal...")
            input_DF_2 = st.session_state.Z_Count
            m = pd.read_csv(st.session_state.outputfile+".csv",index_col=[0, 1],header=[0, 1, 2, 3])
            Plts_rm = pd.DataFrame(columns=['Plate','Condition','Batch','Replicate',
                                            'Normality Percentage','File Name'])
            lst = []
            #goes down row by row and checks if percentage normality is less than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe 
            for row in input_DF_2.iterrows():
                if row[1][8] <  st.session_state.np:
                    #appends the condition, source plate, replicate and batch name to a list
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
            Plts_rm.to_csv(st.session_state.outputfile+"_Normality_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 =  m.columns.tolist()
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            input_DF_2 = st.session_state.MW_cond_Results
            m = n
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Mean Variance P-Value','Mean Variance U-Stat'])
            lst = []
            #goes down row by row and checks if mean variance p-value is larger than the chosen threshold,
            # if so it is appended to a row of Plts_rm dataframe
            for row in input_DF_2.iterrows():
                if row[1][3] > st.session_state.mwct:
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
            Plts_rm.to_csv(st.session_state.outputfile+"_MWC_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            input_DF_2 = st.session_state.Var_cond_Results
            m = n
            #swaps levels so later we can remove by just condition and batch
            m.columns = m.columns.swaplevel(0,1)
            m.columns = m.columns.swaplevel(1, 3)
            Plts_rm = pd.DataFrame(columns=['Condition','Batch',
                                            'Average Variance'])
            lst = []
            for row in input_DF_2.iterrows():
                if row[1][2] > st.session_state.vt:
                    #appends the condition and batch name to a list
                    abc=(str(row[1][0]),row[1][1])
                    lst.append(abc)
                    name = (row[1][0],row[1][1],row[1][2])
                    columns = list(Plts_rm)
                    data = []
                    zipped = zip(columns, name)
                    a_dictionary = dict(zipped)
                    data.append(a_dictionary)
                    Plts_rm = Plts_rm.append(data, True)
            Plts_rm.to_csv(st.session_state.outputfile+"_Var_files_removed.csv")
            #checks if condition and batch in the original 
            # dataset supplied as you can supply a dataset that has had things removed already from another test.
            lst2 = {x[0:2] for x in m.columns}
            lst3 = []
            for i in lst:
                if i in lst2:
                    lst3.append(i)
            #drops columns within the list.
            n = m.drop(columns=lst3, axis=1)
            n.columns = n.columns.swaplevel(3,1)
            n.columns = n.columns.swaplevel(1,0)
            thresholdset = ("Normality Threshold = "+ str(st.session_state.np)+", MW Condition Threshold = "+ str(st.session_state.mwct)+ ", Average Variance Threshold = "+ str(st.session_state.vt))
            n.to_csv(st.session_state.outputfile+"_curated_dataset.csv", sep="\t")
            info1.empty()
            my_expander2 = st.expander(label="New Dataset:")
            element4 = my_expander2.write(n)

        if len(options) > 0:
            if "curated_dataset" not in st.session_state:
                    st.session_state.curated_dataset = n
            elif "curated_dataset" in st.session_state:
                    st.session_state.curated_dataset = n     
            elementn = st.info("Normalising Curated Dataset...")
            m = n
            m = m[sorted(m)]
            m_array = np.array(m)
            mlen = len(m)
            #calculates the size of the inputted plates
            if mlen == 1536:
                rlen = 32
                clen = 48
            if mlen == 384:
                rlen = 16
                clen = 24
            if mlen == 96:
                rlen = 8
                clen = 12
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
            my_expander3 = st.expander(label="Normalised Curated Dataset:")
            element5 = my_expander3.write(ardf)
            if "normalised_curated_dataset" not in st.session_state:
                st.session_state.normalised_curated_dataset = ardf
            elif "normalised_curated_dataset" in st.session_state:
                st.session_state.normalised_curated_dataset = ardf
            ardf.to_csv(st.session_state.outputfile+"_curated_Normalised.csv")
            if "Threshold_vals" not in st.session_state:
                st.session_state.Threshold_vals = thresholdset
            elif "Threshold_vals" in st.session_state:
                st.session_state.Threshold_vals = thresholdset
            with open((os.path.expanduser(st.session_state.outputfile)+'_Chosen_Thresholds.txt'), 'w') as f:
                f.write(st.session_state.Threshold_vals)
            
            # Replicate reproducibility plot
            info2 = st.info("Producing replicate reproducibility plots...")
            def pearsonr_ci(x, y, ci=95, n_boots=100):
                x = np.asarray(x)
                y = np.asarray(y)
                
                # (n_boots, n_observations) paired arrays
                rand_ixs = np.random.randint(0, x.shape[0], size=(n_boots, x.shape[0]))
                x_boots = x[rand_ixs]
                y_boots = y[rand_ixs]
                
                # differences from mean
                x_mdiffs = x_boots - x_boots.mean(axis=1)[:, None]
                y_mdiffs = y_boots - y_boots.mean(axis=1)[:, None]
                
                # sums of squares
                x_ss = np.einsum('ij, ij -> i', x_mdiffs, x_mdiffs)
                y_ss = np.einsum('ij, ij -> i', y_mdiffs, y_mdiffs)
                
                # pearson correlations
                r_boots = np.einsum('ij, ij -> i', x_mdiffs, y_mdiffs) / np.sqrt(x_ss * y_ss)
                
                # upper and lower bounds for confidence interval
                ci_low = np.percentile(r_boots, (100 - ci) / 2)
                ci_high = np.percentile(r_boots, (ci + 100) / 2)
                return ci_low, ci_high
            def replicate_reproducibility(inputnorm,xlimt,ylimt):
                # Reads input file
                rep = inputnorm
                # Finds the unique set of replicate within the dataset
                list1 = {x[2] for x in rep.columns}
                # melts dataset to have each value in single row
                melted_rep = rep.melt(ignore_index=False).reset_index()
                #produces ID from the multiheader column names that were converted by melt into different columns (exluding replicate letter) and the row and column coordinates of each mutant e.g. Sourceplatenumber_condition_row_column
                melted_rep['ID'] = (melted_rep["variable_0"].astype(str)+"_"+melted_rep["variable_1"]+"_"+melted_rep["row"].astype(str)+"_"+melted_rep["column"].astype(str))
                # makes a list of unique pairs of replicates
                combo = [",".join(map(str, comb)) for comb in combinations(list1, 2)]
                # creates new empty dataframe with columns for the ID, the two replicate scores, the replicate letter for replicate 1 and 2.
                df_norm = pd.DataFrame(columns = ["ID","Replicate_1","Replicate_2","Rep_ID1","Rep_ID2"])
                for a in combo:
                    #takes the subset of the melted dataframe that is for the first replicate in the pair
                    rep1 = melted_rep[melted_rep['variable_2'] == a.split(",")[0]]
                    #takes the subset of the melted dataframe that is for the second replicate in the pair
                    rep2 = melted_rep[melted_rep['variable_2'] == a.split(",")[1]]
                    # takes just the ID and value columns from this subset and resets the index
                    rep1 = rep1[['ID',"value"]].reset_index(drop=True)
                    # renames value column as Replicate_1
                    rep1 = rep1.rename(columns={"value": "Replicate_1"})
                    # takes just the ID and value columns from this subset and resets the index
                    rep2 = rep2[['ID',"value"]].reset_index(drop=True)
                    # renames value column as Replicate_2
                    rep2 = rep2.rename(columns={"value": "Replicate_2"})
                    #joins the two subsets togther matching values based on the ID
                    inner_join = pd.merge(rep1, rep2, on ='ID', how ='inner')
                    #assigns the the replicate letters 
                    inner_join["Rep_ID1"] = a.split(",")[0]
                    inner_join["Rep_ID2"] = a.split(",")[1]
                    #concatenates the values to the premade empty dataframe, adding them vertically
                    df_norm = pd.concat([df_norm,inner_join], ignore_index=True)
                    #removes rows with NAs
                    df_norm = df_norm[~df_norm['Replicate_1'].isna()]
                    df_norm = df_norm[~df_norm['Replicate_2'].isna()]
                #sets the colours for the plot
                white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
                    (0, '#ffffff'),
                    (1e-200, '#440053'),
                    (0.001, '#404388'),
                    (0.01, '#2a788e'),
                    (0.1, '#21a784'),
                    (0.3, '#fde624'),
                    (1, 'r'),
                ], N=65535)

                #produces plot
                def using_mpl_scatter_density(fig, x, y):

                    b = fig.add_subplot(1, 1, 1, projection='scatter_density')
                    density = b.scatter_density(x, y, cmap=white_viridis)
                    fig.colorbar(density, label='Number of points per pixel')
                    #calculates correlation coefficient (r)
                    res = pearsonr(list(df_norm["Replicate_1"]),list(df_norm["Replicate_2"]))
                    pci = pearsonr_ci(list(df_norm["Replicate_1"]),list(df_norm["Replicate_2"]), ci=95, n_boots=100)
                    #adds r to the plot
                    b.text(0,ylimt+(ylimt*0.3),("$r = $"+str(round(res[0],3))))
                    b.text(0,ylimt+(ylimt*0.2),("$p value = $"+str('{:.3g}'.format((res[1])))))
                    b.text(0,ylimt+(ylimt*0.1),("95% CI = "+str(round(pci[0],3))+", "+str(round(pci[1],4))))
                sns.set_theme(style="white")
                fig = plt.figure()
                using_mpl_scatter_density(fig, x=df_norm["Replicate_1"], y=df_norm["Replicate_2"])
                plt.xlim(0,xlimt)
                plt.ylim(0,ylimt)
                plt.xlabel("Replicate 1")
                plt.ylabel("Replicate 2")
                return fig
            
            xlimit = st.session_state.normalised_dataset.max().max()
            fig1 = replicate_reproducibility(st.session_state.normalised_dataset,xlimit,xlimit)
            fig2 = replicate_reproducibility(st.session_state.normalised_curated_dataset,xlimit,xlimit)
            my_expander5 = st.expander(label="Replicate Reproducibility Plots:")
            colA1, colA2 = my_expander5.columns((1,1))
            colA1.write("Non-curated")
            colA1.pyplot(fig1)
            colA2.write("Curated")
            colA2.pyplot(fig2)
            info2.empty()
            st.success("Please Continue to Step 3!")
            