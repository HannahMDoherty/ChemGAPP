import streamlit as st
import pandas as pd
import random
from scipy import spatial
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as  mpatches
from sklearn.metrics import auc
import numpy as np
import os
import itertools as it
st.set_page_config(layout="wide")
#st.markdown("# Main page ")
#st.sidebar.markdown("# Main page")
st.title('ChemGAPP Big: Dataset Comparison')
my_expandera = st.expander(label="Instructions:", expanded=True)
my_expandera.markdown(""" 
1- Upload your cluster file. 

File should be in CSV format and consist of two columns; Cluster and Gene. 

The Genes must match the names assigned to the scored dataset in the step 3.


E.g: 

| Cluster | Gene     |                
|  ------ | :------: |      
| 1 | 1_2_PA14_00050 |      
| 2 | 1_3_PA14_00060 |      
| 2 | 1_4_PA14_00070 |      
| 4 | 1_5_PA14_00080 |      
| 5 | 1_6_PA14_00090 |      

Matching to scored dataset:

| Gene | Condition A    | Condition B| 
|  ------| ------ | :------: |     
| 1_2_PA14_00050 | 0.12 | 1.23 |     
| 1_3_PA14_00060 | 0.54 | 1.53 |     
| 1_4_PA14_00070 | 5.62 | -1.23 |     
| 1_5_PA14_00080 | -10.12 | 8.25 |     
| 1_6_PA14_00090 | 0.16 | -5.73 |   

------------------

2- Select if you want to analyse the averaged or non-averaged datasets from Step 3.

---------------

3- Press Begin!
""")

col1, col2 = st.columns((1,1))
gene_names = st.sidebar.file_uploader("Upload Cluster File:", accept_multiple_files=False)
Ave_datasets = st.sidebar.radio(
    "Which datasets would you like to compare:",
    ('Averaged','Non-Averaged'))
complete = st.sidebar.button(label="Begin!")
if complete:
    if gene_names:
        gene_names2 = pd.read_csv(gene_names)
    if "final_DF" in st.session_state and "curated_DF" in st.session_state:

        my_expander2 = st.expander(label="Chosen Thresholds:")
        my_expander2.write(st.session_state.Threshold_vals)
        if gene_names:
            element1 = st.info("Comparing Datasets")
            def Cosine_Similarity(Op_Clus,inpt,outpt,outroc,outdens):
                outroc = os.path.expanduser(outroc)
                outdens = os.path.expanduser(outdens)
                #Open the operon cluster and s-score files
                df = Op_Clus
                df = df[['Gene', 'Cluster']]
                ar1 = np.array(df)
                df_14 = inpt
                df_14 = df_14.reset_index()
                #df_14['Gene'] = [x.split('_')[0].split(' ')[0] for x in df_14["Gene"]]
                df_14 = df_14.rename(columns={"index": "Gene"})
                df_14 = df_14.sort_values('Gene')
                df_14 = df_14[~df_14.Gene.str.contains("empty")]
                df_14 = df_14[~df_14.Gene.str.contains("EMPTY")]
                df_14 = df_14.reset_index(drop=True)
                ar_14 = np.array(df_14)
                #Collect only the clusters for genes in the dataset
                df_Clus = pd.DataFrame(columns=df.columns)

                for i in range(len(ar1)): 
                    x = ar1[i][0]
                    if str(x) in ar_14[:,0]:
                        rows = df.loc[i]
                        df_Clus = df_Clus.append(rows, ignore_index=True)
                df_merge = pd.merge(df_14,df_Clus)
                df_co_sim_R_T = pd.DataFrame()
                ar_merge = np.array(df_merge)
                mainlist = list(range(0,len(df_merge)))
                for i,j in it.combinations(mainlist, 2):
                    #ignores comparison with self
                    if i != j:
                        #if in same cluster list all s-scores
                        if str(ar_merge[i][-1]) == str(ar_merge[j][-1]):
                            p = list(ar_merge[i,1:-1])
                            p2 = list(ar_merge[j,1:-1])
                            data1 = []
                            data2 = []
                            #ensures nans are removed as well as the corresponding value for the other list.
                            for itemp,indexp,itemp2, in zip(p,range(len(p)),p2):
                                if str(itemp) != "nan" and str(itemp2) != "nan":
                                    data1.append(float(itemp))
                                    data2.append(float(itemp2))
                            n = str(ar_merge[i,0:1])
                            n2 = str(ar_merge[j,0:1])
                            #calculates cosine similarity between the phenotypic profiles of s-scores for the two genes and adds to a df.
                            cosine_similarity = 1 - spatial.distance.cosine(data1, data2)
                            name = np.array([n, n2, cosine_similarity, "TRUE"],dtype=object)
                            df_co_sim_R_T = df_co_sim_R_T.append(pd.DataFrame(name).T)
                df_co_sim_R_T.columns=['Gene1','Gene2','Cosine_score','Same Operon']
                df_co_sim_R_T = df_co_sim_R_T.reset_index(drop=True)

                random.seed(2)
                def random_pairs(number_list): 
                    return [number_list[i] for i in random.sample(range(len(number_list)), 2)] 
                #selects twice as many pairs as needed for different operons as likely some
    # will come up as the same and want to later syphon down to same number of comparisions as for same operon.
                n = (2*len(df_co_sim_R_T))
                numbers = list(range(len(df_merge)))
                rp = [random_pairs(numbers) for i in range(n)]
                df_co_sim_R_F = pd.DataFrame()
                for k in range(n):
                        i = rp[k][0]
                        j = rp[k][1]
                        d1 = list(ar_merge[i,-1:])
                        d2 = list(ar_merge[j,-1:])
                        #if operon clusters do not match.
                        if d1 != d2:
                            p = list(ar_merge[i,1:-1])
                            p2 = list(ar_merge[j,1:-1])
                            data1 = []
                            data2 = []
                            for itemp,indexp,itemp2, in zip(p,range(len(p)),p2):
                                if str(itemp) != "nan" and str(itemp2) != "nan":
                                    data1.append(float(itemp))
                                    data2.append(float(itemp2))
                            n = ar_merge[i,0:1]
                            n2 = ar_merge[j,0:1]
                            cosine_similarity = 1 - spatial.distance.cosine(data1, data2)
                            name = np.array([n, n2, cosine_similarity, "FALSE"],dtype=object)
                            df_co_sim_R_F = df_co_sim_R_F.append(pd.DataFrame(name).T)
                df_co_sim_R_F.columns=['Gene1','Gene2','Cosine_score','Same Operon']
                df_co_sim_R_F = df_co_sim_R_F.reset_index(drop=True)

                # drop any NA values
                df_co_sim_R_T.dropna(subset = ["Cosine_score"], inplace=True)
                df_co_sim_R_F.dropna(subset = ["Cosine_score"], inplace=True)
                df_co_sim_R_F = df_co_sim_R_F.reset_index(drop=True)
                df_co_sim_R_T = df_co_sim_R_T.reset_index(drop=True)


                if len(df_co_sim_R_F) < len(df_co_sim_R_T):
                    # pick random gene pairs from the different operons dataframe
                    # selecting same number as there are in the same operon dataframe
                    random.seed(2)
                    r_S2 = random.sample(range(len(df_co_sim_R_T)), len(df_co_sim_R_F))
                    df_S2= pd.DataFrame(columns=df_co_sim_R_T.columns)
                    for i in r_S2:
                        rows = df_co_sim_R_T.loc[i]
                        df_S2 = df_S2.append(rows, ignore_index=True)
                    df_T_F = pd.concat([df_co_sim_R_F, df_S2],ignore_index=True)
                elif len(df_co_sim_R_F) > len(df_co_sim_R_T):
                    # pick random gene pairs from the different operons dataframe
                    # selecting same number as there are in the same operon dataframe
                    random.seed(2)
                    r_S2 = random.sample(range(len(df_co_sim_R_F)), len(df_co_sim_R_T))
                    df_S2= pd.DataFrame(columns=df_co_sim_R_F.columns)
                    for i in r_S2:
                        rows = df_co_sim_R_F.loc[i]
                        df_S2 = df_S2.append(rows, ignore_index=True)
                # Join the subsetted different operon table and same operon table
                    df_T_F = pd.concat([df_co_sim_R_T, df_S2],ignore_index=True)
                # Format the merged table so it can be pivoted to 
                # have cosine similarity scores put into a same or different operon column
                df_T_F['Cosine_score'] = df_T_F['Cosine_score'].astype(float)
                for f in df_T_F.columns:
                    if f != 'Cosine_score':
                        df_T_F[f] = df_T_F[f].astype(str)

                for f in list(['Gene1',"Gene2"]):
                    df_T_F[f] = df_T_F[f].str.replace("\[\'",'', regex = True)
                    df_T_F[f] = df_T_F[f].str.replace("\'\]",'', regex = True)

                for f in df_T_F.columns:
                    if f != 'Cosine_score':    
                        df_T_F[f] = pd.array(df_T_F[f].tolist())
                df_T_F.to_csv(outpt, index=False)
                df_T_F_2 = pd.pivot_table(df_T_F, values='Cosine_score', index=['Gene1','Gene2'],columns=['Same Operon'])
                df_FP = df_T_F_2.reset_index(drop=True)
                F = list(df_FP['FALSE'])
                F = [x for x in F if str(x) != 'nan']
                T = list(df_FP['TRUE'])
                T = [x for x in T if str(x) != 'nan']
                df_thres = pd.DataFrame()
                #calculate true and false positives and negatives at different thresholds for AUC and ROC curve production.
                for i in np.arange(-1, 1.1, 0.1):
                    TP = 0
                    TN = 0
                    FP = 0
                    FN = 0
                    for row in F:
                        if row < i:
                            TN = TN + 1
                        else: 
                            FP = FP + 1
                    for row in T:
                        if row < i:
                            FN = FN + 1
                        else: 
                            TP = TP + 1
                    sens = TP/(TP+FN)
                    spec = TN/(TN+FP)
                    FPR = 1-spec
                    name = np.array([i,sens,spec,FPR],dtype=object)
                    df_thres = df_thres.append(pd.DataFrame(name).T)
                df_thres.columns=['Threshold','Sensitivity' ,'Specificity','False Positive Rate']
                sens_rate = list(df_thres['Sensitivity'])
                spec_rate = list(df_thres['False Positive Rate'])
                fig, ax = plt.subplots(figsize=(6,6))
                ax.plot(spec_rate, sens_rate)
                ax.set_xlabel('False Positive Rate')
                ax.set_ylabel('True Positive Rate')
                ax.set_title('ROC curve for Assigning S-scores to Genes of Same and Different Operons')
                AUC_txt = f'AUC {round(auc(spec_rate, sens_rate),4)}'
                fig.suptitle(AUC_txt, fontsize=12, fontweight='bold')
                return df_T_F_2,fig,ax,sens_rate,spec_rate

            Op_Clus1 = gene_names2
            if Ave_datasets == 'Averaged':
                inpt1 = st.session_state.final_ave_DF
                inpt2 = st.session_state.curated_ave_DF
            if Ave_datasets == 'Non-Averaged':
                inpt1 = st.session_state.final_DF
                inpt2 = st.session_state.curated_DF
            outpt1 = (st.session_state.outputfile+"_Original_Cosine_Similarity_Scores.csv")
            outroc1 = (st.session_state.outputfile+"_Original_ROC_Curve.pdf")
            outdens1 = (st.session_state.outputfile+"_Original_Density_Plot.pdf")
            outpt2 = (st.session_state.outputfile+"_Curated_Cosine_Similarity_Scores.csv")
            outroc2 = (st.session_state.outputfile+"_Curated_ROC_Curve.pdf")
            outdens2 = (st.session_state.outputfile+"_Curated_Density_Plot.pdf")
            densdf1,fig_1,ax_1,sens_rate_1,spec_rate_1 = Cosine_Similarity(Op_Clus1,inpt1,outpt1,outroc1,outdens1)
            densdf2,fig_2,ax_2,sens_rate_2,spec_rate_2 = Cosine_Similarity(Op_Clus1,inpt2,outpt2,outroc2,outdens2)
            sns.set_theme(style="white")
            fig1 = plt.figure()
            # Produce density plot for similarity scores between same operons for the two datasets
            sns.kdeplot(densdf1['TRUE'], shade=False, color="b", label="Original Dataset")
            sns.kdeplot(densdf2['TRUE'], shade=False, color="r", label=("Curated Dataset"))
            plt.xlabel("Cosine Similarity Score")
            plt.ylabel("Density")
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.title('Density Plot of Cosine Similarity Scores for Genes of the Same Operons')
            cola,colb,colc =st.columns((1,3,1))
            colb.pyplot(fig1)
            col3a, col3, col4 = st.columns((1,3,1))
            #produce the ROC plot for the two datasets
            sns.set_theme(style="white")
            fig, ax = plt.subplots()
            ax.plot(spec_rate_1, sens_rate_1,color='r',label="Original")
            ax.plot(spec_rate_2, sens_rate_2,color='b',label="Curated")
            ax.set_xlabel('False Positive Rate')
            ax.set_ylabel('True Positive Rate')
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                                  fancybox=True, shadow=False, ncol=1)
            ax.set_title('ROC curve for Assigning S-scores to Genes of Same and Different Operons')
            AUC_txt = f'Original AUC: {round(auc(spec_rate_1, sens_rate_1),4)}'
            AUC_txt2 = f'Curated AUC: {round(auc(spec_rate_2, sens_rate_2),4)}'
            fig.text(0.15,0.8,AUC_txt,color='r')
            fig.text(0.15,0.75,AUC_txt2,color='b')
            col3.pyplot(fig)
            element1.empty()
    else:
        st.warning("Please Return to Step 3 and Select 'Both' Option")