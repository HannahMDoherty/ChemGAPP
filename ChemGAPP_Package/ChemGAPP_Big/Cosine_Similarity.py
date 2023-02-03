#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import pandas as pd
import random
from scipy import spatial
import re
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as  mpatches
from sklearn.metrics import auc
import numpy as np
import os
import itertools as it 
def get_options():
    parser = argparse.ArgumentParser(description="Calculates the cosine similarity scores for the phenotypic profiles of genes from the same operon and genes from different operons. Produces a density plot of the cosine similarity scores for genes of the same and different operons. Produces an ROC curve testing models ability at different threshold. ",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--InputFile", help="The dataset with gene names added. Output from Add_gene_names.py")
    parser.add_argument("-o", "--OutputFile", help="List of genes compared and the cosine similarity score as well as if they belong to the same operon")
    parser.add_argument("-or", "--Output_ROC_curve", help="Plot of the ROC curve and AUC score.")
    parser.add_argument("-od", "--Output_Density_plot", help="Density plot of the cosine similarity scores for same and different operons.")
    parser.add_argument("-clus", "--Cluster_file", help="A CSV file containing the operon clusters for each gene within the bacterium of interest, where columns = (Cluster,Gene). The Genes must match the names assigned to the scored dataset.")
    return parser.parse_args()

def main():
    options = get_options()
    inpt = options.InputFile
    outpt = options.OutputFile
    Op_Clus = options.Cluster_file
    outroc = options.Output_ROC_curve
    outdens = options.Output_Density_plot
    outroc = os.path.expanduser(outroc)
    outdens = os.path.expanduser(outdens)
    #Open the operon cluster and s-score files
    df = pd.read_csv(Op_Clus)
    df = df[['Gene', 'Cluster']]
    ar1 = np.array(df)
    df_14 = pd.read_table(inpt, index_col=0)
    df_14 = df_14.reset_index()
    df_14 = df_14.rename(columns={"index": "Gene"})
    df_14 = df_14.sort_values('Gene')
    df_14 = df_14[~df_14.Gene.str.contains("empty")]
    df_14 = df_14[~df_14.Gene.str.contains("EMPTY")]
    df_14 = df_14.reset_index(drop=True)
    ar_14 = np.array(df_14)
    
    #Collect only the clusters for genes in the dataset and assigning clusters to the dataset in end column
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
    # Produce density plot for similarity scores between same and different operons
    sns.set(style="darkgrid") 
    fig = sns.kdeplot(df_T_F_2['FALSE'], shade=True, color="r")
    fig = sns.kdeplot(df_T_F_2['TRUE'], shade=True, color="b")
    plt.xlabel("Cosine Similarity Score")
    plt.ylabel("Density")
    handles = [mpatches.Patch(facecolor="r", label="Different Operon"), mpatches.Patch(facecolor="b", label="Same Operon")]
    plt.legend(handles=handles)
    fig.set_title('Density Plot of Cosine Similarity Scores Between Genes of the Same and Different Operons')
    plt.savefig(outdens, bbox_inches='tight')
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
    AUC_txt = f'AUC {auc(spec_rate, sens_rate)}'
    fig.suptitle(AUC_txt, fontsize=12, fontweight='bold')
    plt.savefig(outroc, bbox_inches='tight')

if __name__ == "__main__":
    main()