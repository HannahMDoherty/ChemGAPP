#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import os
import numpy as np
import pandas as pd
import scipy.stats
from scipy.stats import wilcoxon
from scipy.stats import ranksums
import re
import seaborn as sns, matplotlib.pyplot as plt
from bioinfokit.analys import stat
from statannotations.Annotator import Annotator

parser = argparse.ArgumentParser(description="Analyses small scale chemical genomic screen data",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-p", "--PATH", help="Path to folder which contains IRIS files")
parser.add_argument("-o", "--outputfile_prefix", help="Path and prefix for output file")
parser.add_argument("-pf", "--PlateInfoPath", help="The path to the folder containing the plate info files.")
parser.add_argument("-m", "--max_colony_size", help="Maximum colony size allowed, any colony larger than this will be set to this maximum",default=False,type=int)
parser.add_argument("-wt", "--WildType", help="Name of wild type strain within plate info file.")
parser.add_argument("-it", "--IRIS_type", help="Input IRIS morphology to test. Options: size,circularity,opacity", default="size")
parser.add_argument("-col_plot", "--colourpalette", help="Name of Seaborn colour palette to use for the bar and swarm plots.", default="GnBu")
parser.add_argument("-col_heat", "--colourheatmap", help="Name of Seaborn colour palette to use for the heatmap.", default="bwr_r")
parser.add_argument("-wd", "--width", help="Figure width to use for the graphs.", default=5,type=int)
parser.add_argument("-ht", "--height", help="Figure height to use for the graphs.", default=5,type=int)
parser.add_argument("-r", "--rotation", help="X Axis label rotation", default=90,type=int)
parser.add_argument("-cs", "--CircleSize", help="SwarmPlot circle size", default=2.5,type=float)
parser.add_argument("-g", "--group", help="Group bar plots by strain or condition. Options = strain, condition.", default="condition")
parser.add_argument("-pt", "--PlotType", help="Type of Plot. Options: barplot, swarmplot ", default="barplot")
parser.add_argument("-rm","--remove_strain", help="txt file of strain names to remove separated by ';'. Names must match those in plate information file. E.g. mutant1;mutant2;mutant4",default=None)
parser.add_argument("-ymax","--y_max", help="Maximum limit for y axis",default=None,type=float)
parser.add_argument("-ymin","--y_min", help="Minimum limit for y axis",default=None,type=float)
args = vars(parser.parse_args())
PATH1 = args["PATH"]
outputfile1 = args["outputfile_prefix"]
PlateInfo1 = args["PlateInfoPath"]
maximum_size = args["max_colony_size"]
wildtype1 = args["WildType"]
colour_palette = args["colourpalette"]
colour_heatmap = args["colourheatmap"]
IRIS1 = args["IRIS_type"]
rote1 = args["rotation"]
ht1 = args["height"]
wd1 = args["width"]
gr1 = args["group"]
circ_size1 = args["CircleSize"]
plottype1 = args["PlotType"]
removestrains = args["remove_strain"]
ymaximum = args["y_max"]
yminimum = args["y_min"]

def ChemGAPP_Small(PATH,outputfile,PlateInfo,WT,IRIS,width1,height1,rote,group1,circ_size,plottype):  
    m = None
    indir = os.path.expanduser(PATH)
    for batch in sorted(os.listdir(indir)):
        if batch.startswith('Batch'):
            for f in sorted(os.listdir(os.path.join(indir, batch))):
                if f.endswith(".iris"): 
                    #sys.stderr.write(batch + ' ' + f + '\n')
                    #print(os.path.join(indir, batch, f))
                    g = pd.read_csv(os.path.join(indir, batch, f),
                            comment='#',
                            index_col=[0, 1],
                            sep='\t')
                    if m is None:
                        try:
                            m = g[IRIS]
                        except:
                            m = g[IRIS]
                        m.name = (int(f.split('-')[2].split('_')[0]),
                                  '-'.join(f.split('.')[0].split('-')[:2]).replace(",","."),
                                  f.split('.')[0].split('_')[1],
                                  batch)
                        m = m.to_frame()
                    else:
                        try:
                            m1 = g[IRIS]
                        except:
                            m1 = g[IRIS]
                        m1.name = (int(f.split('-')[2].split('_')[0]),
                                  '-'.join(f.split('.')[0].split('-')[:2]).replace(",","."),
                                  f.split('.')[0].split('_')[1],
                                  batch)
                        m1 = m1.to_frame()
                        #sets them in a dataframe grouped by the condition and then columns as the A B C D or E etc
                        m = m.join(m1, how='inner')
    #m.to_csv(outputfile, sep='\t')
    m.to_csv((outputfile+"_"+IRIS+"_Initial_dataset.csv"))
    m = m.apply(pd.to_numeric)
    #m = m.replace({0:np.nan})
    m = m[sorted(m)]
    m_array = np.array(m)
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
    #n = n.replace({0:np.nan})
    n = n[sorted(n)]
    m_ind = m.reset_index()
    mask = (m_ind['row'] != 1) & (m_ind['row'] != rlen) & (m_ind['column'] != 1) & (m_ind['column'] != clen) & (m_ind['row'] != 2) & (m_ind['row'] != (rlen-1)) & (m_ind['column'] != 2) & (m_ind['column'] != (clen-1))
    df_pmm = m_ind[mask]
    df_pmm = df_pmm.set_index(['row','column'])
    df_outer = m_ind[-mask]
    df_outer = df_outer.set_index(['row','column'])
    pmm_array = np.array(df_pmm)
    m_medain = np.nanmedian(m_array)
    df_outer_norm = pd.DataFrame(index=m.index)
    
    n_array = np.array(n)
    ardf = np.zeros((mlen, 1))
    ardf.shape = (mlen,1)
    rounds = 0
    for c1,c2,i in zip(sorted(df_pmm.columns),sorted(df_outer.columns), range(len(m.columns))):
        rounds = rounds + 1
        #print(rounds)
        ar1 = pmm_array[:,i]
        ar2 = n_array[:,i]
        pmm_40 = np.nanpercentile(ar1, 40)
        pmm_60 = np.nanpercentile(ar1, 60)
        mask2= (ar1 >= pmm_40) & (ar1 <= pmm_60) 
        pmm_perc_values = ar1[mask2]
        PMM = pmm_perc_values.mean()
        df1 = df_pmm.xs((c1), axis =1, drop_level=False)
        arA=np.array(df1)
        arA= np.array(arA).flatten()
        df2 = df_outer.xs((c2), axis =1, drop_level=False)
        arB=np.array(df2)
        arB= np.array(arB).flatten()
        w, p = ranksums(arA, arB)
        if p < 0.05:
            #print("Different Dist")
            for ind, j in zip(m.index,range(len(ar2))):
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
                    ar2[j] = (((ar2[j])*(PMM/p_median))*(m_medain/PMM))
                    #if ar2[j] < 1:
                    #    ar2[j] == "X"
                    if maximum_size != False:
                        if ar2[j] > maximum_size:
                            ar2[j] = maximum_size
                elif ind[0] != 1 or ind[0] != rlen or ind[1] != 1 or ind[1] != clen or ind[0] != 2 or ind[0] != (rlen-1) or ind[1] != 2 or ind[1] != (clen-1):
                    ar2[j] = ((ar2[j])*(m_medain/PMM))
                    #if ar2[j] < 1:
                    #    ar2[j] == "X"
                    if maximum_size != False:
                        if ar2[j] > maximum_size:
                            ar2[j] = maximum_size
            ar2.shape = (mlen,1)
            ardf = np.concatenate((ardf,ar2), axis=1)
        else:
            #print(c1,"Same Dist")
            for ind, j in zip(m.index,range(len(ar2))):
                ar2[j] = ((ar2[j])*(m_medain/PMM))
                #if ar2[j] < 1:
                #    ar2[j] == "X"
                if maximum_size != False:
                    if ar2[j] > maximum_size:
                        ar2[j] = maximum_size
            ar2.shape = (mlen,1)
            ardf = np.concatenate((ardf,ar2), axis=1)
    ardf = pd.DataFrame(ardf, index=m.index)
    ardf = ardf.iloc[: , 1:]
    ardf.columns = (pd.MultiIndex.from_tuples(sorted(m.columns)))
    ardf = ardf[sorted(ardf)]
    #ardf = ardf.replace({0:pd.NA})    
    ardf.to_csv(outputfile+"_"+IRIS+"_Normalised_dataset.csv")
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
    # assign directory
    directory = os.path.expanduser(PlateInfo)
    # iterate over files in
    # that directory
    files = []
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if f.endswith(".txt"):
            if os.path.isfile(f):
                ##print(f)
                files.append(f)

    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        return [ atoi(c) for c in re.split(r'(\d+)', text) ]

    files.sort(key=natural_keys)
    plate_DFs = []
    for i in files:
        p = plate_info(i)
        plate_DFs.append(p)
    plates = {x[0] for x in ardf.columns}
    nm4 = ardf.copy(deep=False)
    columns1 = {x[1:4] for x in nm4.columns}
    df_with_strains = pd.DataFrame()
    for a, n in zip(plate_DFs, sorted(plates)):
        df1 = (ardf.xs((n), axis=1, drop_level=True))
        df1.columns = [f'{i}_{j}_{l}' for i,j,l in df1.columns]
        df2 = pd.merge(df1,a,left_index=True, right_index=True)
        df_with_strains = pd.concat([df_with_strains , df2], ignore_index=True)
    df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
    df_with_strains = df_with_strains.set_index('Gene')
    df_with_strains.columns = df_with_strains.columns.str.split('_', expand=True)
    df_with_strains = df_with_strains.sort_index()
    df_with_strains= df_with_strains.reset_index()
    df_with_strains = df_with_strains[df_with_strains['Gene'] != "EMPTY"]
    if removestrains != None:
        with open(removestrains) as f:
            lines = f.readlines()
            for line in lines:
                inclstrains = str(line).split(";")
        for stra in inclstrains:
            print(stra)
            df_with_strains = df_with_strains[df_with_strains['Gene'] != stra]
    df_with_strains = df_with_strains.set_index(['Gene'])
    df_with_strains.to_csv(outputfile+"_"+IRIS+"_Final_dataset.csv")
    #print(df_with_strains)
    group = df_with_strains.groupby(level=0).mean()
    group2 = group.copy(deep=False)
    df3 = pd.DataFrame(index=group.index)
    #print(df3)
    for column in group:
        df1 = group[column]
        df1 = pd.DataFrame(df1)
        ar1 = np.array(df1)
        ar2 = np.array(df1)
        wt_mean = float(df1.loc[WT])
        for ind, j in zip(df1.index,range(len(ar1))):
            ar2[j] = (float(ar1[j])/wt_mean)
        df3 = np.concatenate((df3,ar2), axis=1)    
    df3 = pd.DataFrame(df3, index=group.index, columns=group.columns)
    df3.index.name = None
    df3 = df3.drop(WT,axis=0)
    df3.to_csv(outputfile+"_"+IRIS+"_Scored_Dataset.csv")
    if group1 == "condition":
        if plottype == "barplot":
            conditions = {x[0] for x in df3.columns}
            for c in sorted(conditions):
                df1 = df3.xs((c), axis =1, drop_level=False)
                df2 = df1.melt(ignore_index=False)
                df2= df2.reset_index()
                df2.columns = ["Strain","Condition","Replicate","Batch","Score"]
                sns.set(style="white")
                sns.set_context("paper")
                fig, ax = plt.subplots()
                ax = sns.barplot(x="Strain", y="Score",data=df2,capsize=.5,errwidth=0.8,palette=colour_palette,edgecolor="0.2")
                ax = sns.swarmplot(x="Strain", y="Score",color="0", data=df2, alpha=1,size=3)
                ax.set(xlabel= "",ylabel='Fitness Ratio')
                if ymaximum != None or yminimum != None:
                    ax.set_ylim(yminimum,ymaximum)
                plt.title(c.replace(",","."))
                plt.xticks(rotation=rote)
                fig.set_size_inches(width1, height1)
                fig.set_dpi(1000)
                sns.despine()
                plt.savefig(outputfile1+"_"+c.replace(",",".")+"_"+IRIS+"_Bar_Plot.pdf", bbox_inches='tight')
                plt.clf()
        if plottype == "swarmplot":
            df_swarm1 = df_with_strains.reset_index()
            df_swarm2 = df_swarm1[df_swarm1['Gene'] == WT]
            #df_swarm3 = df_swarm1[df_swarm1['Gene'] != WT]
            df_swarm1['Gene'] = df_swarm1['Gene'].str.replace(WT,("0_"+WT))
            df_swarm3 = df_swarm1.sort_values("Gene") 
            wt_mean = np.array(df_swarm2.mean())
            df_swarm3 = df_swarm3.set_index("Gene")
            df_swarm4 = pd.DataFrame(index=df_swarm3.index)
            for column,mean in zip(df_swarm3,wt_mean):
                df1 = df_swarm3[column]
                df1 = pd.DataFrame(df1)
                ar1 = np.array(df1)
                ar2 = np.array(df1)
                for ind, j in zip(df1.index,range(len(ar1))):
                    ar2[j] = (float(ar1[j])/mean)
                df_swarm4 = np.concatenate((df_swarm4,ar2), axis=1)    
            df_swarm4 = pd.DataFrame(df_swarm4, index=df_swarm3.index, columns=df_swarm3.columns)
            df_swarm4.index.name = None
            conditions = {x[0] for x in df_swarm4.columns}
            for c,i in zip(sorted(conditions),range(len(conditions))):
                df1 = df_swarm4.xs((c), axis =1, drop_level=False)
                df2 = df1.melt(ignore_index=False)
                df2= df2.reset_index()
                df2.columns = ["Strain","Condition","Replicate","Batch","Score"]
                res = stat()
                res.tukey_hsd(df=df2, res_var='Score', xfac_var='Strain', anova_model='Score ~ C(Strain)')
                stats1 = res.tukey_summary[res.tukey_summary['group1'] == ("0_"+WT)]
                pairs = list(zip(stats1.group2, stats1.group2))
                pvalues = list(stats1['p-value'])
                df2 = df2[df2.Strain != ("0_"+WT)]
                subcat_palette = sns.dark_palette("#8BF", reverse=True, n_colors=5)
                plotting_parameters = {
                    'data':    df2,
                    'x':       'Strain',
                    'y':       'Score'}
                sns.set(style="white")
                sns.set_context("paper")
                fig, ax = plt.subplots()
                ax = sns.swarmplot(x="Strain", y="Score",hue="Strain", data=df2, alpha=1,size=circ_size,palette=colour_palette)
                ax = sns.boxplot(showmeans=True,
                        meanline=True,
                        meanprops={'color': 'k', 'ls': '-', 'lw': 1},
                        medianprops={'visible': False},
                        whiskerprops={'visible': False},
                        zorder=10,
                        x="Strain", y="Score",data=df2,
                        showfliers=False,
                        showbox=False,
                        showcaps=False,
                        ax=ax)
                annotator = Annotator(ax, pairs, **plotting_parameters)
                annotator.set_pvalues(pvalues)
                annotator.configure(line_width = 0)
                annotator.annotate()
                for patch in ax.artists:
                    r, g, b, a = patch.get_facecolor()
                    patch.set_facecolor((r, g, b, 0))    
                ax.set(xlabel= "Strain",ylabel='Fitness Ratio')
                if ymaximum != None or yminimum != None:
                    ax.set_ylim(yminimum,ymaximum)
                plt.title(c.replace(",",".").replace("-"," "))
                plt.xticks(rotation=rote)
                fig.set_size_inches(width1, height1)
                fig.set_dpi(1000)
                sns.despine()
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                              fancybox=True, shadow=False, ncol=1)
                plt.savefig(outputfile1+"_"+c.replace(",",".")+"_"+IRIS+"_Swarm_Plot.pdf", bbox_inches='tight')
                plt.clf()
    if group1 == "strain":
        if plottype == "barplot":
            b = df3.T
            conditions = {x for x in b.columns}
            for c in sorted(conditions):
                df1 = b.xs((c), axis=1, drop_level=False)
                df1 = df1.T
                df2 = df1.melt(ignore_index=False)
                df2.columns = ["Condition","Replicate","Batch","Score"]
                #for j in range(len(df2.Condition)):
                #    df2.Condition[j] = df2.Condition[j].replace(",",".")
                sns.set(style="white")
                sns.set_context("paper")
                fig, ax = plt.subplots()
                ax = sns.barplot(x="Condition", y="Score",data=df2,capsize=.5,errwidth=0.8,palette=colour_palette,edgecolor="0.2")
                ax = sns.swarmplot(x="Condition", y="Score",color="0", data=df2, alpha=1,size=3)
                ax.set(xlabel= "",ylabel='Fitness Ratio')
                if ymaximum != None or yminimum != None:
                    ax.set_ylim(yminimum,ymaximum)
                plt.xticks(rotation=rote)
                plt.title(c)
                fig.set_size_inches(width1, height1)
                fig.set_dpi(1000)
                sns.despine()
                plt.savefig(outputfile1+"_"+c+"_"+IRIS+"_Bar_Plot.pdf", bbox_inches='tight')
                plt.clf()
        if plottype == "swarmplot":
            df_swarm1 = df_with_strains.reset_index()
            df_swarm2 = df_swarm1[df_swarm1['Gene'] == WT]
            #df_swarm3 = df_swarm1[df_swarm1['Gene'] != WT]
            df_swarm1['Gene'] = df_swarm1['Gene'].str.replace(WT,("0_"+WT))
            df_swarm3 = df_swarm1.sort_values("Gene") 
            wt_mean = np.array(df_swarm2.mean())
            df_swarm3 = df_swarm3.set_index("Gene")
            df_swarm4 = pd.DataFrame(index=df_swarm3.index)
            for column,mean in zip(df_swarm3,wt_mean):
                df1 = df_swarm3[column]
                df1 = pd.DataFrame(df1)
                ar1 = np.array(df1)
                ar2 = np.array(df1)
                for ind, j in zip(df1.index,range(len(ar1))):
                    ar2[j] = (float(ar1[j])/mean)
                df_swarm4 = np.concatenate((df_swarm4,ar2), axis=1)    
            df_swarm4 = pd.DataFrame(df_swarm4, index=df_swarm3.index, columns=df_swarm3.columns)
            df_swarm4.index.name = None
            
            conditions = {x[0] for x in df_swarm4.columns}
            genes = {x for x in df_swarm4.index}
            leng = len(genes)-1
            anovas = np.zeros((leng, 3),dtype=object)
            anovas.shape = (leng,3)
            for c in sorted(conditions):
                df1 = df_swarm4.xs((c), axis =1, drop_level=False)
                df2 = df1.melt(ignore_index=False)
                df2 = df2.reset_index()
                df2.columns = ["Strain","Condition","Replicate","Batch","Score"]
                res = stat()
                res.tukey_hsd(df=df2, res_var='Score', xfac_var='Strain', anova_model='Score ~ C(Strain)')
                stats1 = res.tukey_summary[res.tukey_summary['group1'] == ("0_"+WT)]
                #print(stats1)
                cond = ([c] * leng)
                pairs = list(zip(stats1.group2, stats1.group2))
                pvalues = list(stats1['p-value'])
                grouped = list(zip(pairs,cond,pvalues))
                #print(grouped)
                grouped = np.array(grouped,dtype=object)
                anovas = np.append(anovas,grouped, axis=0)
            #print(anovas)
            anovas = anovas[leng:len(anovas)]
            df_swarm4 = df_swarm4[df_swarm4.index != ("0_"+WT)]
            b = df_swarm4.T
            conditions = {x for x in b.columns}
            for c,i in zip(sorted(conditions),range(len(conditions))):
                df1 = b.xs((c), axis=1, drop_level=False)
                df1 = df1.T
                df2 = df1.melt(ignore_index=False)
                df2.columns = ["Condition","Replicate","Batch","Score"]
                m = [row for row in anovas if c == row[0][0]]
                pair1 = [(row[1],row[1]) for row in m]
                pvalue1 = [row[2] for row in m]
                #print(pair1,pvalue1)
                df2 = df2[df2.Condition != "WT"]
                subcat_palette = sns.dark_palette("#8BF", reverse=True, n_colors=5)
                plotting_parameters = {
                'data':    df2,
                'x':       'Condition',
                'y':       'Score',
                'palette': subcat_palette[1:]}
                sns.set(style="white")
                sns.set_context("paper")
                fig, ax = plt.subplots()
                ax = sns.swarmplot(x="Condition", y="Score",hue="Condition", data=df2, alpha=1,size=circ_size,palette=colour_palette)
                ax = sns.boxplot(showmeans=True,
                        meanline=True,
                        meanprops={'color': 'k', 'ls': '-', 'lw': 1},
                        medianprops={'visible': False},
                        whiskerprops={'visible': False},
                        zorder=10,
                        x="Condition", y="Score",data=df2,
                        showfliers=False,
                        showbox=False,
                        showcaps=False,
                        ax=ax)
                annotator = Annotator(ax, pair1, **plotting_parameters)
                annotator.set_pvalues(pvalue1)
                annotator.configure(line_width = 0)
                annotator.annotate()
                for patch in ax.artists:
                    r, g, b, a = patch.get_facecolor()
                    patch.set_facecolor((r, g, b, 0))    
                ax.set(xlabel= "Condition",ylabel='Fitness Ratio')
                if ymaximum != None or yminimum != None:
                    ax.set_ylim(yminimum,ymaximum)
                plt.title(c.replace(",",".").replace("-"," "))
                plt.xticks(rotation=rote)
                fig.set_size_inches(width1, height1)
                fig.set_dpi(1000)
                sns.despine()
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                              fancybox=True, shadow=False, ncol=1)
                plt.savefig(outputfile1+"_"+c+"_"+IRIS+"_Swarm_Plot.pdf", bbox_inches='tight')
                plt.clf()
    b = df3.T
    c = b.groupby(level=0).mean()
    x = c.T
    for i in range(len(x.columns)):
        x = x.rename(columns={x.columns[i]: x.columns[i].replace(",",".").replace("-"," ")})
    c = x.T
    sns.set(style="white")
    sns.set_context("paper")
    colp = sns.color_palette(colour_heatmap, as_cmap=True)
    #vmax1 = c.max().max()
    fig, ax = plt.subplots()
    ax = sns.heatmap(c,center=1,cmap=colp,square=False,annot=True,annot_kws={"size": 4},fmt=".4f",linewidths=.1,linecolor='0',vmin=0, vmax=2)
    fig.set_size_inches(width1, height1)
    plt.savefig(outputfile1+"_"+IRIS+"_Heatmap.pdf", bbox_inches='tight')
    plt.clf()
    return df3
ChemGAPP_Small(PATH1,outputfile1,PlateInfo1,wildtype1,IRIS1,wd1,ht1,rote1,gr1,circ_size1,plottype1)