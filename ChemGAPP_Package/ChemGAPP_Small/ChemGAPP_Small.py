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

def get_options():
    parser = argparse.ArgumentParser(description="Analyses small scale chemical genomic screen data",
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--PATH", help="Path to folder which contains IRIS files")
    parser.add_argument("-o", "--outputfile_prefix", help="Path and prefix for output file")
    parser.add_argument("-pf", "--PlateInfoPath", help="The path to the folder containing the plate info files.")
    parser.add_argument("-m", "--max_colony_size", help="Maximum colony size allowed, any colony larger than this will be removed",default=False,type=int)
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
    return parser.parse_args()


def ChemGAPP_Small():  
    options = get_options()
    PATH = options.PATH
    outputfile = options.outputfile_prefix
    PlateInfo = options.PlateInfoPath
    maximum_size = options.max_colony_size
    WT = options.WildType
    colour_palette = options.colourpalette
    colour_heatmap = options.colourheatmap
    IRIS = options.IRIS_type
    rote = options.rotation
    height1 = options.height
    width1 = options.width
    group1 = options.group
    circ_size = options.CircleSize
    plottype = options.PlotType
    removestrains = options.remove_strain
    ymaximum = options.y_max
    yminimum = options.y_min
    m = None
    indir = os.path.expanduser(PATH)
    # cycles through iris files and uses filename to produce column headers.
    for f in sorted(os.listdir(indir)):
        if f.endswith(".iris"): 
            g = pd.read_csv(os.path.join(indir, f),
                    comment='#',
                    index_col=[0, 1],
                    sep='\t')
            if m is None:
                m = g[IRIS]
                # replaces commas to periods, since these are not allowed within the filename and dashes to spaces
                m.name = (int(f.split('-')[2].split('_')[0]),
                          '-'.join(f.split('.')[0].split('-')[:2]).replace(",",".").replace("-"," "),
                          f.split('.')[0].split('_')[1])
                m = m.to_frame()
            else:
                m1 = g[IRIS]
                m1.name = (int(f.split('-')[2].split('_')[0]),
                          '-'.join(f.split('.')[0].split('-')[:2]).replace(",",".").replace("-"," "),
                          f.split('.')[0].split('_')[1])
                m1 = m1.to_frame()
                #sets files into a dataframe grouped by the condition and then columns as the A B C D or E etc
                m = m.join(m1, how='inner')
    # saves this initial dataset
    m.to_csv((outputfile+"_"+IRIS+"_Initial_dataset.csv"))
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
    
    n_array = np.array(n)
    ardf = np.zeros((mlen, 1))
    ardf.shape = (mlen,1)
    rounds = 0
    for c1,c2,i in zip(sorted(df_pmm.columns),sorted(df_outer.columns), range(len(m.columns))):
        rounds = rounds + 1
        #runs through each plate individually matching the plates for the outer and inner dataframes
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
            #adds the column of adjusted scores to the 'normalised' dataframe 
            ardf = np.concatenate((ardf,ar2), axis=1)
        # if outer edge and inner colonies are not signficantly different 
        # then just scales entire plate such that the PMM is equal to the median colony size of entire dataset
        else:
            
            for ind, j in zip(m.index,range(len(ar2))):
                ar2[j] = ((ar2[j])*(m_medain/PMM))

            ar2.shape = (mlen,1)
            ardf = np.concatenate((ardf,ar2), axis=1)
    ardf = pd.DataFrame(ardf, index=m.index)
    #removes first column of zero values
    ardf = ardf.iloc[: , 1:]
    # adds in the multiheader columns
    ardf.columns = (pd.MultiIndex.from_tuples(sorted(m.columns)))
    ardf = ardf[sorted(ardf)]  
    ardf.to_csv(outputfile+"_"+IRIS+"_Normalised_dataset.csv")

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
    # takes the specified directory for the plate information files
    directory = os.path.expanduser(PlateInfo)
    # iterate over plate files in that directory 
    files = []
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if f.endswith(".txt"):
            if os.path.isfile(f):
                ##print(f)
                files.append(f)

    # looks for digits within the file names so that order of 
    # plate info files is sorted plate1, plate2, plate3 etc. not plate1, plate10, plate2 etc.
    def atoi(text):
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        return [ atoi(c) for c in re.split(r'(\d+)', text) ]

    files.sort(key=natural_keys)
    #adds the plate files to a list
    plate_DFs = []
    for i in files:
        p = plate_info(i)
        plate_DFs.append(p)
    
    # makes set of all source plate numbers within dataset
    plates = {x[0] for x in ardf.columns}
    nm4 = ardf.copy(deep=False)
    df_with_strains = pd.DataFrame()
    # iterates through each source plate number and makes sub-dataset for each such that gene names can be appended.
    for a, n in zip(plate_DFs, sorted(plates)):
        df1 = (ardf.xs((n), axis=1, drop_level=True))
        df1.columns = [f'{i}_{j}' for i,j in df1.columns]
        #merges sub-dataset with the gene names, matching based on row and column location.
        df2 = pd.merge(df1,a,left_index=True, right_index=True)
        #adds merged dataset to new dataframe, such that all source plate named columns are within dataset.
        df_with_strains = pd.concat([df_with_strains , df2], ignore_index=True)
    df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
    df_with_strains = df_with_strains.set_index('Gene')
    # splits columns back into multiheaders
    df_with_strains.columns = df_with_strains.columns.str.split('_', expand=True)
    df_with_strains = df_with_strains.sort_index()
    if maximum_size != False:
        df_with_strains[df_with_strains > float(maximum_size)] = np.nan
    df_with_strains= df_with_strains.reset_index()
    # removes rows named "EMPTY"
    df_with_strains = df_with_strains[df_with_strains['Gene'] != "EMPTY"]
    # if specified to remove certain strains will read the input and split by ";" to remove rows with specified names.
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
    
    # For the bar plots:
    # averages the scores of rows with the same gene name.
    group = df_with_strains.groupby(level=0).mean()
    group2 = group.copy(deep=False)
    df3 = pd.DataFrame(index=group.index)
   
    # iterates through each condition
    for column in group:
        df1 = group[column]
        df1 = pd.DataFrame(df1)
        ar1 = np.array(df1)
        ar2 = np.array(df1)
        wt_mean = float(df1.loc[WT])
        # divides each averaged colony size by the mean wildtype colonysize for that condition
        for ind, j in zip(df1.index,range(len(ar1))):
            ar2[j] = (float(ar1[j])/wt_mean)
        #adds scored column to new dataframe
        df3 = np.concatenate((df3,ar2), axis=1)    
    df3 = pd.DataFrame(df3, index=group.index, columns=group.columns)
    df3.index.name = None
    # drops the WT row from the scored dataset.
    df3 = df3.drop(WT,axis=0)
    df3.to_csv(outputfile+"_"+IRIS+"_Scored_Dataset.csv")

    # if selected to group by condition, and bar plot then produces bar plot for each condition. 
    if group1 == "condition":
        if plottype == "barplot":
            conditions = {x[0] for x in df3.columns}
            for c in sorted(conditions):
                df1 = df3.xs((c), axis =1, drop_level=False)
                df2 = df1.melt(ignore_index=False)
                df2= df2.reset_index()
                df2.columns = ["Strain","Condition","Replicate","Score"]
                sns.set(style="white")
                sns.set_context("paper")
                fig, ax = plt.subplots()
                # plots bar plot and a swarm plot so the fitness scores for each replicates can be seen on the bar chart.
                ax = sns.barplot(x="Strain", y="Score",data=df2,capsize=.5,errwidth=0.8,palette=colour_palette,edgecolor="0.2")
                ax = sns.swarmplot(x="Strain", y="Score",color="0", data=df2, alpha=1,size=3)
                ax.set(xlabel= "",ylabel='Fitness Ratio')
                if ymaximum != None or yminimum != None:
                    ax.set_ylim(yminimum,ymaximum)
                plt.title(c)
                plt.xticks(rotation=rote)
                fig.set_size_inches(width1, height1)
                fig.set_dpi(1000)
                sns.despine()
                plt.savefig(outputfile1+"_"+c.replace(",",".")+"_"+IRIS+"_Bar_Plot.pdf", bbox_inches='tight')
                plt.clf()
        # if swarmplot selected produces swarmplots instead. Here the data is scored differently. 
        if plottype == "swarmplot":
            df_swarm1 = df_with_strains.reset_index()
            #makes df including just the wildtype values
            df_swarm2 = df_swarm1[df_swarm1['Gene'] == WT]
            # df_swarm3 = df_swarm1[df_swarm1['Gene'] != WT]
            # adds "0_" to WT column name such that it is always sorted first 
            # when sorting alphabetically, necessary for the ANOVA tests.
            df_swarm1['Gene'] = df_swarm1['Gene'].str.replace(WT,("0_"+WT))
            df_swarm3 = df_swarm1.sort_values("Gene") 
            # calculates the wildtype colony size mean for each condition plate.
            wt_mean = np.array(df_swarm2.mean())
            df_swarm3 = df_swarm3.set_index("Gene")
            df_swarm4 = pd.DataFrame(index=df_swarm3.index)
            # iterates through each condition plate and associated WT mean
            for column,mean in zip(df_swarm3,wt_mean):
                df1 = df_swarm3[column]
                df1 = pd.DataFrame(df1)
                ar1 = np.array(df1)
                ar2 = np.array(df1)
                #divides each individual colony size by the WT_mean, including for the WT values.
                for ind, j in zip(df1.index,range(len(ar1))):
                    ar2[j] = (float(ar1[j])/mean)
                #adds to a new dataframe
                df_swarm4 = np.concatenate((df_swarm4,ar2), axis=1)    
            df_swarm4 = pd.DataFrame(df_swarm4, index=df_swarm3.index, columns=df_swarm3.columns)
            df_swarm4.index.name = None
            conditions = {x[0] for x in df_swarm4.columns}
            #iterates through conditions and makes sub-dataset including all replicate plates for same condition.
            for c,i in zip(sorted(conditions),range(len(conditions))):
                df1 = df_swarm4.xs((c), axis =1, drop_level=False)
                #metls dataset so all scores are in sigular column
                df2 = df1.melt(ignore_index=False)
                df2= df2.reset_index()
                df2.columns = ["Strain","Condition","Replicate","Score"]
                # performs anova and tukey-hsd but only takes the values for comarison to the WT for each gene,
                # such that significance can be plotted on swarm plots
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
                # adds mean bars
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
                #adds anova asterisk type annotations for sigificance. 
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
                plt.title(c)
                plt.xticks(rotation=rote)
                fig.set_size_inches(width1, height1)
                fig.set_dpi(1000)
                sns.despine()
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                              fancybox=True, shadow=False, ncol=1)
                plt.savefig(outputfile1+"_"+c+"_"+IRIS+"_Swarm_Plot.pdf", bbox_inches='tight')
                plt.clf()
    
    #here groups by strain and not condition.
    if group1 == "strain":
        if plottype == "barplot":
            #transposes dataset such that it can be grouped by strains.
            b = df3.T
            # makes set of all strains. 
            strain_list = {x for x in b.columns}
            #splits and iterates by strain.
            for c in sorted(strain_list):
                df1 = b.xs((c), axis=1, drop_level=False)
                df1 = df1.T
                #melts dataset such that all scores in one column with key being for the conditions
                df2 = df1.melt(ignore_index=False)
                df2.columns = ["Condition","Replicate","Score"]
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
            #as above with condition grouped swarmplots, 
            # produces dataset with each individual colonzy size divided by WT mean for that condition. 
            # then calulates anova and tukey-hsd and takes p values for WT vs each strain. 
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
            # produces zeros matrix that can be filled for [(gene1,gene1),condition,p-value], 
            # need to compare pairs to themselves since plots sorted by strain then plot based on the conditions split by gene name
            anovas = np.zeros((leng, 3),dtype=object)
            anovas.shape = (leng,3)
            for c in sorted(conditions):
                df1 = df_swarm4.xs((c), axis =1, drop_level=False)
                df2 = df1.melt(ignore_index=False)
                df2 = df2.reset_index()
                df2.columns = ["Strain","Condition","Replicate","Score"]
                res = stat()
                res.tukey_hsd(df=df2, res_var='Score', xfac_var='Strain', anova_model='Score ~ C(Strain)')
                #takes only those compared to the WT
                stats1 = res.tukey_summary[res.tukey_summary['group1'] == ("0_"+WT)]
                # makes list of the current condition same length as number of genes - 1. 
                # Allows for zipping to the pairs and p-values for the matrix.
                cond = ([c] * leng)
                pairs = list(zip(stats1.group2, stats1.group2))
                pvalues = list(stats1['p-value'])
                grouped = list(zip(pairs,cond,pvalues))
                grouped = np.array(grouped,dtype=object)
                anovas = np.append(anovas,grouped, axis=0)
            anovas = anovas[leng:len(anovas)]
            df_swarm4 = df_swarm4[df_swarm4.index != ("0_"+WT)]
            #transposes dataset such that can be split by strain instead of condition.
            b = df_swarm4.T
            conditions = {x for x in b.columns}
            for c,i in zip(sorted(conditions),range(len(conditions))):
                df1 = b.xs((c), axis=1, drop_level=False)
                df1 = df1.T
                df2 = df1.melt(ignore_index=False)
                df2.columns = ["Condition","Replicate","Score"]
                #finds the annotation pair (condition,condition) and pvalue for the current strain for annotation of significance scores.
                m = [row for row in anovas if c == row[0][0]]
                pair1 = [(row[1],row[1]) for row in m]
                pvalue1 = [row[2] for row in m]
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
                # adds p-value annotations for each condition.
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
                plt.title(c)
                plt.xticks(rotation=rote)
                fig.set_size_inches(width1, height1)
                fig.set_dpi(1000)
                sns.despine()
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                              fancybox=True, shadow=False, ncol=1)
                plt.savefig(outputfile1+"_"+c+"_"+IRIS+"_Swarm_Plot.pdf", bbox_inches='tight')
                plt.clf()
    #produces heatmap based on the same dataset used for the barplots
    # transposes data such that rows are conditions
    b = df3.T
    # averages replicate scores
    c = b.groupby(level=0).mean()
    #transposes back for plotting of the heatmap
    x = c.T
    for i in range(len(x.columns)):
        x = x.rename(columns={x.columns[i]: x.columns[i]})
    c = x.T
    sns.set(style="white")
    sns.set_context("paper")
    colp = sns.color_palette(colour_heatmap, as_cmap=True)
    fig, ax = plt.subplots()
    ax = sns.heatmap(c,center=1,cmap=colp,square=False,annot=True,annot_kws={"size": 4},fmt=".4f",linewidths=.1,linecolor='0',vmin=0, vmax=2)
    fig.set_size_inches(width1, height1)
    plt.savefig(outputfile1+"_"+IRIS+"_Heatmap.pdf", bbox_inches='tight')
    plt.clf()
    return df3

if __name__ == "__main__":
    main()