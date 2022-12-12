#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
parser = argparse.ArgumentParser(description="Produces barplot displaying genetic interaction type.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-p", "--PATH", help="Path to Interaction Score Files")
parser.add_argument("-o", "--OutputFile", help="A PDF file of the final bar plot.")
parser.add_argument("-g", "--PrimaryGene", help="The primary interacting gene being compared")
args = vars(parser.parse_args())
outputfile1 = args["OutputFile"]
path1 = args["PATH"]
primary_gene1 = args["PrimaryGene"]
def GI_BarPlot(Path,outputfile,primary_gene):
    directory = Path
    files = []
    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        # checking if it is a file
        if f.endswith("Interaction_Scores.csv"):
            if os.path.isfile(f):
                print((f.split('/')[-1].split("_Interaction_Scores.csv")[0]))
                files.append(f)
    m_s = None          
    for a in files:
        if m_s is None:
            m_s = pd.read_csv(a,index_col=[0])
            m_s['Gene'] = a.split('/')[-1].split("_Interaction_Scores.csv")[0]
        else:
            m_s2 = pd.read_csv(a,index_col=[0])
            m_s2['Gene'] = a.split('/')[-1].split("_Interaction_Scores.csv")[0]
            m_s = pd.concat([m_s,m_s2], axis=0)
    m_s
    sns.set_theme(style="white")
    g = sns.catplot(data=m_s,kind="bar", alpha=0.8, height=6,col="Gene",col_wrap=15,aspect=.4,capsize=.15,errwidth=0.5)
    g.set_xticklabels([])
    def change_width(ax, new_value) :
        for patch in ax.patches :
            current_width = patch.get_width()
            diff = current_width - new_value

            # we change the bar width
            patch.set_width(new_value)

            # we recenter the bar
            patch.set_x(patch.get_x() + diff * .5)

    for ax in g.axes:
        change_width(ax, 1)
        for p,i in zip(ax.patches,range(len(ax.patches))):
                if i == 0:
                    p.set_color('coral')
                if i == 1:
                    p.set_color('peachpuff')
                if i == 2:
                    p.set_color('olive')
                if i == 3:
                    p.set_color('lightseagreen')
    g.set_titles("{col_name}",size=15,y=-0.08)
    eight = mlines.Line2D([], [], color='coral', marker='s', ls='', label=primary_gene,ms=15)
    nine = mlines.Line2D([], [], color='peachpuff', marker='s', ls='', label="Secondary Gene",ms=15)
    ten = mlines.Line2D([], [], color='olive', marker='s', ls='', label="Double Expected",ms=15)
    eleven = mlines.Line2D([], [], color='lightseagreen', marker='s', ls='', label="Double Observed",ms=15)
    g.axes.flat[0].set_ylabel("Fitness Ratio", fontsize=15)
    plt.legend(handles=[eight,nine,ten,eleven],
               loc='center left', bbox_to_anchor=(1, 0.5),
          fancybox=True, shadow=False, ncol=1,fontsize=14)
    plt.subplots_adjust(hspace=0.4, wspace=0.05)
    #g.fig.suptitle((primary_gene+" Genetic Interactions"),size=15)
    plt.savefig(outputfile, bbox_inches='tight')
GI_BarPlot(path1,outputfile1,primary_gene1)