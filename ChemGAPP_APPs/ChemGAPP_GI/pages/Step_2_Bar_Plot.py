
from re import M
import streamlit as st
import pandas as pd
import numpy as np
import os
import io
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from io import BytesIO, StringIO
from bioinfokit.analys import stat
from statannotations.Annotator import Annotator

st.set_page_config(layout="wide")
st.title('ChemGAPP GI: Bar Plot')
#path1 = st.sidebar.text_input("Path to the folder where your interaction score files are stored:")
uploaded_files = st.sidebar.file_uploader("Upload X_Interaction_Scores.csv files", accept_multiple_files=True)
primary_gene1 = st.sidebar.text_input("Name of the primary gene common to all files:")
color1 = st.sidebar.color_picker('Click to select Primary Gene colour', '#F19D7B')
color2 = st.sidebar.color_picker('Click to select Secondary Gene colour', '#FBE3CA')
color3 = st.sidebar.color_picker('Click to select Double Expected colour', '#9A9945')
color4 = st.sidebar.color_picker('Click to select Double Observed colour', '#6EBFBA')
my_expandera = st.expander(label="Instructions:", expanded=True)
my_expandera.markdown("""
1. Upload the desired interaction score files. 
```
These are the X_Interaction_Scores.csv files from the previous step.
Ensure the file names end with 'Interaction_Scores.csv'. 

E.g.

📂Interaction_Score_Files
 ┣ 📜 A_Interaction_Scores.csv
 ┣ 📜 B_Interaction_Scores.csv
 ┗ 📜 C_Interaction_Scores.csv
 
```

-----

2. Input the name of the Primary Gene. This is the gene common to all the tested gene pairs.
    
```
E.g.

MexA::MexB
MexA::OmpM
MexA::MexY

MexA would be the Primary Gene

```

-----

3. Click on the coloured squares to select the colour of the bars within the bar plot.

-----

4. Press 'Begin!'

-----

5. Click 'Download Image'. This will store the file within your downloads folder. 

-----

""")
if "plot" not in st.session_state:
        st.session_state.plot = []
if st.session_state.plot:
    col1, col2,col3 = st.columns((1,3,1))
    element = col2.pyplot(st.session_state.plot)
complete = st.sidebar.button(label="Begin!")
if complete:
    #if already run once, will rerun again after you click download.
    if st.session_state.plot:
        element.empty()
        def GI_BarPlot(primary_gene):
            m_s = None 
            #adds the name of the secondary gene into a column.         
            for a in uploaded_files:
                if m_s is None:
                    m_s = pd.read_csv(a,index_col=[0])
                    m_s['Gene'] = a.name.split('/')[-1].split("_Interaction_Scores.csv")[0]
                else:
                    m_s2 = pd.read_csv(a,index_col=[0])
                    m_s2['Gene'] = a.name.split('/')[-1].split("_Interaction_Scores.csv")[0]
                    m_s = pd.concat([m_s,m_s2], axis=0)
            
            secondary_genes = set(m_s.Gene.values)
            # performs ANOVA and Tukey HSD between the expected and experimental double knockout fitness ratios
            p_value_list = []
            pairs_list = []
            for i in sorted(secondary_genes):
                df1 = m_s[m_s['Gene'] == i]
                df1_doub = df1[["Double Expected","Double Observed"]]
                df1_doub = df1_doub.melt()
                res = stat()
                res.tukey_hsd(df=df1_doub, res_var='value', xfac_var='variable', anova_model='value ~ C(variable)')
                stats1 = res.tukey_summary
                pairs = list(zip(stats1.group1, stats1.group2))
                pairs_list.append(pairs)
                pvalues = list(stats1['p-value'])
                p_value_list.append(pvalues)
            #plots the bar plots, split by the secondary gene
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
                change_width(ax, 0.975)
                #sets colours of bars based on colour values picked.
                for p,i in zip(ax.patches,range(len(ax.patches))):
                        if i == 0:
                            p.set_color(color1)
                        if i == 1:
                            p.set_color(color2)
                        if i == 2:
                            p.set_color(color3)
                        if i == 3:
                            p.set_color(color4)
            #annotates the anova p-values for each secondary gene plot.
            for ax,pairs,pvalues in zip(g.axes,pairs_list,p_value_list):
                boxax = sns.boxplot(showmeans=False,
                        meanline=False,
                        medianprops={'visible': False},
                        whiskerprops={'visible': False},whis=0,
                        data=m_s,
                        showfliers=False,
                        showbox=False,
                        showcaps=False,
                        ax=ax)
                plt.xticks([])
                annotator = Annotator(boxax, pairs, data=m_s)
                annotator.set_pvalues(pvalues)
                annotator.configure(line_width = 1)
                annotator.annotate()
            g.set_titles("{col_name}",size=15,y=-0.08)
            #adds legend and matches colours to slider colours selected.
            eight = mlines.Line2D([], [], color=color1, marker='s', ls='', label=primary_gene,ms=15)
            nine = mlines.Line2D([], [], color=color2, marker='s', ls='', label="Secondary Gene",ms=15)
            ten = mlines.Line2D([], [], color=color3, marker='s', ls='', label="Double Expected",ms=15)
            eleven = mlines.Line2D([], [], color=color4, marker='s', ls='', label="Double Observed",ms=15)
            g.axes.flat[0].set_ylabel("Fitness Ratio", fontsize=15)
            plt.legend(handles=[eight,nine,ten,eleven],
                       loc='center left', bbox_to_anchor=(1, 0.5),
                  fancybox=True, shadow=False, ncol=1,fontsize=14)
            plt.subplots_adjust(hspace=0.4, wspace=0.05)
            fn = (primary_gene1+"_GI_Plot.pdf")
            img = io.BytesIO()
            plt.savefig(img, format='pdf', bbox_inches='tight')
            btn = st.download_button(
               label="Download image",
               data=img,
               file_name=fn,
               mime="image/pdf"
            )
            col1, col2,col3 = st.columns((1,2,1))
            col2.pyplot(g)
            return g
        g = GI_BarPlot(primary_gene1)
        if "plot" not in st.session_state:
                st.session_state.plot = g
        elif "plot" in st.session_state:
                st.session_state.plot = g
    else:
        #repeats above code such that plots always load after downloading plot.
        def GI_BarPlot(primary_gene):
            m_s = None          
            for a in uploaded_files:
                if m_s is None:
                    m_s = pd.read_csv(a,index_col=[0])
                    m_s['Gene'] = a.name.split('/')[-1].split("_Interaction_Scores.csv")[0]
                else:
                    m_s2 = pd.read_csv(a,index_col=[0])
                    m_s2['Gene'] = a.name.split('/')[-1].split("_Interaction_Scores.csv")[0]
                    m_s = pd.concat([m_s,m_s2], axis=0)
            secondary_genes = set(m_s.Gene.values)
            p_value_list = []
            pairs_list = []
            for i in sorted(secondary_genes):
                df1 = m_s[m_s['Gene'] == i]
                df1_doub = df1[["Double Expected","Double Observed"]]
                df1_doub = df1_doub.melt()
                res = stat()
                res.tukey_hsd(df=df1_doub, res_var='value', xfac_var='variable', anova_model='value ~ C(variable)')
                stats1 = res.tukey_summary
                pairs = list(zip(stats1.group1, stats1.group2))
                pairs_list.append(pairs)
                pvalues = list(stats1['p-value'])
                p_value_list.append(pvalues)
            
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
            
            for ax,pairs,pvalues in zip(g.axes,pairs_list,p_value_list):
                boxax = sns.boxplot(showmeans=False,
                        meanline=False,
                        medianprops={'visible': False},
                        whiskerprops={'visible': False},whis=0,
                        data=m_s,
                        showfliers=False,
                        showbox=False,
                        showcaps=False,
                        ax=ax)
                plt.xticks([])
                annotator = Annotator(boxax, pairs, data=m_s)
                annotator.set_pvalues(pvalues)
                annotator.configure(line_width = 1)
                annotator.annotate()
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
            fn = (primary_gene1+"_GI_Plot.pdf")
            img = io.BytesIO()
            plt.savefig(img, format='pdf', bbox_inches='tight')
            btn = st.download_button(
               label="Download image",
               data=img,
               file_name=fn,
               mime="image/pdf"
            )
            col1, col2,col3 = st.columns((1,2,1))
            col2.pyplot(g)
            return g
        g = GI_BarPlot(primary_gene1)
        if "plot" not in st.session_state:
                st.session_state.plot = g
        elif "plot" in st.session_state:
                st.session_state.plot = g