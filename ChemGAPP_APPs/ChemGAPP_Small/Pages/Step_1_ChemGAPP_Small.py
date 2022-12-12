import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import scipy.stats
from scipy.stats import wilcoxon
from scipy.stats import ranksums
import re
import io 
from matplotlib.pyplot import figure
from bioinfokit.analys import stat
from statannotations.Annotator import Annotator
from tempfile import NamedTemporaryFile
import pickle


st.set_page_config(layout="wide")
st.title('ChemGAPP Small')
my_expandera = st.expander(label="Instructions:", expanded=True)
my_expandera.markdown("""

1- First upload all iris files you wish to include.
Ensure IRIS file names are in the format: `CONDITION-concentration-platenumber_replicate.JPG.iris`

E.g. `AMPICILLIN-50mM-6_B.JPG.iris`

Where concentrations have decimals, use a comma instead of a period:

E.g. `AMPICILLIN-0,5mM-6_B.JPG.iris`

-----

2- Enter a path to the folder you would like to save the output files to. Ensure you include a prefix which will be added to the start of all output file names e.g: ~/Desktop/project/Test

Output files include:

`Test_Intial_dataset.csv`

`Test_Normalised_dataset.csv`

`Test_Scored_Dataset.csv`

`Test_Final_dataset.csv`

-----

3- Upload your plate information files.

These should be txt files with the format:

| Row | Column | Strain |
|  ------------- | :-------------: | :-------------: |
| 1 | 1 | WT |
| 1 | 2 | WT  |
| 1 | 3 | Mutant1 |
| 1 | 4 | Mutant1  |
| ... |...| ... |
| 16 | 23 | Mutantx |
| 16 | 24 | Mutantx |

-----

4- Enter the name of the wild type strain. 

This must match the name given in the plate information file.
E.g. `WT`

-----

5- Optionally enter the name of a Seaborn colour palette to change the colour of the output plots. 

The default palette is `Spectral`.

-----

6- Select how you want the bar plots to be group. Either by strain or by condition.

-----

7- Press Begin!

-----

8- To save bar plot images, click on the `Download image` button beneath the plot. This will save the image as a pdf file.

    * Files will be downloaded to your Downloads Folder

    * Do not try save an image until all plots are produced.

""")
def new_click():
    if 'figures' not in st.session_state:
        pass
    else:
        del st.session_state["figures"]
        
def del_pickle(fnm):
    if st.session_state.figures:
        del st.session_state["figures"]
    dir = os.path.abspath(os.getcwd())
    for fa in fnm:
        path = os.path.join(dir, fa)
        os.remove(path)
par_dir = os.path.abspath(os.getcwd())
new_dir = "Temp"
temp_path = os.path.join(par_dir, new_dir)
try:
    os.mkdir(temp_path)
except:
    pass

if 'figures' not in st.session_state:
        st.session_state.figures = []
        for f in sorted(os.listdir(temp_path)):
            if f.endswith(".pdf"): 
                path = os.path.join(temp_path, f)
                os.remove(path)
if st.session_state.figures:
    my_expandera = st.expander(label="Scored Dataset", expanded=False)
    my_expandera.write(st.session_state.scored_dataset, header=[0,1,2])
    elwarn2 = st.warning("Preparing plots for download, do not download until ready.")
    for fig,i,c,fna in zip(st.session_state.figures, range(len(st.session_state.figures)),st.session_state.conds,st.session_state.fnames):
        if i == 0:
            myexp = st.expander(label="Heatmap:", expanded=False)
            col1a,col1b,col1c = myexp.columns((1,4,1))
            col1b.pyplot(fig)
            pathfna = os.path.join("Temp", fna)
            with open(pathfna,'rb') as fid:
                ax = pickle.load(fid)
                imga = io.BytesIO()
            plt.savefig(imga, format='pdf', bbox_inches='tight')
            fn =("ChemGAPP_Heatmap.pdf")
            myexp.download_button(
                            label="Download image",
                            data=imga,
                            file_name=fn,
                            mime="image/pdf"
                         )
            st.write(str(st.session_state.plot_type)+"s")   
        else:
            if st.session_state.width_type < 10:
                a = i
                if (a % 2) != 0:
                    cola, colb = st.columns((1,1))
                    myexpand = cola.expander(label=c.replace(",","."), expanded=False)
                    myexpand.pyplot(fig)
                    pathfna = os.path.join("Temp", fna)
                    with open(pathfna,'rb') as fid:
                        ax = pickle.load(fid)
                        imga = io.BytesIO()
                    fn =("ChemGAPP_"+c+"_"+str(st.session_state.iris_type)+"_"+str(st.session_state.plot_type)+".pdf")
                    plt.savefig(imga, format='pdf', bbox_inches='tight')
                    myexpand.download_button(
                        label="Download image",
                        data=imga,
                        file_name=fn,
                        mime="image/pdf")
                else:
                    myexpand = colb.expander(label=c.replace(",","."), expanded=False)
                    myexpand.pyplot(fig)
                    pathfna = os.path.join("Temp", fna)
                    with open(pathfna,'rb') as fid:
                        ax = pickle.load(fid)
                        imga = io.BytesIO()
                    fn = ("ChemGAPP_"+c+"_"+str(st.session_state.iris_type)+"_"+str(st.session_state.plot_type)+".pdf")
                    plt.savefig(imga, format='pdf', bbox_inches='tight')
                    myexpand.download_button(
                         label="Download image",
                         data=imga,
                         file_name=fn,
                         mime="image/pdf"
                      )
            if st.session_state.width_type >= 10:
                myexpand = st.expander(label=c.replace(",","."), expanded=False)
                myexpand.pyplot(fig)
                pathfna = os.path.join("Temp", fna)
                with open(pathfna,'rb') as fid:
                    ax = pickle.load(fid)
                    imga = io.BytesIO()
                imga = io.BytesIO()
                fn =("ChemGAPP_"+c+"_"+str(st.session_state.iris_type)+"_"+str(st.session_state.plot_type)+".pdf")
                plt.savefig(imga, format='pdf', bbox_inches='tight')
                myexpand.download_button(
                    label="Download image",
                    data=imga,
                    file_name=fn,
                    mime="image/pdf")
    elwarn2.empty()
    st.success("Images ready to download!")
if "iris_type" not in st.session_state:
        st.session_state.iris_type = []
def new_iris_changed():
        if st.session_state.new_iris_type:
            st.session_state.iris_type = st.session_state.new_iris_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "plot_type" not in st.session_state:
        st.session_state.plot_type = []
def new_plot_changed():
        if st.session_state.new_plot_type:
            st.session_state.plot_type = st.session_state.new_plot_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "group_type" not in st.session_state:
        st.session_state.group_type = []
def new_group_changed():
        if st.session_state.new_group_type:
            st.session_state.group_type = st.session_state.new_group_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "width_type" not in st.session_state:
        st.session_state.width_type = []
def new_width_changed():
        if st.session_state.new_width_type:
            st.session_state.width_type = st.session_state.new_width_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "height_type" not in st.session_state:
        st.session_state.height_type = []
def new_height_changed():
        if st.session_state.new_height_type:
            st.session_state.height_type = st.session_state.new_height_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "colour_type" not in st.session_state:
        st.session_state.colour_type = []
def new_colour_changed():
        if st.session_state.new_colour_type:
            st.session_state.colour_type = st.session_state.new_colour_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "in_type" not in st.session_state:
        st.session_state.in_type = []
def new_in_changed():
        if st.session_state.new_in_type:
            st.session_state.in_type = st.session_state.new_in_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "out_type" not in st.session_state:
        st.session_state.out_type = []
def new_out_changed():
        if st.session_state.new_out_type:
            st.session_state.out_type = st.session_state.new_out_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]


if "info_type" not in st.session_state:
        st.session_state.info_type = []
def new_info_changed():
        if st.session_state.new_info_type:
            st.session_state.info_type = st.session_state.new_info_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "WT_type" not in st.session_state:
        st.session_state.WT_type = []
def new_WT_changed():
        if st.session_state.new_WT_type:
            st.session_state.WT_type = st.session_state.new_WT_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "rote_type" not in st.session_state:
        st.session_state.rote_type = []
def new_rote_changed():
        if st.session_state.new_rote_type:
            st.session_state.rote_type = st.session_state.new_rote_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]


if "circ_type" not in st.session_state:
        st.session_state.circ_type = []
def new_circ_changed():
        if st.session_state.new_circ_type:
            st.session_state.circ_type = st.session_state.new_circ_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "xlab_type" not in st.session_state:
        st.session_state.xlab_type = []
def new_xlab_changed():
        if st.session_state.new_xlab_type:
            st.session_state.xlab_type = st.session_state.new_xlab_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "ylim_type" not in st.session_state:
        st.session_state.ylim_type = []
def new_ylim_changed():
        if st.session_state.new_ylim_type:
            st.session_state.ylim_type = st.session_state.new_ylim_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "ylow_type" not in st.session_state:
        st.session_state.ylow_type = []
def new_ylow_changed():
        if st.session_state.new_ylow_type:
            st.session_state.ylow_type = st.session_state.new_ylow_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "yhigh_type" not in st.session_state:
        st.session_state.yhigh_type = []
def new_yhigh_changed():
        if st.session_state.new_yhigh_type:
            st.session_state.yhigh_type = st.session_state.new_yhigh_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "pal_type" not in st.session_state:
        st.session_state.pal_type = "icefire"
def new_pal_changed():
        if st.session_state.new_pal_type:
            st.session_state.pal_type = st.session_state.new_pal_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "pal2_type" not in st.session_state:
        st.session_state.pal2_type = "bwr_r"
def new_pal2_changed():
        if st.session_state.new_pal2_type:
            st.session_state.pal2_type = st.session_state.new_pal2_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "remove_type" not in st.session_state:
        st.session_state.remove_type = []
def new_remove_changed():
        if st.session_state.new_remove_type:
            st.session_state.remove_type = st.session_state.new_remove_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

if "order_type" not in st.session_state:
        st.session_state.order_type = []
def new_order_changed():
        if st.session_state.new_order_type:
            st.session_state.order_type = st.session_state.new_order_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]

uploaded_files = st.sidebar.file_uploader("Upload multiple IRIS files", accept_multiple_files=True)
plate_info_files = st.sidebar.file_uploader("Upload multiple plate information files", accept_multiple_files=True)
outputfile1 = st.sidebar.text_input("Enter Path and Prefix for Output Files:", on_change=new_out_changed, key="new_out_type")
wildtype1 = st.sidebar.text_input("Name of the Wildtype Strain:", on_change=new_WT_changed, key="new_WT_type")
iris_type = st.sidebar.radio(
    "Select what the IRIS morphology to test:",
    ('size', 'circularity','opacity'), on_change=new_iris_changed, key="new_iris_type")



if iris_type == 'size':
    max_size = st.sidebar.text_input('Optional: Max value for colony '+iris_type)
    if max_size:
        max_size = max_size
    else:
        max_size = None
if iris_type == 'circularity' or iris_type == 'opacity':
    max_size = None

plot_type = st.sidebar.radio(
    "Select what the plots are grouped by:",
    ('Strain', 'Condition'), on_change=new_group_changed, key="new_group_type")
bar_swarm = st.sidebar.radio(
    "Select type of plot:",
    ('Bar Plot', 'Swarm Plot'), on_change=new_plot_changed, key="new_plot_type")

st.sidebar.title("Customisation:")
remove_y_n = st.sidebar.radio(
    "Would you like to remove strains?",
    ('No','Yes'), on_change=new_remove_changed, key="new_remove_type")
if remove_y_n == 'Yes':
    txt = st.sidebar.text_area('Write strains you wish to remove, seperate by \';\'') 
order_type = st.sidebar.radio(
    "Would you like to order your X-axis labels?",
    ('No','Yes'), on_change=new_order_changed, key="new_order_type")
if order_type == 'No':
    order_list = None
if order_type == 'Yes':
    order_set = st.sidebar.text_area('Write order of labels, seperate by \';\'')
colour_palette2 = st.sidebar.selectbox('Colour Palette for Heatmap:',('bwr_r','Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 
 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 
 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 
 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 
 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 
 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 
 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 
 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 
 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 
 'bone_r', 'brg', 'brg_r', 'bwr', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r',
 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r',
 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 
 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r',
 'hot', 'hot_r', 'hsv', 'hsv_r', 'icefire', 'icefire_r', 'inferno', 
 'inferno_r', 'magma', 'magma_r', 'mako', 'mako_r', 
 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r',
 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r',
 'rocket', 'rocket_r', 'seismic', 'seismic_r', 'spring', 'spring_r',
 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b',
 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'twilight',
 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'vlag', 'vlag_r', 'winter', 'winter_r'), on_change=new_pal2_changed, key="new_pal2_type")
colour_palette = st.sidebar.selectbox('Colour Palette for Bar/Swarm Plots:',('icefire','Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 
 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 
 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 
 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 
 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 
 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 
 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 
 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 
 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 
 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r',
 'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r',
 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 
 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r',
 'hot', 'hot_r', 'hsv', 'hsv_r', 'icefire_r', 'inferno', 
 'inferno_r', 'magma', 'magma_r', 'mako', 'mako_r', 
 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r',
 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r',
 'rocket', 'rocket_r', 'seismic', 'seismic_r', 'spring', 'spring_r',
 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b',
 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'twilight',
 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'vlag', 'vlag_r', 'winter', 'winter_r'), on_change=new_pal_changed, key="new_pal_type")

width1 = st.sidebar.slider('Figure Width', 1, 20, 5, on_change=new_width_changed, key="new_width_type")
height1 = st.sidebar.slider('Figure Height', 1, 20, 5,on_change=new_height_changed, key="new_height_type")
if bar_swarm == 'Bar Plot':
    rote = st.sidebar.slider('X Label Rotation', 0, 180, 90,on_change=new_rote_changed, key="new_rote_type")
if bar_swarm == 'Swarm Plot':
    circ_size = st.sidebar.slider('Swarm Plot Circle Size', 1.0, 10.0, 3.0,on_change=new_circ_changed, key="new_circ_type")
    xlabes = st.sidebar.radio(
    "Include X Label Names:",
    ('No', 'Yes'), on_change=new_xlab_changed, key="new_xlab_type")
    if xlabes == 'Yes':
        rote = st.sidebar.slider('X Label Rotation', 0, 180, 90,on_change=new_rote_changed, key="new_rote_type")
sety = st.sidebar.radio(
    "Adjust y axis limits:",
    ('No', 'Yes'), on_change=new_ylim_changed, key="new_ylim_type")
if sety == 'Yes':
    ylow = st.sidebar.slider('Y Axis Lower Limit', 0.0, 10.0, 0.0, on_change=new_ylow_changed, key="new_ylow_type")
    yhigh = st.sidebar.slider('Y Axis Higher Limit', 0.0, 10.0, 2.0, on_change=new_yhigh_changed, key="new_yhigh_type")

complete = st.sidebar.button(label="Begin!",on_click=new_click)

if complete:
    del st.session_state["figures"]
    if remove_y_n == 'Yes':
        inclstrains = txt.split(";")
    if order_type == 'Yes':
        order_list = order_set.split(";")
    if st.session_state.iris_type == []:
        st.session_state.iris_type = iris_type
    if st.session_state.plot_type == []:
        st.session_state.plot_type = bar_swarm
    m = None
    for f in uploaded_files:
        g = pd.read_csv(f,
                comment='#',
                index_col=[0, 1],
                sep='\t')
        if m is None:
            try:
                m = g[iris_type]
            except:
                m = g[iris_type]
            m.name = (int(f.name.split('-')[2].split('_')[0]),
                      '-'.join(f.name.split('.')[0].split('-')[:2]).replace(",",".").replace("-"," "),
                      f.name.split('.')[0].split('_')[1])
            m = m.to_frame()
        else:
            try:
                m1 = g[iris_type]
            except:
                m1 = g[iris_type]
            m1.name = (int(f.name.split('-')[2].split('_')[0]),
                      '-'.join(f.name.split('.')[0].split('-')[:2]).replace(",",".").replace("-"," "),
                      f.name.split('.')[0].split('_')[1])
            m1 = m1.to_frame()
            #sets them in a dataframe grouped by the condition and then columns as the A B C D or E etc
            m = m.join(m1, how='inner')
    m.to_csv((outputfile1+"_Initial_dataset.csv"))
    m = m.apply(pd.to_numeric)
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
    for c1,c2,i in zip(sorted(df_pmm.columns),sorted(df_outer.columns), range(len(m.columns))):
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
                elif ind[0] != 1 or ind[0] != rlen or ind[1] != 1 or ind[1] != clen or ind[0] != 2 or ind[0] != (rlen-1) or ind[1] != 2 or ind[1] != (clen-1):
                    ar2[j] = ((ar2[j])*(m_medain/PMM))
            ar2.shape = (mlen,1)
            ardf = np.concatenate((ardf,ar2), axis=1)
        else:
            for ind, j in zip(m.index,range(len(ar2))):
                ar2[j] = ((ar2[j])*(m_medain/PMM))
            ar2.shape = (mlen,1)
            ardf = np.concatenate((ardf,ar2), axis=1)
    ardf = pd.DataFrame(ardf, index=m.index)
    ardf = ardf.iloc[: , 1:]
    ardf.columns = (pd.MultiIndex.from_tuples(sorted(m.columns)))
    ardf = ardf[sorted(ardf)]
    ardf.to_csv(outputfile1+"_Normalised_dataset.csv")
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
        return [ atoi(c) for c in re.split(r'(\d+)', text) ]

    plate_DFs = []
    for i in plate_info_files:
        p = plate_info(i)
        plate_DFs.append(p)
    plates = {x[0] for x in ardf.columns}
    nm4 = ardf.copy(deep=False)
    columns1 = {x[1:4] for x in nm4.columns}
    df_with_strains = pd.DataFrame()
    for a, n in zip(plate_DFs, sorted(plates)):
        df1 = (ardf.xs((n), axis=1, drop_level=True))
        df1.columns = [f'{i}_{j}' for i,j in df1.columns]
        df2 = pd.merge(df1,a,left_index=True, right_index=True)
        df_with_strains = pd.concat([df_with_strains , df2], ignore_index=True)
    df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
    df_with_strains = df_with_strains.set_index('Gene')
    df_with_strains.columns = df_with_strains.columns.str.split('_', expand=True)
    df_with_strains = df_with_strains.sort_index()
    if max_size != None:
        df_with_strains[df_with_strains > float(max_size)] = np.nan
    
    df_with_strains= df_with_strains.reset_index()
    df_with_strains = df_with_strains[df_with_strains['Gene'] != "EMPTY"]
    if remove_y_n == 'Yes':
        for stra in inclstrains:
            df_with_strains = df_with_strains[df_with_strains['Gene'] != stra]
    df_with_strains = df_with_strains.set_index(['Gene'])

    df_with_strains.to_csv(outputfile1+"_Final_dataset.csv")
    group = df_with_strains.groupby(level=0).mean()
    group2 = group.copy(deep=False)
    df3 = pd.DataFrame(index=group.index)
    for column in group:
        df1 = group[column]
        df1 = pd.DataFrame(df1)
        ar1 = np.array(df1)
        ar2 = np.array(df1)
        wt_mean = float(df1.loc[wildtype1])
        for ind, j in zip(df1.index,range(len(ar1))):
            ar2[j] = (float(ar1[j])/wt_mean)
        df3 = np.concatenate((df3,ar2), axis=1)    
    df3 = pd.DataFrame(df3, index=group.index, columns=group.columns)
    df3.index.name = None
    df3 = df3.drop(wildtype1,axis=0)
    my_expandera = st.expander(label="Scored Dataset", expanded=False)
    my_expandera.write(df3, header=[0,1,2])
    df3.to_csv(outputfile1+"_Scored_Dataset.csv")
    if "scored_dataset" not in st.session_state:
        st.session_state.scored_dataset = df3
    elif "scored_dataset" in st.session_state:
        st.session_state.scored_dataset = df3
    
    elwarn = st.warning("Preparing plots for download, do not download until ready.")
    conditions = {x[0] for x in df3.columns}
    clen = len(conditions)
    slen = len(df3)
    figs = []
    imgs = []
    fns = []
    cs = []
    b = df3.T
    c = b.groupby(level=0).mean()
    x = c.T
    for i in range(len(x.columns)):
        x = x.rename(columns={x.columns[i]: x.columns[i].replace(",",".")})
    c = x.T
    #c = c.reset_index()
    #c = pd.melt(c, id_vars =['index'])
    #if order_type == "Yes":
    #    c['value'] = pd.Categorical(c['value'],dtype='category', categories=order_list, ordered=True)
    #c = c.pivot('index', 'variable', 'value')
    #st.write(c)
    sns.set(style="white")
    sns.set_context("paper")
    fig, ax = plt.subplots()
    colp = sns.color_palette(st.session_state.pal2_type, as_cmap=True)
    ax = sns.heatmap(c,center=1,cmap=colp,square=False,annot=True,annot_kws={"size": 4},fmt=".4f",linewidths=.1,linecolor='0',vmin=0, vmax=2)
    myexp = st.expander(label="Heatmap:", expanded=False)
    col1a,col1b,col1c = myexp.columns((1,4,1))
    col1b.pyplot(fig)
    figs.append(fig)
    c="_"
    cs.append(c)
    img = io.BytesIO()
    fn =("ChemGAPP_Heatmap.pdf")
    temp = os.path.join(new_dir,fn)
    with open(temp,'wb') as fid:
        pickle.dump(fig, fid)
    fns.append(fn)
    plt.savefig(img, format='pdf', bbox_inches='tight')
    myexp.download_button(
                         label="Download image",
                         data=img,
                         file_name=fn,
                         mime="image/pdf"
                      )
    
    st.write(bar_swarm+"s")
    if plot_type == 'Condition':
        if bar_swarm == 'Bar Plot':
            for c,i in zip(sorted(conditions),range(len(conditions))):
                df1 = df3.xs((c), axis =1, drop_level=False)
                df2 = df1.melt(ignore_index=False)
                df2= df2.reset_index()
                df2.columns = ["Strain","Condition","Replicate","Score"]
                sns.set(style="white")
                sns.set_context("paper")
                fig, ax = plt.subplots()
                ax = sns.barplot(x="Strain", y="Score",data=df2,capsize=.5,errwidth=0.8,palette=st.session_state.pal_type,edgecolor="0.2", order=order_list)
                ax = sns.swarmplot(x="Strain", y="Score",color="0", data=df2, alpha=1,size=4,order=order_list)
                ax.set(xlabel= "",ylabel='Fitness Ratio')
                if sety == "Yes":
                    ax.set_ylim(ylow,yhigh)
                fig.set_size_inches(width1, height1)
                fig.set_dpi(1000)
                plt.title(c.replace(",","."))
                plt.xticks(rotation=rote)
                sns.despine()
                a= i+1
                if (a % 2) != 0:
                    cola, colb = st.columns((1,1))
                    myexpand = cola.expander(label=c.replace(",","."), expanded=False)
                    myexpand.pyplot(fig)
                    figs.append(fig)
                    cs.append(c)
                    img = (io.BytesIO())
                    
                    plt.savefig(img, format='pdf', bbox_inches='tight')
                    fn =("ChemGAPP_"+c+"_"+iris_type+"_"+bar_swarm+".pdf")
                    temp = os.path.join(new_dir,fn)
                    with open(temp,'wb') as fid:
                        pickle.dump(fig, fid)
                    fns.append(fn)
                    but1 = myexpand.download_button(
                         label="Download image",
                         data=img,
                         file_name=fn,
                         mime="image/pdf"
                      )
                else:
                    myexpand = colb.expander(label=c.replace(",","."), expanded=False)
                    myexpand.pyplot(fig)
                    figs.append(fig)
                    cs.append(c)
                    img = (io.BytesIO())
                    
                    plt.savefig(img, format='pdf', bbox_inches='tight')
                    fn =("ChemGAPP_"+c+"_"+iris_type+"_"+bar_swarm+".pdf")
                    temp = os.path.join(new_dir,fn)
                    with open(temp,'wb') as fid:
                        pickle.dump(fig, fid)
                    fns.append(fn)
                    but2 = myexpand.download_button(
                         label="Download image",
                         data=img,
                         file_name=fn,
                         mime="image/pdf"
                      )
                
        if bar_swarm == 'Swarm Plot':
            df_swarm1 = df_with_strains.reset_index()
            df_swarm2 = df_swarm1[df_swarm1['Gene'] == wildtype1]
            df_swarm1['Gene'] = df_swarm1['Gene'].str.replace(wildtype1,("0_"+wildtype1))
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
                df2.columns = ["Strain","Condition","Replicate","Score"]
                res = stat()
                res.tukey_hsd(df=df2, res_var='Score', xfac_var='Strain', anova_model='Score ~ C(Strain)')
                stats1 = res.tukey_summary[res.tukey_summary['group1'] == ("0_"+wildtype1)]
                pairs = list(zip(stats1.group2, stats1.group2))
                pvalues = list(stats1['p-value'])
                df2 = df2[df2.Strain != ("0_"+wildtype1)]
                subcat_palette = sns.dark_palette("#8BF", reverse=True, n_colors=5)
                plotting_parameters = {
                    'data':    df2,
                    'x':       'Strain',
                    'y':       'Score'}
                sns.set(style="white")
                sns.set_context("paper")
                fig, ax1 = plt.subplots()
                ax1 = sns.swarmplot(x="Strain", y="Score",hue="Strain", data=df2, alpha=1,size=circ_size,palette=st.session_state.pal_type, order=order_list)
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
                        ax=ax1,order=order_list)
                annotator = Annotator(ax, pairs, **plotting_parameters,order=order_list)
                annotator.set_pvalues(pvalues)
                annotator.configure(line_width = 0)
                annotator.annotate()
                for patch in ax.artists:
                    r, g, b, a = patch.get_facecolor()
                    patch.set_facecolor((r, g, b, 0))    
                ax.set(xlabel= "Strain",ylabel='Fitness Ratio')
                if sety == "Yes":
                    ax1.set_ylim(ylow,yhigh)
                plt.title(c.replace(",",".").replace("-"," "))
                if xlabes == 'No':
                    plt.xticks([])
                if  xlabes == 'Yes':
                    plt.xticks(rotation=rote)
                fig.set_size_inches(width1, height1)
                
                sns.despine()
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                              fancybox=True, shadow=False, ncol=1)
                if width1 >= 10:
                    myexpand = st.expander(label=c, expanded=False)
                    myexpand.pyplot(fig)
                    figs.append(fig)
                    cs.append(c)
                    img =(io.BytesIO())
                    
                    plt.savefig(img, format='pdf', bbox_inches='tight')
                    fn = ("ChemGAPP_"+c+"_"+iris_type+"_"+bar_swarm+".pdf")
                    temp = os.path.join(new_dir,fn)
                    with open(temp,'wb') as fid:
                        pickle.dump(fig, fid)
                    fns.append(fn)
                    but3 = myexpand.download_button(
                         label="Download image",
                         data=img,
                         file_name=fn,
                         mime="image/pdf"
                      )
                if width1 < 10:
                    a= i+1
                    if (a % 2) != 0:
                        cola, colb = st.columns((1,1))
                        myexpand = cola.expander(label=c, expanded=False)
                        myexpand.pyplot(fig)
                        figs.append(fig)
                        cs.append(c)
                        img =(io.BytesIO())
                        
                        
                        plt.savefig(img, format='pdf', bbox_inches='tight')
                        fn = ("ChemGAPP_"+c+"_"+iris_type+"_"+bar_swarm+".pdf")
                        temp = os.path.join(new_dir,fn)
                        with open(temp,'wb') as fid:
                            pickle.dump(fig, fid)
                        fns.append(fn)
                        but4 = myexpand.download_button(
                         label="Download image",
                         data=img,
                         file_name=fn,
                         mime="image/pdf"
                      )
                    else:
                        myexpand = colb.expander(label=c, expanded=False)
                        myexpand.pyplot(fig)
                        figs.append(fig)
                        cs.append(c)
                        img =(io.BytesIO())
                        
                        
                        plt.savefig(img, format='pdf', bbox_inches='tight')
                        fn = ("ChemGAPP_"+c+"_"+iris_type+"_"+bar_swarm+".pdf")
                        temp = os.path.join(new_dir,fn)
                        with open(temp,'wb') as fid:
                            pickle.dump(fig, fid)
                        fns.append(fn)
                        but5 = myexpand.download_button(
                         label="Download image",
                         data=img,
                         file_name=fn,
                         mime="image/pdf"
                      )
                
    if plot_type == 'Strain':
        if bar_swarm == 'Bar Plot':
            b = df3.T
            conditions = {x for x in b.columns}
            for c,i in zip(sorted(conditions),range(len(conditions))):
                df1 = b.xs((c), axis=1, drop_level=False)
                df1 = df1.T
                df2 = df1.melt(ignore_index=False)
                df2.columns = ["Condition","Replicate","Score"]
                sns.set(style="white")
                sns.set_context("paper")
                fig, ax = plt.subplots()
                ax = sns.barplot(x="Condition", y="Score",data=df2,capsize=.5,errwidth=0.8,palette=st.session_state.pal_type,edgecolor="0.2", order=order_list)
                ax = sns.swarmplot(x="Condition", y="Score",color="0", data=df2, alpha=1,size=4,order=order_list)
                ax.set(xlabel= "",ylabel='Fitness Ratio')
                if sety == "Yes":
                    ax.set_ylim(ylow,yhigh)
                plt.xticks(rotation=rote)
                plt.title(c)
                fig.set_size_inches(width1, height1)
                
                sns.despine()
                a= i+1
                if (a % 2) != 0:
                    cola, colb = st.columns((1,1))
                    myexpand = cola.expander(label=c, expanded=False)
                    myexpand.pyplot(fig)
                    figs.append(fig)
                    cs.append(c)
                    img =(io.BytesIO())
                    
                    
                    plt.savefig(img, format='pdf', bbox_inches='tight')
                    fn =("ChemGAPP_"+c+"_"+iris_type+"_"+bar_swarm+".pdf")
                    temp = os.path.join(new_dir,fn)
                    with open(temp,'wb') as fid:
                        pickle.dump(fig, fid)
                    fns.append(fn)
                    but6 = myexpand.download_button(
                         label="Download image",
                         data=img,
                         file_name=fn,
                         mime="image/pdf"
                      )
                else:
                    myexpand = colb.expander(label=c, expanded=False)
                    myexpand.pyplot(fig)
                    figs.append(fig)
                    cs.append(c)
                    img =(io.BytesIO())
                    
                    
                    plt.savefig(img, format='pdf', bbox_inches='tight')
                    fn =("ChemGAPP_"+c+"_"+iris_type+"_"+bar_swarm+".pdf")
                    temp = os.path.join(new_dir,fn)
                    with open(temp,'wb') as fid:
                        pickle.dump(fig, fid)
                    fns.append(fn)
                    but7 = myexpand.download_button(
                         label="Download image",
                         data=img,
                         file_name=fn,
                         mime="image/pdf"
                      )
                
        if bar_swarm == "Swarm Plot":
            df_swarm1 = df_with_strains.reset_index()
            df_swarm2 = df_swarm1[df_swarm1['Gene'] == wildtype1]
            df_swarm1['Gene'] = df_swarm1['Gene'].str.replace(wildtype1,("0_"+wildtype1))
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
                df2.columns = ["Strain","Condition","Replicate","Score"]
                res = stat()
                res.tukey_hsd(df=df2, res_var='Score', xfac_var='Strain', anova_model='Score ~ C(Strain)')
                stats1 = res.tukey_summary[res.tukey_summary['group1'] == ("0_"+wildtype1)]
                print(stats1)
                cond = ([c] * leng)
                pairs = list(zip(stats1.group2, stats1.group2))
                pvalues = list(stats1['p-value'])
                grouped = list(zip(pairs,cond,pvalues))
                print(grouped)
                grouped = np.array(grouped,dtype=object)
                anovas = np.append(anovas,grouped, axis=0)
            print(anovas)
            anovas = anovas[leng:len(anovas)]
            df_swarm4 = df_swarm4[df_swarm4.index != ("0_"+wildtype1)]
            b = df_swarm4.T
            conditions = {x for x in b.columns}
            for c,i in zip(sorted(conditions),range(len(conditions))):
                df1 = b.xs((c), axis=1, drop_level=False)
                df1 = df1.T
                df2 = df1.melt(ignore_index=False)
                df2.columns = ["Condition","Replicate","Score"]
                m = [row for row in anovas if c == row[0][0]]
                pair1 = [(row[1],row[1]) for row in m]
                pvalue1 = [row[2] for row in m]
                print(pair1,pvalue1)
                df2 = df2[df2.Condition != wildtype1]
                subcat_palette = sns.dark_palette("#8BF", reverse=True, n_colors=5)
                plotting_parameters = {
                'data':    df2,
                'x':       'Condition',
                'y':       'Score',
                'palette': subcat_palette[1:]}
                sns.set(style="white")
                sns.set_context("paper")
                fig, ax = plt.subplots()
                ax = sns.swarmplot(x="Condition", y="Score",hue="Condition", data=df2, alpha=1,size=circ_size,palette=st.session_state.pal_type, order=order_list)
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
                        ax=ax,order=order_list)
                annotator = Annotator(ax, pair1, **plotting_parameters,order=order_list)
                annotator.set_pvalues(pvalue1)
                annotator.configure(line_width = 0)
                annotator.annotate()
                for patch in ax.artists:
                    r, g, b, a = patch.get_facecolor()
                    patch.set_facecolor((r, g, b, 0))    
                ax.set(xlabel= "Condition",ylabel='Fitness Ratio')
                if sety == "Yes":
                    ax.set_ylim(ylow,yhigh)
                plt.title(c.replace(",",".").replace("-"," "))
                if xlabes == 'No':
                    plt.xticks([])
                if  xlabes == 'Yes':
                    plt.xticks(rotation=rote)
                fig.set_size_inches(width1, height1)
                
                sns.despine()
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                              fancybox=True, shadow=False, ncol=1)
                if width1 >= 10:
                    myexpand = st.expander(label=c, expanded=False)
                    myexpand.pyplot(fig)
                    figs.append(fig)
                    cs.append(c)
                    img =(io.BytesIO())
                    
                    plt.savefig(img, format='pdf', bbox_inches='tight')
                    fn = ("ChemGAPP_"+c+"_"+iris_type+"_"+bar_swarm+".pdf")
                    temp = os.path.join(new_dir,fn)
                    with open(temp,'wb') as fid:
                        pickle.dump(fig, fid)
                    fns.append(fn)
                    but8 = myexpand.download_button(
                         label="Download image",
                         data=img,
                         file_name=fn,
                         mime="image/pdf"
                      )
                    
                if width1 < 10:
                    a= i+1
                    if (a % 2) != 0:
                        cola, colb = st.columns((1,1))
                        myexpand = cola.expander(label=c, expanded=False)
                        myexpand.pyplot(fig)
                        figs.append(fig) 
                        cs.append(c)
                        img =(io.BytesIO())
                        
                        plt.savefig(img, format='pdf', bbox_inches='tight')
                        fn =("ChemGAPP_"+c+"_"+iris_type+"_"+bar_swarm+".pdf")
                        temp = os.path.join(new_dir,fn)
                        with open(temp,'wb') as fid:
                            pickle.dump(fig, fid)
                        fns.append(fn)
                        but9 = myexpand.download_button(
                         label="Download image",
                         data=img,
                         file_name=fn,
                         mime="image/pdf"
                      )
                        
                    else:
                        myexpand = colb.expander(label=c, expanded=False)
                        myexpand.pyplot(fig)
                        figs.append(fig)
                        cs.append(c)
                        img =(io.BytesIO())
                        
                        
                        plt.savefig(img, format='pdf', bbox_inches='tight')
                        fn =("ChemGAPP_"+c+"_"+iris_type+"_"+bar_swarm+".pdf")
                        temp = os.path.join(new_dir,fn)
                        with open(temp,'wb') as fid:
                            pickle.dump(fig, fid)
                        fns.append(fn)
                        but10 = myexpand.download_button(
                         label="Download image",
                         data=img,
                         file_name=fn,
                         mime="image/pdf"
                      )
                        
    if "figures" not in st.session_state:
            st.session_state.figures = figs
    elif "figures" in st.session_state:
            st.session_state.figures = figs
    if "conds" not in st.session_state:
            st.session_state.conds = cs
    elif "conds" in st.session_state:
            st.session_state.conds = cs
    if "fnames" not in st.session_state:
            st.session_state.fnames = fns
    elif "fnames" in st.session_state:
            st.session_state.fnames = fns
    elwarn.empty()
    st.success("Images ready to download!")