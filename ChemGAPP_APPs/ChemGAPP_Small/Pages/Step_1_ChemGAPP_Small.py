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

4- Decide whether you wish to compare mutants to a Wildype within the same condition or compare mutants to themselves in a control condition.

If you choose `Wildtype`:

Enter the name of the wild type strain. 

This must match the name given in the plate information file.
E.g. `WT`

If you choose `Control Condition`:

Enter the name of the control condition. 

This must match the iris file name after it has been adjusted for the datasets. 

E.g. for `AMPICILLIN-50mM-6_B.JPG.iris` you would input `AMPICILLIN 50mM`.

-----

5- Select which IRIS phenotype you would like to analyse. If size is selected optionally input a maximum conlony size value.

-----

6- Select how you want the plots to be grouped. Either by strain or by condition.

-----

7- Select the type of plot you wish to produce; Bar plots or Swarm plots. 

Bar plots will display 95 % confidence intervals as the statistic. Swarm plots will display ANOVA significance annotations. 

-----

8- Select your customisation options.

------

9- Click `Begin!`

------
10- To save bar plot images, click on the `Download image` button beneath the plot. This will save the image as a pdf file.

    * Files will be downloaded to your Downloads Folder

    * Do not try save an image until all plots are produced.

""")

# makes sure all figures saved in session state deleted so program starts from scratch when options are changed
def new_click():
    if 'figures' not in st.session_state:
        pass
    else:
        del st.session_state["figures"]
        del st.session_state["conds"]
        del st.session_state["fnames"]

# makes temp folder for the figure files so that they can be downloaded as pdfs and can be reloaded when other files downloaded.         
par_dir = os.path.abspath(os.getcwd())
new_dir = "Temp"
temp_path = os.path.join(par_dir, new_dir)
try:
    os.mkdir(temp_path)
except:
    pass

# checks if streamlit needs to reload figures after clicking the download button 
# following first full run through, as clicking the download button makes 
# streamlit re-run from the top and we dont want to lose our files/images
# if programme hasnt been run yet then there will be no figures in the session state,
# therefore all files in the temp folder will be deleted as these may be from a previous run.

# must set the session states as empty if not yet, set so streamlit has them intialised for later mentions
if 'figures' not in st.session_state:
        st.session_state.figures = []
        st.session_state.conds = []
        st.session_state.fnames = []
        for f in sorted(os.listdir(temp_path)):
            if f.endswith(".pdf"): 
                path = os.path.join(temp_path, f)
                os.remove(path)
# if it has been run then it reloads the data as it was displayed before clicking a 'download' button.
if st.session_state.figures:
    my_expandera = st.expander(label="Scored Dataset", expanded=False)
    my_expandera.write(st.session_state.scored_dataset, header=[0,1,2])
    elwarn2 = st.warning("Preparing plots for download, do not click 'Download' until ready.")
    for fig,i,c,fna in zip(st.session_state.figures, range(len(st.session_state.figures)),st.session_state.conds,st.session_state.fnames):
        # i == 0 is the heatmap as this is the first displayed and has different headers and layout
        if i == 0:
            myexp = st.expander(label="Heatmap:", expanded=False)
            col1a,col1b,col1c = myexp.columns((1,4,1))
            col1b.pyplot(fig)
            #finds the path to the temporary files such that they can be saved as an io.BytesIO() meaning they can be downloaded as a PDF.
            pathfna = os.path.join("Temp", fna)
            with open(pathfna,'rb') as fid:
                ax = pickle.load(fid)
                imga = io.BytesIO()
            plt.savefig(imga, format='pdf', bbox_inches='tight')
            fn =("ChemGAPP_Heatmap.pdf")
            #produces a download button to download image as pdf
            myexp.download_button(
                            label="Download image",
                            data=imga,
                            file_name=fn,
                            mime="image/pdf"
                         )
            st.write(str(st.session_state.plot_type)+"s")   
        # everything else is the bar or swarm plots
        else:
            # sets into two column layout if width of figures in smaller than 10.
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
            # set in one column layout if width of figure >= 10
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

# all these new definitions until next comment make it such that if any setting is changed
# the session states for the figures are deleted and the programme run from scratch, without losing the other settings already set.
if "iris_type" not in st.session_state:
        st.session_state.iris_type = []
def new_iris_changed():
        if st.session_state.new_iris_type:
            st.session_state.iris_type = st.session_state.new_iris_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "comp_type" not in st.session_state:
        st.session_state.comp_type = []
def new_comp_changed():
        if st.session_state.new_comp_type:
            st.session_state.comp_type = st.session_state.new_comp_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "plot_type" not in st.session_state:
        st.session_state.plot_type = []
def new_plot_changed():
        if st.session_state.new_plot_type:
            st.session_state.plot_type = st.session_state.new_plot_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "group_type" not in st.session_state:
        st.session_state.group_type = []
def new_group_changed():
        if st.session_state.new_group_type:
            st.session_state.group_type = st.session_state.new_group_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "width_type" not in st.session_state:
        st.session_state.width_type = []
def new_width_changed():
        if st.session_state.new_width_type:
            st.session_state.width_type = st.session_state.new_width_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "heatsize_type" not in st.session_state:
        st.session_state.heatsize_type = []
def new_heatsize_changed():
        if st.session_state.new_heatsize_type:
            st.session_state.heatsize_type = st.session_state.new_heatsize_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "height_type" not in st.session_state:
        st.session_state.height_type = []
def new_height_changed():
        if st.session_state.new_height_type:
            st.session_state.height_type = st.session_state.new_height_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "heat_width_type" not in st.session_state:
        st.session_state.heat_width_type = []
def new_heatwidth_changed():
        if st.session_state.new_heatwidth_type:
            st.session_state.heat_width_type = st.session_state.new_heatwidth_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "heat_height_type" not in st.session_state:
        st.session_state.heat_height_type = []
def new_heatheight_changed():
        if st.session_state.new_heatheight_type:
            st.session_state.heat_height_type = st.session_state.new_heatheight_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "colour_type" not in st.session_state:
        st.session_state.colour_type = []
def new_colour_changed():
        if st.session_state.new_colour_type:
            st.session_state.colour_type = st.session_state.new_colour_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "in_type" not in st.session_state:
        st.session_state.in_type = []
def new_in_changed():
        if st.session_state.new_in_type:
            st.session_state.in_type = st.session_state.new_in_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "out_type" not in st.session_state:
        st.session_state.out_type = []
def new_out_changed():
        if st.session_state.new_out_type:
            st.session_state.out_type = st.session_state.new_out_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]


if "info_type" not in st.session_state:
        st.session_state.info_type = []
def new_info_changed():
        if st.session_state.new_info_type:
            st.session_state.info_type = st.session_state.new_info_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "WT_type" not in st.session_state:
        st.session_state.WT_type = []
def new_WT_changed():
        if st.session_state.new_WT_type:
            st.session_state.WT_type = st.session_state.new_WT_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "vscond_type" not in st.session_state:
        st.session_state.vscond_type = []
def new_vscond_changed():
        if st.session_state.new_vscond_type:
            st.session_state.vscond_type = st.session_state.new_vscond_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "rote_type" not in st.session_state:
        st.session_state.rote_type = []
def new_rote_changed():
        if st.session_state.new_rote_type:
            st.session_state.rote_type = st.session_state.new_rote_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]


if "circ_type" not in st.session_state:
        st.session_state.circ_type = []
def new_circ_changed():
        if st.session_state.new_circ_type:
            st.session_state.circ_type = st.session_state.new_circ_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "xlab_type" not in st.session_state:
        st.session_state.xlab_type = []
def new_xlab_changed():
        if st.session_state.new_xlab_type:
            st.session_state.xlab_type = st.session_state.new_xlab_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "ylim_type" not in st.session_state:
        st.session_state.ylim_type = []
def new_ylim_changed():
        if st.session_state.new_ylim_type:
            st.session_state.ylim_type = st.session_state.new_ylim_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "ylow_type" not in st.session_state:
        st.session_state.ylow_type = []
def new_ylow_changed():
        if st.session_state.new_ylow_type:
            st.session_state.ylow_type = st.session_state.new_ylow_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "yhigh_type" not in st.session_state:
        st.session_state.yhigh_type = []
def new_yhigh_changed():
        if st.session_state.new_yhigh_type:
            st.session_state.yhigh_type = st.session_state.new_yhigh_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "pal_type" not in st.session_state:
        st.session_state.pal_type = "icefire"
def new_pal_changed():
        if st.session_state.new_pal_type:
            st.session_state.pal_type = st.session_state.new_pal_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "pal2_type" not in st.session_state:
        st.session_state.pal2_type = "bwr_r"
def new_pal2_changed():
        if st.session_state.new_pal2_type:
            st.session_state.pal2_type = st.session_state.new_pal2_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "remove_type" not in st.session_state:
        st.session_state.remove_type = []
def new_remove_changed():
        if st.session_state.new_remove_type:
            st.session_state.remove_type = st.session_state.new_remove_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "order_type" not in st.session_state:
        st.session_state.order_type = []
def new_order_changed():
        if st.session_state.new_order_type:
            st.session_state.order_type = st.session_state.new_order_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

if "max_size" not in st.session_state:
        st.session_state.max_type = []
def new_max_changed():
        if st.session_state.new_max_type:
            st.session_state.max_type = st.session_state.new_max_type
            if 'figures' not in st.session_state:
                pass
            else:
                del st.session_state["figures"]
                del st.session_state["conds"]
                del st.session_state["fnames"]

# various input buttons to upload files or make customisations.
uploaded_files = st.sidebar.file_uploader("Upload multiple IRIS files", accept_multiple_files=True)
plate_info_files = st.sidebar.file_uploader("Upload multiple plate information files", accept_multiple_files=True)
outputfile1 = st.sidebar.text_input("Enter Path and Prefix for Output Files:", on_change=new_out_changed, key="new_out_type")
comp = st.sidebar.radio(
    "Would you like to compare mutants to the wildtype within the same condition, or compare mutants to themselves within a control condition?",
    ("Wildtype","Control Condition"), on_change=new_comp_changed, key="new_comp_type")
if comp == "Wildtype":
    wildtype1 = st.sidebar.text_input("Name of the wildtype Strain:", on_change=new_WT_changed, key="new_WT_type")
if comp == "Control Condition":
    vscond = st.sidebar.text_input("Name of the control condition:", on_change=new_vscond_changed, key="new_vscond_type")
iris_type = st.sidebar.radio(
    "Select what the IRIS morphology to test:",
    ('size', 'circularity','opacity'), on_change=new_iris_changed, key="new_iris_type")



if iris_type == 'size':
    max_size = st.sidebar.text_input('Optional: Max value for colony '+iris_type, on_change=new_max_changed, key="new_max_type")
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
heatsize = st.sidebar.slider('Heatmap Label Size', 0, 10, 4, on_change=new_heatsize_changed, key="new_heatsize_type")
heat_width = st.sidebar.slider('Heatmap Width', 1, 30, 10, on_change=new_heatwidth_changed, key="new_heatwidth_type")
heat_height = st.sidebar.slider('Heatmap Height', 1, 30, 10,on_change=new_heatheight_changed, key="new_heatheight_type")
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

# sets off code only once all settings have ben set and use presses begin.
if complete:
    #removes anything saved from previous runs
    del st.session_state["figures"]
    del st.session_state["conds"]
    del st.session_state["fnames"]
    # if user has set which strains to remove this produces a list of the strains split by ;
    if remove_y_n == 'Yes':
        inclstrains = txt.split(";")
    # if user has set order of strains this produces a list of the strains in prefered order split by ;
    if order_type == 'Yes':
        order_list = order_set.split(";")
    #sets the chosen iris types and plot type in session state
    if st.session_state.iris_type == []:
        st.session_state.iris_type = iris_type
    if st.session_state.plot_type == []:
        st.session_state.plot_type = bar_swarm
    m = None
    # cycles through iris files and uses filename to produce column headers.
    for f in uploaded_files:
        g = pd.read_csv(f,
                comment='#',
                index_col=[0, 1],
                sep='\t')
        if m is None:
            m = g[iris_type]
            # replaces commas to periods, since these are not allowed within the filename and dashes to spaces
            m.name = (int(f.name.split('-')[2].split('_')[0]),
                      '-'.join(f.name.split('.')[0].split('-')[:2]).replace(",",".").replace("-"," "),
                      f.name.split('.')[0].split('_')[1])
            m = m.to_frame()
        else:
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
    df_outer_norm = pd.DataFrame(index=m.index)

    n_array = np.array(n)
    ardf = np.zeros((mlen, 1))
    ardf.shape = (mlen,1)
    for c1,c2,i in zip(sorted(df_pmm.columns),sorted(df_outer.columns), range(len(m.columns))):
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
                # to the median colony size of the entire datatset. 
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
    ardf.columns = (pd.MultiIndex.from_tuples(sorted(m.columns)))
    ardf = ardf[sorted(ardf)]
    ardf.to_csv(outputfile1+"_Normalised_dataset.csv")

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

    # looks for digits within the file names so that order of 
    # plate info files is sorted plate1, plate2, plate3 etc. not plate1, plate10, plate2 etc.
    def atoi(text):
        return int(text) if text.isdigit() else text
    def natural_keys(text):
        return [ atoi(c) for c in re.split(r'(\d+)', text) ]

    #adds the plate files to a list
    plate_DFs = []
    for i in plate_info_files:
        p = plate_info(i)
        plate_DFs.append(p)
    # makes set of all source plate numbers within dataset
    plates = {x[0] for x in ardf.columns}
    nm4 = ardf.copy(deep=False)
    columns1 = {x[1:4] for x in nm4.columns}
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
    if max_size != None:
        df_with_strains[df_with_strains > float(max_size)] = np.nan
    
    df_with_strains= df_with_strains.reset_index()
    # removes rows named "EMPTY"
    df_with_strains = df_with_strains[df_with_strains['Gene'] != "EMPTY"]
    # if specified to remove certain strains will read the input and split by ";" to remove rows with specified names.
    if remove_y_n == 'Yes':
        for stra in inclstrains:
            df_with_strains = df_with_strains[df_with_strains['Gene'] != stra]
    df_with_strains = df_with_strains.set_index(['Gene'])
    df_with_strains.to_csv(outputfile1+"_Final_dataset.csv")
    figs = []
    imgs = []
    fns = []
    cs = []
    if comp == "Wildtype":
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
            wt_mean = float(df1.loc[wildtype1])
            # divides each averaged colony size by the mean wildtype colonysize for that condition
            for ind, j in zip(df1.index,range(len(ar1))):
                ar2[j] = (float(ar1[j])/wt_mean)
            #adds scored column to new dataframe
            df3 = np.concatenate((df3,ar2), axis=1)    
        df3 = pd.DataFrame(df3, index=group.index, columns=group.columns)
        df3.index.name = None
        # drops the WT row from the scored dataset.
        df3 = df3.drop(wildtype1,axis=0)
        my_expandera = st.expander(label="Scored Dataset", expanded=False)
        my_expandera.write(df3, header=[0,1,2])
        df3.to_csv(outputfile1+"_Scored_Dataset.csv")
        if "scored_dataset" not in st.session_state:
            st.session_state.scored_dataset = df3
        elif "scored_dataset" in st.session_state:
            st.session_state.scored_dataset = df3
        elwarn = st.warning("Preparing plots for download, do not click 'Download' until ready.")
        conditions = {x[0] for x in df3.columns}
        clen = len(conditions)
        slen = len(df3)

        #preps lists for the appending various attributes such that they can be reloaded if a plot is downloaded.

        #produces heatmap based on the same dataset used for the barplots
        # transposes data such that rows are conditions
        b = df3.T
        # averages replicate scores
        c = b.groupby(level=0).mean()
        #transposes back for plotting of the heatmap
        x = c.T
        for i in range(len(x.columns)):
            x = x.rename(columns={x.columns[i]: x.columns[i].replace(",",".")})
        c = x.T
        sns.set(style="white")
        sns.set_context("paper")
        fig, ax = plt.subplots()
        fig.set_size_inches(heat_width, heat_height)
        colp = sns.color_palette(st.session_state.pal2_type, as_cmap=True)
        if heatsize == 0:
            ax = sns.heatmap(c,center=1,cmap=colp,square=False,annot=False,fmt=".2f",linewidths=.1,linecolor='0',vmin=0, vmax=2)
        else:
            ax = sns.heatmap(c,center=1,cmap=colp,square=False,annot=True,annot_kws={"size": heatsize},fmt=".2f",linewidths=.1,linecolor='0',vmin=0, vmax=2)
        myexp = st.expander(label="Heatmap:", expanded=False)
        col1a,col1b,col1c = myexp.columns((1,4,1))
        col1b.pyplot(fig)
        figs.append(fig)
        c="_"
        cs.append(c)
        img = io.BytesIO()
        fn =("ChemGAPP_Heatmap.pdf")
        temp = os.path.join(new_dir,fn)
        #writes a temporary file into pickle of the figure such that it can be reloaded later
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
        # if selected to group by condition, and bar plot then produces bar plot for each condition. 
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
            # if swarmplot selected produces swarmplots instead. Here the data is scored differently.          
            if bar_swarm == 'Swarm Plot':
                df_swarm1 = df_with_strains.reset_index()
                #makes df including just the wildtype values
                df_swarm2 = df_swarm1[df_swarm1['Gene'] == wildtype1]
                # df_swarm3 = df_swarm1[df_swarm1['Gene'] != WT]
                # adds "0_" to WT column name such that it is always sorted first 
                # when sorting alphabetically, necessary for the ANOVA tests.
                df_swarm1['Gene'] = df_swarm1['Gene'].str.replace(wildtype1,("0_"+wildtype1))
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
                    for ind, j in zip(df1.index,range(len(ar1))):
                        #divides each individual colony size by the WT_mean, including for the WT values.
                        ar2[j] = (float(ar1[j])/mean)
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
                            ax=ax1,order=order_list)
                    #adds anova asterisk type annotations for sigificance. 
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
                    plt.title(c)
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
        #here groups by strain and not condition.
        if plot_type == 'Strain':
            if bar_swarm == 'Bar Plot':
                #transposes dataset such that it can be grouped by strains.
                b = df3.T
                conditions = {x for x in b.columns}
                #splits and iterates by strain.
                for c,i in zip(sorted(conditions),range(len(conditions))):
                    df1 = b.xs((c), axis=1, drop_level=False)
                    df1 = df1.T
                    #melts dataset such that all scores in one column with key being for the conditions
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
                #as above with condition grouped swarmplots, 
                # produces dataset with each individual colonzy size divided by WT mean for that condition. 
                # then calulates anova and tukey-hsd and takes p values for WT vs each strain. 
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
                # produces zeros matrix that can be filled for [(gene1,gene1),condition,p-value], 
                # need to compare pairs to themselves since plots sorted by strain then plot based on the conditions split by gene name
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
                    #takes only those compared to the WT
                    stats1 = res.tukey_summary[res.tukey_summary['group1'] == ("0_"+wildtype1)]
                    # makes list of the current condition same length as number of genes - 1. 
                    # Allows for zipping to the pairs and p-values for the matrix.
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
                #transposes dataset such that can be split by strain instead of condition.
                b = df_swarm4.T
                conditions = {x for x in b.columns}
                for c,i in zip(sorted(conditions),range(len(conditions))):
                    df1 = b.xs((c), axis=1, drop_level=False)
                    df1 = df1.T
                    df2 = df1.melt(ignore_index=False)
                    df2.columns = ["Condition","Replicate","Score"]
                    #finds the annotation pair (condition,condition) 
                    # and pvalue for the current strain for annotation of significance scores.
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
                    plt.title(c)
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
    if comp == "Control Condition":
        df_grouped = df_with_strains.groupby(level=0).mean()
        if vscond in df_grouped.columns:
            vscond = vscond
        else:
            vscond = (vscond+" ")
            st.session_state.vscond_type = vscond
        df_grouped["versus"] = np.mean(df_grouped[vscond],axis=1)
        cols = df_grouped.columns
        df_bar = pd.DataFrame(index=df_grouped.index)
        for c in cols:     
            df1 = pd.DataFrame(df_grouped.xs((c), axis=1, drop_level=True))
            df2 = pd.DataFrame(df_grouped.xs(("versus"), axis=1, drop_level=True))
            df1.columns = [f'{i}_{j}' for i,j in df1.columns]
            df2 = df2.rename(columns={df2.columns[0]: df1.columns[0]})
            df5 = pd.DataFrame(df1/df2, index = df_grouped.index)
            df_bar = pd.concat([df_bar,df5],axis=1)

        df_bar.columns = df_bar.columns.str.split('_', expand=True)    
        for x in list(df_bar.xs((vscond), axis=1, drop_level=False).columns):
            df_bar = df_bar.drop(x, axis=1)
        df_bar = df_bar.drop("versus", axis=1)
        my_expandera = st.expander(label="Scored Dataset", expanded=False)
        my_expandera.write(df_bar, header=[0,1,2])
        df_bar.to_csv(outputfile1+"_Scored_Dataset.csv")
        if "scored_dataset" not in st.session_state:
            st.session_state.scored_dataset = df_bar
        elif "scored_dataset" in st.session_state:
            st.session_state.scored_dataset = df_bar
        
        elwarn = st.warning("Preparing plots for download, do not click 'Download' until ready.")
        b = df_bar.T
        # averages replicate scores
        c = b.groupby(level=0).mean()
        #transposes back for plotting of the heatmap
        x = c.T
        for i in range(len(x.columns)):
            x = x.rename(columns={x.columns[i]: x.columns[i]})
        c = x.T   
        sns.set(style="white")
        sns.set_context("paper")
        fig, ax = plt.subplots()
        fig.set_size_inches(heat_width, heat_height)
        colp = sns.color_palette(st.session_state.pal2_type, as_cmap=True)
        if heatsize == 0:
            ax = sns.heatmap(c,center=1,cmap=colp,square=False,annot=False,fmt=".2f",linewidths=.1,linecolor='0',vmin=0, vmax=2)
        else:
            ax = sns.heatmap(c,center=1,cmap=colp,square=False,annot=True,annot_kws={"size": heatsize},fmt=".2f",linewidths=.1,linecolor='0',vmin=0, vmax=2)
        myexp = st.expander(label="Heatmap:", expanded=False)
        col1a,col1b,col1c = myexp.columns((1,4,1))
        col1b.pyplot(fig)
        figs.append(fig)
        c="_"
        cs.append(c)
        img = io.BytesIO()
        fn =("ChemGAPP_Heatmap.pdf")
        temp = os.path.join(new_dir,fn)
        #writes a temporary file into pickle of the figure such that it can be reloaded later
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


        cols = df_with_strains.columns
        df_with_strains["versus"] = np.mean(df_with_strains[vscond],axis=1)
        means = df_with_strains["versus"].groupby(level=0).mean()
        df_swarm4 = pd.DataFrame(index=df_with_strains.index)
        df_with_strains["means"] = means
        for c in cols:     
            df1 = pd.DataFrame(df_with_strains.xs((c), axis=1, drop_level=True))
            df2 = pd.DataFrame(df_with_strains.xs(("means"), axis=1, drop_level=True))
            df1.columns = [f'{i}_{j}' for i,j in df1.columns]
            df2 = df2.rename(columns={df2.columns[0]: df1.columns[0]})
            df5 = pd.DataFrame(df1/df2, index = df_with_strains.index)
            df_swarm4 = pd.concat([df_swarm4,df5],axis=1)
        df_swarm4.columns = df_swarm4.columns.str.split('_', expand=True)
        df_swarm4.index.name = None
        df_swarm4 = df_swarm4.rename(columns={vscond:("0_"+vscond)},level=0)
        
        st.write(bar_swarm+"s")
        if plot_type == 'Condition':
            if bar_swarm == 'Bar Plot':
                conditions = {x[0] for x in df_bar.columns}
                for c,i in zip(sorted(conditions),range(len(conditions))):
                    df1 = df_bar.xs((c), axis =1, drop_level=False)
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
                    
            if bar_swarm == "Swarm Plot":
                t_swarm = df_swarm4.T
                t_swarm = t_swarm.sort_index()
                conditions = {x[0] for x in df_swarm4.columns}
                genes = {x for x in df_swarm4.index}
                leng = len(conditions)-1
                # produces zeros matrix that can be filled for [(gene1,gene1),condition,p-value], 
                # need to compare pairs to themselves since plots sorted by strain then plot based on the conditions split by gene name
                anovas = np.zeros((leng, 3),dtype=object)
                anovas.shape = (leng,3)
                for c in sorted(genes):
                    df1 = t_swarm.xs((c), axis =1, drop_level=False)
                    df2 = df1.melt(ignore_index=False)
                    df2 = df2.reset_index()

                    df2.columns = ["Condition","Replicate","Strain","Score"]
                    res = stat()
                    res.tukey_hsd(df=df2, res_var='Score', xfac_var='Condition', anova_model='Score ~ C(Condition)')
                    #takes only those compared to the WT
                    stats1 = res.tukey_summary[res.tukey_summary['group1'] == ("0_"+vscond)]
                    # makes list of the current condition same length as number of genes - 1. 
                    # Allows for zipping to the pairs and p-values for the matrix.
                    genes = ([c] * leng)
                    pairs = list(zip(genes,genes))
                    pvalues = list(stats1['p-value'])
                    cond = list(stats1.group2)
                    grouped = list(zip(pairs,cond,pvalues))
                    grouped = np.array(grouped,dtype=object)
                    anovas = np.append(anovas,grouped, axis=0)
                anovas = anovas[leng:len(anovas)]
                print(anovas)

                for x in list(df_swarm4.xs(("0_"+vscond), axis=1, drop_level=False).columns):
                    df_swarm4 = df_swarm4.drop(x, axis=1)
                conditions = {x[0] for x in df_swarm4.columns}
                #iterates through conditions and makes sub-dataset including all replicate plates for same condition.
                for c,i in zip(sorted(conditions),range(len(conditions))):
                    df1 = df_swarm4.xs((c), axis =1, drop_level=False)
                    #metls dataset so all scores are in sigular column
                    df2 = df1.melt(ignore_index=False)
                    df2= df2.reset_index()
                    df2.columns = ["Strain","Condition","Replicate","Score"]
                    print(c)
                    m = [row for row in anovas if c == row[1]]
                    pair1 = [(row[0][0],row[0][1]) for row in m]
                    pvalue1 = [row[2] for row in m]
                    subcat_palette = sns.dark_palette("#8BF", reverse=True, n_colors=5)
                    plotting_parameters = {
                                'data':    df2,
                                'x':       'Strain',
                                'y':       'Score',
                                'palette': subcat_palette[1:]}
                    sns.set(style="white")
                    sns.set_context("paper")
                    fig, ax = plt.subplots()
                    ax = sns.swarmplot(x="Strain", y="Score",hue="Strain", data=df2, alpha=1,palette=st.session_state.pal_type, order=order_list)
                    ## adds mean bars
                    ax = sns.boxplot(showmeans=True,
                            meanline=True,
                            meanprops={'color': 'k', 'ls': '-', 'lw': 1},
                            medianprops={'visible': False},
                            whiskerprops={'visible': False},
                            zorder=10,
                            x="Strain", y="Score",data=df2,
                            showfliers=False,
                            showbox=False,
                            showcaps=False,order=order_list,
                            ax=ax)
                    #adds anova asterisk type annotations for sigificance. 
                    annotator = Annotator(ax, pair1, **plotting_parameters)
                    annotator.set_pvalues(pvalue1)
                    annotator.configure(line_width = 0)
                    annotator.annotate()
                    for patch in ax.artists:
                        r, g, b, a = patch.get_facecolor()
                        patch.set_facecolor((r, g, b, 0))    
                    ax.set(xlabel= "Strain",ylabel='Fitness Ratio')
                    if sety == "Yes":
                        ax.set_ylim(ylow,yhigh)
                    plt.title(c)
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
                #transposes dataset such that it can be grouped by strains.
                b = df_bar.T
                # makes set of all strains. 
                strain_list = {x for x in b.columns}
                #splits and iterates by strain.
                for c,i in zip(sorted(strain_list),range(len(strain_list))):
                    df1 = b.xs((c), axis=1, drop_level=False)
                    df1 = df1.T
                    #melts dataset such that all scores in one column wit
                    # h key being for the conditions
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
                t_swarm = df_swarm4.T
                t_swarm = t_swarm.sort_index()
                conditions = {x[0] for x in df_swarm4.columns}
                genes = {x for x in df_swarm4.index}
                leng = len(conditions)-1
                # produces zeros matrix that can be filled for [(gene1,gene1),condition,p-value], 
                # need to compare pairs to themselves since plots sorted by strain then plot based on the conditions split by gene name
                anovas = np.zeros((leng, 3),dtype=object)
                anovas.shape = (leng,3)
                for c in sorted(genes):
                    df1 = t_swarm.xs((c), axis =1, drop_level=False)
                    df2 = df1.melt(ignore_index=False)
                    df2 = df2.reset_index()

                    df2.columns = ["Condition","Replicate","Strain","Score"]
                    res = stat()
                    res.tukey_hsd(df=df2, res_var='Score', xfac_var='Condition', anova_model='Score ~ C(Condition)')
                    #takes only those compared to the WT
                    stats1 = res.tukey_summary[res.tukey_summary['group1'] == ("0_"+vscond)]
                    # makes list of the current condition same length as number of genes - 1. 
                    # Allows for zipping to the pairs and p-values for the matrix.
                    genes = ([c] * leng)
                    pairs = list(zip(stats1.group2,stats1.group2))
                    pvalues = list(stats1['p-value'])
                    cond = list(genes)
                    grouped = list(zip(pairs,cond,pvalues))
                    grouped = np.array(grouped,dtype=object)
                    anovas = np.append(anovas,grouped, axis=0)
                anovas = anovas[leng:len(anovas)]
                for x in list(df_swarm4.xs(("0_"+vscond), axis=1, drop_level=False).columns):
                    df_swarm4 = df_swarm4.drop(x, axis=1)
                b = df_swarm4.T
                conditions = {x for x in b.columns}
                for c,i in zip(sorted(conditions),range(len(conditions))):
                    df1 = b.xs((c), axis=1, drop_level=False)
                    df1 = df1.T
                    df2 = df1.melt(ignore_index=False)
                    df2.columns = ["Condition","Replicate","Score"]
                    m = [row for row in anovas if c == row[1]]
                    pair1 = [(row[0][0],row[0][1]) for row in m]
                    pvalue1 = [row[2] for row in m]
                    subcat_palette = sns.dark_palette("#8BF", reverse=True, n_colors=5)
                    plotting_parameters = {
                                'data':    df2,
                                'x':       'Condition',
                                'y':       'Score',
                                'palette': subcat_palette[1:]}
                    sns.set(style="white")
                    sns.set_context("paper")
                    fig, ax = plt.subplots()
                    ax = sns.swarmplot(x="Condition", y="Score",hue="Condition", data=df2, alpha=1,palette=st.session_state.pal_type, order=order_list)
                    ax = sns.boxplot(showmeans=True,
                            meanline=True,
                            meanprops={'color': 'k', 'ls': '-', 'lw': 1},
                            medianprops={'visible': False},
                            whiskerprops={'visible': False},
                            zorder=10,
                            x="Condition", y="Score",data=df2,
                            showfliers=False,
                            showbox=False,
                            showcaps=False,order=order_list,
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
                    if sety == "Yes":
                        ax.set_ylim(ylow,yhigh)
                    plt.title(c)
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