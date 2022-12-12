import streamlit as st
import pandas as pd
import numpy as np
import os
import string
import scipy.stats
import re
import io
import streamlit as st
from PIL import Image

st.set_page_config(layout="wide")
col1,col2,col3 = st.columns((1,9,1))
image = Image.open('title.png')
col2.image(image, use_column_width=True,output_format='PNG')
st.markdown("""
---
This step aids in the renaming of files for use within ChemGAPP.

Upload a set of replicates which are for the same condition, concentration and source plate number. This page will reformat their name to the correct format for ChemGAPP.

--- """)
cola,colb = st.columns((1,1))
inp = cola.text_input('Path to folder containing inputted files:')
inp = os.path.expanduser(inp)
uploaded_files = colb.file_uploader("Upload files for replicates", accept_multiple_files=True)
st.markdown("""
---
##### Input name fields:
 """)

c1,c2,c3 = st.columns((1,1,1))
condition = c1.text_input('Condition:')
concentration = c2.text_input('Concentration:')
platenumber = c3.text_input('Source plate number:')
complete = st.button(label="Begin!")
alphabet = list(string.ascii_uppercase)
names = []
for file in uploaded_files:
    names.append(file.name)

if complete:
    element1 = cola.warning("Renaming Files")
    for name, replicate in zip(names,alphabet):
        old_name = (inp+"/"+name)
        new_name = (inp+"/"+condition+"-"+concentration.replace(".",",")+"-"+platenumber+"_"+replicate+".JPG.iris")
        os.rename(old_name, new_name)
    element1.empty()
    element2 = cola.success("Files Renamed")

