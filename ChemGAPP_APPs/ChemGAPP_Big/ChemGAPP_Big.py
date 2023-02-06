import streamlit as st
import pandas as pd
from PIL import Image

st.set_page_config(layout="wide")
#st.markdown("# Main page ")
#st.sidebar.markdown("# Main page")

image = Image.open('Logo.png')
st.image(image, use_column_width=True)
# Page title
st.markdown("""
### Thank you for using the ChemGAPP Big app!

##### ChemGAPP: *Chem*ical *G*enomic *A*nalysis and *P*henotypic *P*rofiling

This app allows for the quality control analysis of chemical genomic screen data. This app aims to improve the quality of the inputted data. It achieves this via the normalisation of plate data and by performing a series of statistical analyses for the removal of detrimental replicates or conditions at chosen thresholds. Following this, it is able to score data and assign fitness scores (S-scores). The statistical analyses included are: the `Z-score test`, the `Mann-Whitney test`, and `Condition Variance`.

--- """)
st.markdown("""
### Go to Step 1 to begin!

---

**Credits**
- App built in `Python` + `Streamlit` by [Hannah Doherty](https://github.com/HannahMDoherty/ChemGAPP)

- This work is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
---
""")

