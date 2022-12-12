import streamlit as st
from PIL import Image

st.set_page_config(layout="wide")
#st.markdown("# Main page ")
#st.sidebar.markdown("# Main page")

image = Image.open('logo.png')
st.image(image, use_column_width=True,output_format='PNG')
# Page title
st.markdown("""
### Thank you for using the ChemGAPP Small app!

##### ChemGAPP Small

This app is an extension of the ChemGAPP app, for the analysis of small scale chemical genomic screens.

ChemGAPP Small compares the mean colony size of within plate replicates to the mean colony size of the within plate wildtype replicates.  

--- """)
st.markdown("""
### Go to Step 1 to begin!

---

**Credits**
- App built in `Python` + `Streamlit` by [Hannah Doherty](https://github.com/HannahMDoherty/ChemGAPP)
- This work is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
---
""")