import streamlit as st
from PIL import Image

st.set_page_config(layout="wide")
#st.markdown("# Main page ")
#st.sidebar.markdown("# Main page")

image = Image.open('logo.png')
st.image(image, use_column_width=True)
# Page title
st.markdown("""
### Thank you for using the ChemGAPP GI app!

##### ChemGAPP GI 

This app is an extension of the ChemGAPP package, allowing for the calculation and visualisation of fitness ratios within genetic interaction studies.

--- """)
st.markdown("""
### Go to Step 1 to begin!

---

**Credits**
- App built in `Python` + `Streamlit` by [Hannah Doherty](https://github.com/HannahMDoherty/ChemGAPP)
- This work is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
---
""")

