import streamlit as st
import pandas as pd
import numpy as np
import os
import re

st.set_page_config(layout="wide")
st.title('ChemGAPP GI: Interaction Scores')
inptpath = st.sidebar.text_input("Path to the folder where you wish to save your files:")
if "inputpath" not in st.session_state:
            st.session_state.inputpath = inptpath
elif "inputpath" in st.session_state:
            st.session_state.inputpath = inptpath
#col1, col2 = st.columns((1,1))
my_expandera = st.expander(label="Instructions:", expanded=True)
my_expandera.markdown("""

1. Input a path to your desired output folder. This is where your files will be stored. 
```
    E.g. ~/Desktop/GI_Files
```

* Saved files will be:
    
    * X_Colony_sizes.csv
    
    * X_Interaction_Scores.csv
        
        * *Where X = Secondary Gene Name*
        * If using multiple sets in one plate, seperate files will be produced for each set. 

```
If you had:

MexA::MexB
MexA::MexY
MexA::OmpM

MexA would be the Primary Gene
MexB, MexY or OmpM would be the Secondary Gene.

```
-----

2. Upload your IRIS files.

Ensure IRIS file names are in the format: `SecondaryGeneName_replicate.JPG.iris`

E.g. `MexB_A.JPG.iris`


-----

3. Upload plate information files. 

* These should be txt files with the format:
    
| Row | Column | Strain | Replicate | Order | Set |
|  ------------- | :-------------: | :-------------: | :-------------: | :-------------: |:-------------: |
| 1 | 1 |WT | 1 | 0 | 1 |
| 1 | 2 |Primary Gene1 | 1 | 1 | 1 |
| 1 | 3 |WT | 2 | 0 | 1 |
| 1 | 4 |Primary Gene1 | 2 | 1 | 1 |
| ... |...|... |
| 2 | 1 |Secondary Gene | 1 | 2 | 1 |
| 2 | 2 |Primary Gene1::Secondary Gene | 1 | 3 | 1 |
| 2 | 3 |Secondary Gene| 2 | 2 | 1 |
| 2 | 4 |Primary Gene1::Secondary Gene | 2 | 3 | 1 |
| ... |...|... |


* The `Replicate` column tells the software which of the four different strains to group together as a replicate. 
    
    * Each of the four mutants you wish to group must be assigned the same number. 
    * If you have multiple sets, make sure the replicate number starts at 1 for each of the sets. 
    * If you have across plate replicates, ensure that the replicate numbers continue on from each other. 
    
    ```
    E.g
    📦Project
    ┣ 📂 MexB
    ┃ ┣ 📜 MexB_A.JPG.iris
    ┃ ┣ 📜 MexB_B.JPG.iris
    ┃ ┣ 📜 MexB_C.JPG.iris
    ┃ ┗ 📂 Plate_info
        ┣ plateA.txt
        ┣ plateB.txt
        ┗ plateC.txt
    
    plateA = Replicates 1-96
    plateB = Replicates 97-192
    plateC = Replicates 193-288  

    ```

* The `Order` column tells the software which strain is the wildtype, primary gene, secondary gene, and double knockout.


    * The `Order` column must match these values:

        * Wildtype = 0
        * Primary gene = 1
        * Secondary gene = 2
        * Double knockout = 3

* The `Set` column tells the software which areas of the plate are designated for different secondary genes.
    * If just using one set, apply `1` to all rows.
    * If using multiple sets, E.g:
        * MexA::MexB = `1`
        * MexA::MexY = `2`
        * MexA::OmpM = `3`
        
-----

4. Press Begin. 
""")

#inputfile1 = st.sidebar.file_uploader("Upload Iris File:", accept_multiple_files=False)
#ident1 = st.sidebar.file_uploader("Upload Plate Info File:", accept_multiple_files=False)
#indir1 = st.sidebar.text_input("Path to IRIS Files:")
#infodir1 = st.sidebar.text_input("Path to Plate Information Files:")
uploaded_files = st.sidebar.file_uploader("Upload IRIS files", accept_multiple_files=True)
plate_info_files = st.sidebar.file_uploader("Upload plate information files", accept_multiple_files=True)

complete = st.sidebar.button(label="Begin!")
if complete:
            def GI_dataset(PATH):   
                m = None
                # cycles through iris files and uses filename to produce column headers.
                for f in uploaded_files:
                    g = pd.read_csv(f,
                            comment='#',
                            index_col=[0, 1],
                            sep='\t')
                    if m is None:
                        try:
                            m = g['colony size']
                        except:
                            m = g['size']
                            m.name = (f.name.split('_')[0],f.name.split('.')[0].split('_')[1])
                            m = m.to_frame()
                    else:
                        try:
                            m1 = g['colony size']
                        except:
                            m1 = g['size']
                            m1.name = (f.name.split('_')[0], f.name.split('.')[0].split('_')[1])
                            m1 = m1.to_frame()
                        #sets them in a dataframe grouped by the secondary gene name and then replicates.
                        m = m.join(m1, how='inner')
                # renames the plate info file columns and adds row and column to index before sorting by index.
                def plate_info(file):
                    p = pd.read_table(file)
                    if len(p.columns) != 6:
                        print("ERROR: Info file formatting incorrect! Should be 6 columns: Row,Column,Strain,Replicate,Order,Set")
                    if len(p.columns) == 6:
                        p.columns = ['row','column','strain','Replicate','Order','Set']
                        p = p.set_index(['row','column'])
                        p = p.sort_index()
                    return p
                
                # looks for digits within the file names so that order of 
                # plate info files is sorted plate1, plate2, plate3 etc. not plate1, plate10, plate2 etc.        
                def atoi(text):
                    return int(text) if text.isdigit() else text

                def natural_keys(text):
                    return [ atoi(c) for c in re.split(r'(\d+)', text.name) ]
                
                plate_info_files.sort(key=natural_keys)
                #adds the plate files to a list
                plate_DFs = []
                for i in plate_info_files:
                    p = plate_info(i)
                    plate_DFs.append(p)
                plates = {x[0:2] for x in m.columns}
                m2 = pd.DataFrame(m.copy(deep=False))
                columns1 = [x[0] for x in m2.columns]
                df_with_strains = pd.DataFrame(columns = [columns1[0],'strain','Replicate','Order','Set'])
                # splits by secondary gene name 
                # adds the gene names, replicate number, set and order to the colony size data based on row and column
                for a, n in zip(plate_DFs, sorted(plates)):
                    df1 = (m.xs((n), axis =1, drop_level=True))
                    df1 = pd.DataFrame(df1)
                    df2 = pd.merge(df1,a,left_index=True, right_index=True)
                    df2.columns = [n[0],'strain','Replicate','Order','Set']
                    df_with_strains = pd.concat([df_with_strains, df2], ignore_index=True)
                df_with_strains = df_with_strains.rename(columns={'strain': 'Gene'})
                df2 = df_with_strains.set_index(['Set','Replicate','Order','Gene'])
                set = {x[0] for x in df2.index}
                #splits by set, thus different gene sets have their own output files. 
                for s in sorted(set):
                    df = pd.DataFrame()
                    df1 = df2.xs((s), axis=0, drop_level=True)
                    replicate = {x[0] for x in df1.index}

                    for r in sorted(replicate):
                        dfb = df1.xs((r), axis=0, drop_level=False)
                        dfb = dfb.sort_index(level="Order")
                        if dfb.iloc[0][0] == 0:
                            dfb.iloc[0][0] = np.nan
                        #calculates the fitness ratios and double expected ratios and adds them to the dataset
                        name = (dfb.iloc[1][0]/dfb.iloc[0][0],dfb.iloc[2][0]/dfb.iloc[0][0],(dfb.iloc[1][0]/dfb.iloc[0][0])*(dfb.iloc[2][0]/dfb.iloc[0][0]),dfb.iloc[3][0]/dfb.iloc[0][0])
                        columns = list([dfb.iloc[1].name[2],"Secondary Gene","Double Expected","Double Observed"])
                        data = []
                        zipped = zip(columns, name)
                        a_dictionary = dict(zipped)
                        data.append(a_dictionary)
                        df = df.append(data, True)
                    df1.to_csv((PATH+"/"+dfb.iloc[2].name[2]+"_Colony_sizes.csv"))
                    df.to_csv((PATH+"/"+dfb.iloc[2].name[2]+"_Interaction_Scores.csv"))
                #return df2,df
            GI_dataset(inptpath)
            st.success("Files saved to: "+inptpath)
            st.success("Upload another set or continue to step 2!")