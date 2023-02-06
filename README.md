# ChemGAPP: A Package for *Chem*ical *G*enomic *A*nalysis and *P*henotypic *P*rofiling.

## Table of Contents
1. [Introduction](#content)
2. [License](#license)
3. [Installation](#installation)
4. [Manual](#manual)
5. [Test files](#package)
6. [Contact](#contact)


----

### Introduction <a name="content"></a>

ChemGAPP encompasses three modules, each with a dedicated Streamlit App.



#### ChemGAPP Big
This package is for the quality control analysis of large chemical genomic screen data. Since many issues can arise during chemical genomic screens, such as pinning mistakes and edge effects, this package aims to improve the quality of the inputted data. It achieves this via the normalisation of plate data and by performing a series of statistical analyses for the removal of detrimental replicates or conditions. Following this, it is able to score data and assign fitness scores (S-scores). The statistical analyses included are: the Z-score test, the Mann-Whitney test, and condition variance analysis. 

#### ChemGAPP Small

ChemGAPP_Small has been produced to deal with small scale chemical genomic screens where replicates are within the plate. This differs from large chemical genomic screen where replicates are across multiple plates. ChemGAPP Small produces three types of plots, a heatmap, bar plots and swarm plots. For the bar plot and heatmap, ChemGAPP Small compares the mean colony size of within plate replicates to the mean colony size of the within plate wildtype replicates, producing a fitness ratio. The bar plots are then optionally grouped by strain or by condition. The heatmap displays all conditions and strains. For the swarm plots each mutant colony size is divided by the mean colony size of the wildtype, to produce the fitness ratio. A one-way ANOVA and Tukey-HSD analysis determines the significance in difference between each mutant fitness ratio distribution and the wildtype fitness ratio distribution. colony size can be substituted for any IRIS phenotype.


#### ChemGAPP GI

ChemGAPP GI focuses on the analysis of genetic interaction studies. ChemGAPP GI calculates the fitness ratio of two single mutant strains and a double knockout in comparison to the wildtype. It also calculates the expected double knockout fitness ratio for comparison to the observed fitness ratio. ChemGAPP GI displays this data as a bar plot. 

------

### License <a name="license"></a>

This work is licensed under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

----

### Installation <a name="installation"></a>

There are two ways to run the tool, via the package files or by running the streamlit applications.

### Package:

- The easiest way to install the package is via `pip` :

```
pip install ChemGAPP
```

- The package module files can also be downloaded and run as a python or bin files. To download the files run:
```
git clone https://github.com/HannahMDoherty/ChemGAPP
```
- Then run:

```
pip install -r requirements.txt
```
or 
```
pip3 install -r requirements.txt
```
### Streamlit applications: 
  
  - The graphical interfaces can be run using streamlit. 
    - First download the APP files
    - Ensure you have Python 3. This can be easily downloaded by installing Anaconda-Naviagator (available here: https://www.anaconda.com/products/distribution)
  - For Mac:
    - Open terminal to navigate into the `ChemGAPP_APPs` directory and follow the below instructions.
  - For Windows:
    - Within Anaconda-Navigator click on the `Environments` tab, then for `base (root)` click the green play button and select `Open Terminal`. Then navigate into the `ChemGAPP_APPs` directory and follow the below instructions.


  - For ChemGAPP Big, you need to navigate into the `ChemGAPP_APPs` directory and then the `ChemGAPP_Big` directory and launch the application, using the following command in terminal:

```
pip install -r requirements.txt
streamlit run ChemGAPP_Big.py
```

The commmand provides a link to the following front web application:
<p align="center">
<img src="https://github.com/HannahMDoherty/ChemGAPP/blob/main/ChemGAPP_APP.png" width="500" />
</p>

  - For ChemGAPP_Small, you need to navigate into the `ChemGAPP_APPs` directory and then the `ChemGAPP_Small` directory and launch the application, using the following command in terminal:

```
pip install -r requirements.txt
streamlit run chemgapp_small
```

The commmand provides a link to the following front web application:
<p align="center">
<img src="https://github.com/HannahMDoherty/ChemGAPP/blob/main/ChemGAPP_Small.png" width="500" />
</p>

  - For ChemGAPP_GI, you need to navigate into the `ChemGAPP_APPs` directory and then the `ChemGAPP_GI` directory and launch the application, using the following command in terminal:

```
pip install -r requirements.txt
streamlit run ChemGAPP_GI.py
```
The commmand provides a link to the following front web application:
<p align="center">
<img src="https://github.com/HannahMDoherty/ChemGAPP/blob/main/ChemGAPP_GI.png" width="500" />
</p>

----

### Manual <a name="manual"></a>

#### Python Modules

##### ChemGAPP Big
1. [iris_to_dataset](#Iris_to_dataset)
2. [check_normalisation](#Check_Normalisation)
3. [z_score](#Z_Score)
4. [z_score_count](#Z_Score_Count)
5. [mw_plate_level](#Mann_Whitney_Plate_Level)
6. [mw_condition_level](#Mann_Whitney_Condition_Level)
7. [condition_variance](#Condition_Variance)
8. [pass_fail_conditions](#Pass_Fail_Conditions)
9. [pass_fail_plates](#Pass_Fail_Plates)
10. [bar_plot_plates](#Bar_Plot_Plate)
11. [bar_plot_conditions](#Bar_Plot_Conditions)
12. [mw_plates_to_remove](#MW_Plates_to_Remove)
13. [z_plates_to_remove](#Z_Plates_to_Remove)
14. [MW_Conditions_to_Remove](#MW_Conditions_to_Remove)
15. [variance_conditions_to_remove](#Variance_Conditions_to_Remove)
16. [s_scores](#S_Scores)
17. [add_gene_names](#Add_Gene_Names)
18. [cosine_similarity](#Cosine_Similarity)

##### ChemGAPP Small

1. [chemgapp_small](#ChemGAPP_Small)

##### ChemGAPP GI

1. [gi_dataset](#GI_Dataset)
2. [gi_barplot](#GI_BarPlot)

#### Streamlit APPs

##### ChemGAPP Big

1. [Step_1_Normalisation.py](#Step_1_Normalisation)
2. [Step_2_Threshold_Selector.py](#Step_2_Threshold_Selector)
3. [Step_3_S_Score_Calculator.py](#Step_3_S_Score_Calculator)
4. [Step_4_Dataset_Comparison.py](#Step_4_Dataset_Comparison)

##### ChemGAPP Small

1. [Step_1_chemgapp_small](#Step_1_ChemGAPP_Small)

##### ChemGAPP GI

1. [Step_1_Interaction_Scores.py](#Step_1_Interaction_Scores)
2. [Step_2_Bar_Plot.py](#Step_2_Bar_Plot)


###### If downloaded via pip commands can be initiated from any folder. The help instruction is called using -h option. E.g:

```
iris_to_dataset [-h] [-p PATH] [-o OUTPUTFILE]
```

###### Python files are initiated using the python command. The help instruction is called using -h option. E.g:
```
python Iris_to_Dataset.py [-h] [-p PATH] [-o OUTPUTFILE]
```

###### Bin files are intiated by specifying the path to the file. E.g, if within the files' directory:

```
./Iris_to_Dataset [-h] [-p PATH] [-o OUTPUTFILE]
```

`Colony Size is stated as the phenotype within the below examples for ease. However, any Iris phenotype (e.g opacity, circularity etc) can be analysed`

#### iris_to_dataset <a name="Iris_to_dataset"></a>

Takes a directory of Iris files and turns them into the combined .csv dataset used for normalisation. 

Input files do not have to be from IRIS as long as they are in the IRIS format and named in the format below.

The IRIS file format is a tab delimited tabular file starting with columns for plate locations, followed by measured phenotypes. E.g:

| row |	column|	size | circularity | opacity |
|  ------------- | ------------- | ------------- | ------------- | :-------------: |
| 1 |	1	| 12348	 | 0.549 | 512559 |
| 1 |	2	| 11786	 | 0.572 | 509877 |
| 1 |	3	| 11265	 | 0.578 | 488846 |

Iris file names MUST be in the format: 

`CONDITION-concentration-platenumber-batchnumber_replicate.JPG.iris`

E.g. `AMPICILLIN-50 mM-6-1_B.JPG.iris`

Where concentrations have decimals, use a comma instead of a period:

E.g. `AMPICILLIN-0,5 mM-6-1_B.JPG.iris`

Where a concentration is not relevant, put two dashes between condition and plate number:

E.g. `LB--1-2_A.JPG.iris`

If only one source plate and/or only one batch was produced, assign 1 for these:

E.g. `AMPICILLIN-0,5 mM-1-1_B.JPG.iris`

`platenumber` refers to the source plate number, i.e which mutants are on the experiment plate. This will match the plate information file number in later steps.


```
usage: iris_to_dataset [-h] [-p PATH] [-o OUTPUTFILE]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -p PATH, --PATH PATH  
                        Path to folder which contains IRIS files, IRIS file names should be in the format: CONDITION-concentration-platenumber-batchnumber_replicate.JPG.iris (default:None)

  -o OUTPUTFILE, --outputfile OUTPUTFILE
                        Name of output file, should be a .csv file (default:None)

  -it IRISPHENOTYPE, --IRISphenotype IRISPHENOTYPE
                        IRIS phenotype from IRIS files you wish to analyse (default: size)

```
    
#### check_normalisation <a name="Check_Normalisation"></a>

Checks each plate individually to see if outer-edge normalisation is required due to plate effects. The module uses the wilcoxson rank sum test to determine if the distribution of outer edge colony sizes, e.g colony size, are the same as the inner colony sizes. If the distributions differ, the outer edge is normalised such that the row or column median of each outer edge colony is equal to the Plate Middle Mean (PMM). The PMM is equal to the mean colony size of all colonies within the middle of the plate within the 40th to 60th percentile of size. Following this, all plates are normalised such that all colonies are scaled to adjust the PMM to the medain colony size of all colonies within the dataset. 

False zero values, are also changed to NaNs, false zero values are values where a colony has a size of zero but other replicates within the condition are not. This is likely due to pinning defects. 

```
usage: check_normalisation [-h] [-i INPUTFILE] [-o OUTPUTFILE]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        The txt file from Iris_to_dataset.py of the colony dataset with conditions across the top and row/column coordinates downwards (default: None)
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        CSV file of the normalised colony sizes. (default: None)
```

#### z_score <a name="Z_Score"></a>

Compares each replicate colony to find outliers within colony size for each plate. Outliers include, colonies smaller than the mean of the replicates (S), colonies bigger than the mean of the replicates (B) and NaN values (X). 

```
usage: z_score [-h] [-i INPUTFILE] [-o OUTPUTFILE]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        The normalised csv file from check_normalisation (default: None)
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        A CSV file of the dataset where colony sizes are replaced with the colony type values. (default: None)

```

#### z_score_count <a name="Z_Score_Count"></a>

Counts the number of each colony type within each plate and the percentage of each colony type. 

```
usage: z_score_count [-h] [-i INPUTFILE] [-o OUTPUTFILE]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        The CSV output file from z_score (default: None)
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        A CSV file with the counts and percentages of the each colony type. (default: None)
```

#### mw_plate_level <a name="Mann_Whitney_Plate_Level"></a>

Compares the distributions of the colony sizes of replicate plates of the same condition and determines if replicate plates have the same distribution based on the p value of the Mann whitney test. A p-value < ⍺ indicates that the two distributions differ with statistical signifcance. The mean p-value is then averaged for each replicate, e.g average(A vs B, A vs C, A vs D) = replicate mean of A. 

```
usage: mw_plate_level [-h] [-i INPUTFILE] [-o OUTPUTFILE] [-o2 OUTPUTFILE_MEAN]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        The normalised csv file from check_normalisation output (default:None)
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        A CSV file with the u statistics and p-values for each comparison. (default: None)
  -o2 OUTPUTFILE_MEAN, --OutputFile_Mean OUTPUTFILE_MEAN
                        A CSV file with the mean u statistics and p-values for each replicate (default: None)

``` 

#### mw_condition_level <a name="Mann_Whitney_Condition_Level"></a>

The variance of the replicate means for each condition is calculated and then the average of these variance is calculated for each plate within that conditions, i.e the variance between replicate plate A,B,C,D for plate 1 of condition A, and then the average of plate 1, 2, 3 etc. for condition A. 

```
usage: mw_condition_level [-h] [-i INPUTFILE] [-o OUTPUTFILE]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        The CSV file with the mean u statistics and p values for each replicate from mw_plate_level (default: None)
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        A CSV file with the mean variance values for the u statistic and p values of the mann-whitney test for each condition. (default: None)
```
#### condition_variance <a name="Condition_Variance"></a>

The variance of replicate colony sizes is calculated for each plate and these variance values are averaged for each plate within a condition. 

```
usage: condition_variance [-h] [-i INPUTFILE] [-o OUTPUTFILE]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        The normalised csv file from check_normalisation (default: None)
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        A CSV file of the average variances for each condition. (default: None)
```
#### pass_fail_conditions <a name="Pass_Fail_Conditions"></a>

The output files of the Mann-Whitney condition level analysis and the condition variance analysis are inputted. The files are tested to see which conditions fail at certain thresholds of variance and Mann-Whitney p value. 

```
usage: pass_fail_conditions [-h] [-iv INPUTFILE_VARIANCE] [-imwc INPUTFILE_MWC] [-ov OUTPUTFILE_VARIANCE] [-omwc OUTPUTFILE_MWC]

optional arguments:
  -h, --help            show this help message and exit
  -iv INPUTFILE_VARIANCE, --InputFile_Variance INPUTFILE_VARIANCE
                        Output file from condition_variance (default: None)
  -imwc INPUTFILE_MWC, --InputFile_MWC INPUTFILE_MWC
                        Output file from mw_condition_level (default: None)
  -ov OUTPUTFILE_VARIANCE, --OutputFile_Variance OUTPUTFILE_VARIANCE
                        A CSV file showing the conditions and the thresholds at which they pass and fail. Here variances which are greater than the threshold tested fail. (default: None)
  -omwc OUTPUTFILE_MWC, --OutputFile_MWC OUTPUTFILE_MWC
                        A CSV file showing the conditions and the thresholds at which they pass and fail. Here p values which are lower than the threshold tested fail. (default: None)
```

#### pass_fail_plates <a name="Pass_Fail_Plates"></a>

The output files of the Mann-Whitney plate level analysis and the Z score analysis are inputted. The files are tested to see which conditions fail at certain thresholds of Normality and Mann-Whitney p value.

```
usage: pass_fail_plates [-h] [-iz INPUTFILE_Z_SCORE] [-imwp INPUTFILE_MWP] [-oz OUTPUTFILE_Z_SCORE] [-omwp OUTPUTFILE_MWP] [-mo MERGED_OUTPUTFILE]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -iz INPUTFILE_Z_SCORE, --InputFile_Z_Score INPUTFILE_Z_SCORE
                        output file from z_score_count (default: None)
  -imwp INPUTFILE_MWP, --InputFile_MWP INPUTFILE_MWP
                        output file from mw_plate_level (default: None)
  -oz OUTPUTFILE_Z_SCORE, --OutputFile_Z_Score OUTPUTFILE_Z_SCORE
                        A CSV file showing the plates and the thresholds at which they pass and fail for the Z-score test. Here normality percentages which are lower than the threshold tested fail. (default: None)
  -omwp OUTPUTFILE_MWP, --OutputFile_MWP OUTPUTFILE_MWP
                        A CSV file showing the plates and the thresholds at which they pass and fail for the Mann-Whitney test. Here p values which are lower than the threshold tested fail. (default: None)
  -mo MERGED_OUTPUTFILE, --Merged_Outputfile MERGED_OUTPUTFILE 
                        A CSV file showing the plates and the thresholds at which they pass and fail for both. (default: None)
```

#### bar_plot_plates <a name="Bar_Plot_Plate"></a>

Produces a bar plot showing the counts of conditions with a certain number of plates lost at different thresholds of normality (z-score) and Mann-Whitney p-value.

```
usage: bar_plot_plates [-h] [-i INPUTFILE] [-o OUTPUTPLOT]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        The input file for the module. Uses the merged output file from pass_fail_plates (default: None)
  -o OUTPUTPLOT, --OutputPlot OUTPUTPLOT
                        Name of output file, a PDF of the bar chart (default: None)

```

#### bar_plot_conditions <a name="Bar_Plot_Conditions"></a>

Produces a bar plot showing the counts of conditions with a certain number of plates lost at different thresholds of Variance and Mann-Whitney mean p-value variance .

```
usage: bar_plot_conditions [-h] [-i INPUTFILE] [-o OUTPUTPLOT]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        The input file for the module. Uses the output files from pass_fail_conditions, either Variance or Mann_Whitney. (default: None)
  -o OUTPUTPLOT, --OutputPlot OUTPUTPLOT
                        Name of output file, a PDF of the bar chart (default: None)
```

#### mw_plates_to_remove <a name="MW_Plates_to_Remove"></a>

Outputs a list of plates which were removed at a certain chosen threshold for the Mann-Whitney test. Also outputs a new dataset to go back into the process of normalisation and scoring, but with detrimental plates removed. 

```
usage: mw_plates_to_remove [-h] [-i INPUTFILE] [-o OUTPUTFILE] [-od ORIGINAL_DATASET] [-or OUTPUT_REMOVED] [-t THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        Input file is the mean output from mw_plate_level (default: None)
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        A CSV file with the name of the plates that were removed and their file names. (default: None)
  -od ORIGINAL_DATASET, --Original_Dataset ORIGINAL_DATASET
                        The original .csv dataset used in the first stage or the output of z_plates_to_remove to remove more plates (default: None)
  -or OUTPUT_REMOVED, --Output_removed OUTPUT_REMOVED
                        A .csv dataset with detrimental plates removed. (default: None)
  -t THRESHOLD, --Threshold THRESHOLD
                        A chosen threshold, usually based off of the bar chart produced by bar_plot_plates. (default: None)
```

#### z_plates_to_remove <a name="Z_Plates_to_Remove"></a>
    
Outputs a list of plates which were removed at a certain chosen threshold for the Z-score test. Also outputs a new dataset to go back into the process of normalisation and scoring, but with detrimental plates removed. 

```
usage: z_plates_to_remove [-h] [-i INPUTFILE] [-o OUTPUTFILE] [-od ORIGINAL_DATASET] [-or OUTPUT_REMOVED] [-t THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        output from z_score_count (default: None)
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        A CSV file with the name of the plates that were removed and their file names. (default: None)
  -od ORIGINAL_DATASET, --Original_Dataset ORIGINAL_DATASET
                        The original .csv dataset used in the first stage or the output of mw_plates_to_remove to remove more plates (default: None)
  -or OUTPUT_REMOVED, --Output_removed OUTPUT_REMOVED
                        A .csv dataset with detrimental plates removed. (default: None)
  -t THRESHOLD, --Threshold THRESHOLD
                        A chosen threshold, usually based off of the bar chart produced by bar_plot_plates. (default: None)
```

#### mw_conditions_to_remove <a name="MW_Conditions_to_Remove"></a>

Outputs a list of conditions which were removed at a certain chosen threshold for the Mann Whitney Condition Level test. Also outputs a new dataset to go back into the process of normalisation and scoring, but with detrimental plates removed.

```
usage: mw_conditions_to_remove [-h] [-i INPUTFILE] [-o OUTPUTFILE] [-od ORIGINAL_DATASET] [-or OUTPUT_REMOVED] [-t THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        output from mw_condition_level (default: None)
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        A CSV file with the name of the plates that were removed and their file names. (default: None)
  -od ORIGINAL_DATASET, --Original_Dataset ORIGINAL_DATASET
                        The original .csv dataset used in the first stage or the output of mw_plates_to_remove or z_plates_to_remove or variance_conditions_to_remove to remove more plates (default: None)
  -or OUTPUT_REMOVED, --Output_removed OUTPUT_REMOVED
                        A .csv dataset with detrimental plates removed. (default: None)
  -t THRESHOLD, --Threshold THRESHOLD
                        A chosen threshold, usually based off of the bar chart produced by Bar_plot_Condition.py. (default: None)
```


#### variance_conditions_to_remove <a name="Variance_Conditions_to_Remove"></a>

Outputs a list of conditions which were removed at a certain chosen threshold for the variance test. Also outputs a new dataset to go back into the process of normalisation and scoring, but with detrimental plates removed.

```
usage: variance_conditions_to_remove [-h] [-i INPUTFILE] [-o OUTPUTFILE] [-od ORIGINAL_DATASET] [-or OUTPUT_REMOVED] [-t THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        output from condition_variance (default: None)
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        A CSV file with the name of the plates that were removed and their file names. (default: None)
  -od ORIGINAL_DATASET, --Original_Dataset ORIGINAL_DATASET
                        The original .csv dataset used in the first stage or the output of mw_plates_to_remove or z_plates_to_remove to remove more plates (default: None)
  -or OUTPUT_REMOVED, --Output_removed OUTPUT_REMOVED
                        A .csv dataset with detrimental plates removed. (default: None)
  -t THRESHOLD, --Threshold THRESHOLD
                        A chosen threshold, usually based off of the bar chart produced by Bar_plot_Condition.py. (default: None)

```

#### s_scores <a name="S_Scores"></a>

Computes the S-scores from the normalised dataset. 

```
usage: s_scores [-h] [-i INPUTFILE] [-o OUTPUTFILE]

Computes the S-scores from the normalised dataset.

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        The normalised csv file from check_normalisation (default: None)
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        A CSV file of the dataset as S-scores (default: None)
```

#### add_gene_names <a name="Add_Gene_Names"></a>

Add the gene names from the plate info files to make the final dataset. 

The plate info files must be in a folder by themselves and should be .txt files. Files such also be numbered e.g:

```
📂Plate_info
 ┣ 📜plat1.txt
 ┗ 📜plat2.txt
```

Plate info files should be formatted as such:

| Row	          | Column      | Strain | 
| ------------- |:-------------:| :-----: |
| 1          | 1 |  PA1230   |
| 1        | 2 |  PA2543   |

```
usage: add_gene_names [-h] [-i INPUTFILE] [-o OUTPUTFILE] [-p PATH]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        The CSV output of s_scores (default: None)
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        A CSV file of the final dataset. (default: None)
  -p PATH, --PATH PATH  
                        The path to the folder containing the plate info files. (default: None)
```

#### cosine_similarity <a name="Cosine_Similarity"></a>

Calculates the cosine similarity scores for the phenotypic profiles of genes from the same operon and genes from different operons. Produces a density plot of the cosine similarity scores for genes of the same and different operons. Produces an ROC curve testing models ability at different threshold. 

```
usage: cosine_similarity [-h] [-i INPUTFILE] [-o OUTPUTFILE] [-or OUTPUT_ROC_CURVE] [-od OUTPUT_DENSITY_PLOT] [-clus CLUSTER_FILE]

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -i INPUTFILE, --InputFile INPUTFILE
                        The dataset with gene names added. Output from add_gene_names (default: None)
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        List of genes compared and the cosine similarity score as well as if they belong to the same operon (default: None)
  -or OUTPUT_ROC_CURVE, --Output_ROC_curve OUTPUT_ROC_CURVE
                        Plot of the ROC curve and AUC score. (default: None)
  -od OUTPUT_DENSITY_PLOT, --Output_Density_plot OUTPUT_DENSITY_PLOT
                        Density plot of the cosine similarity scores for same and different operons. (default: None)
  -clus CLUSTER_FILE, --Cluster_file CLUSTER_FILE
                        A CSV file containing the operon clusters for each gene within the bacterium of interest, where columns = (Cluster,Gene). (default: None)
```

#### chemgapp_small <a name="ChemGAPP_Small"></a>

ChemGAPP Small is an extension within ChemGAPP for the analysis of small scale chemical genomic screens. ChemGAPP Small produces three types of plots, a heatmap, bar plots and swarm plots. For the bar plot and heatmap, ChemGAPP Small compares the mean colony size of within plate replicates to the mean colony size of the within plate wildtype replicates, producing a fitness ratio. The bar plots are then optionally grouped by strain or by condition. The heatmap displays all conditions and strains. For the swarm plots each mutant colony size is divided by the mean colony size of the wildtype, to produce the fitness ratio. A one-way ANOVA and Tukey-HSD analysis determines the significance in difference between each mutant fitness ratio distribution and the wildtype fitness ratio distribution.

Ensure IRIS file names are in the format: `CONDITION-concentration-platenumber_replicate.JPG.iris`

E.g. `AMPICILLIN-50mM-6_B.JPG.iris`

Where concentrations have decimals, use a comma instead of a period:

E.g. `AMPICILLIN-0,5mM-6_B.JPG.iris`

```
usage: chemgapp_small [-h] [-p PATH] [-o OUTPUTFILE_PREFIX] [-pf PLATEINFOPATH] [-m MAX_COLONY_SIZE] [-wt WILDTYPE] [-it IRIS_TYPE] [-col_plot COLOURPALETTE] [-col_heat COLOURHEATMAP] [-wd WIDTH] [-ht HEIGHT] [-r ROTATION] [-cs CIRCLESIZE] [-g GROUP] [-pt PLOTTYPE]
                         [-rm REMOVE_STRAIN] [-ymax Y_MAX] [-ymin Y_MIN]

Analyses small scale chemical genomic screen data

optional arguments:
  -h, --help            show this help message and exit
  -p PATH, --PATH PATH  Path to folder which contains IRIS files (default: None)
  -o OUTPUTFILE_PREFIX, --outputfile_prefix OUTPUTFILE_PREFIX
                        Path and prefix for output file (default: None)
  -pf PLATEINFOPATH, --PlateInfoPath PLATEINFOPATH
                        The path to the folder containing the plate info files. (default: None)
  -m MAX_COLONY_SIZE, --max_colony_size MAX_COLONY_SIZE
                        Maximum colony size allowed, any colony larger than this will be set to this maximum (default: False)
  -wt WILDTYPE, --WildType WILDTYPE
                        Name of wild type strain within plate info file. (default: None)
  -it IRIS_TYPE, --IRIS_type IRIS_TYPE
                        Input IRIS morphology to test. Options: size,circularity,opacity (default: size)
  -col_plot COLOURPALETTE, --colourpalette COLOURPALETTE
                        Name of Seaborn colour palette to use for the bar and swarm plots. (default: GnBu)
  -col_heat COLOURHEATMAP, --colourheatmap COLOURHEATMAP
                        Name of Seaborn colour palette to use for the heatmap. (default: bwr_r)
  -wd WIDTH, --width WIDTH
                        Figure width to use for the graphs. (default: 5)
  -ht HEIGHT, --height HEIGHT
                        Figure height to use for the graphs. (default: 5)
  -r ROTATION, --rotation ROTATION
                        X Axis label rotation (default: 90)
  -cs CIRCLESIZE, --CircleSize CIRCLESIZE
                        SwarmPlot circle size (default: 2.5)
  -g GROUP, --group GROUP
                        Group bar plots by strain or condition. Options = strain, condition. (default: condition)
  -pt PLOTTYPE, --PlotType PLOTTYPE
                        Type of Plot. Options: barplot, swarmplot (default: barplot)
  -rm REMOVE_STRAIN, --remove_strain REMOVE_STRAIN
                        txt file of strain names to remove separated by ';'. Names must match those in plate information file. E.g. mutant1;mutant2;mutant4 (default: None)
  -ymax Y_MAX, --y_max Y_MAX
                        Maximum limit for y axis (default: None)
  -ymin Y_MIN, --y_min Y_MIN
                        Minimum limit for y axis (default: None)

```     

#### gi_dataset <a name="GI_Dataset"></a>

gi_dataset calculates the fitness ratio of two single mutant strains and a double knockout in comparison to the wildtype. It also calculates the expected double knockout fitness ratio for comparison to the observed fitness ratio. This outputs a Colony_Size.csv file and Interaction_Score.csv file for each secondary gene within the pair.

Ensure IRIS file names are in the format: `SecondaryGeneName_replicate.JPG.iris`

E.g. `MexB_A.JPG.iris`

```

usage: gi_dataset [-h] [-i INPUTFILE] [-p PATH] [-n NAMEINFOFILE]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputfile INPUTFILE
                        Input IRIS file. (default: None)
  -p PATH, --PATH PATH  Path for the output files. (default: None)
  -n NAMEINFOFILE, --nameinfofile NAMEINFOFILE
                        The plate information file. Plate info files should be txt files, with the columns: Row, Column, Strain, Replicate, Order, Set. (default: None)

```

#### gi_barplot <a name="GI_Barplot"></a>

GI_Barplot produces a grouped bar plot for all Interaction_score.csv files within the given directory.  

```

usage: gi_barplot [-h] [-p PATH] [-o OUTPUTFILE] [-g PRIMARYGENE]

optional arguments:
  -h, --help            show this help message and exit
  -p PATH, --PATH PATH  Path to Interaction Score Files (default: None)
  -o OUTPUTFILE, --OutputFile OUTPUTFILE
                        A PDF file of the final bar plot. (default: None)
  -g PRIMARYGENE, --PrimaryGene PRIMARYGENE
                        The primary interacting gene being compared (default: None)

```

### ChemGAPP Big

#### Step_1_Normalisation.py <a name="Step_1_Normalisation"></a>

1- Upload your IRIS files.

Ensure IRIS file names are in the format:

`CONDITION-concentration-platenumber-batchnumber_replicate.JPG.iris`

E.g. `AMPICILLIN-50 mM-6-1_B.JPG.iris`

Where concentrations have decimals, use a comma instead of a period:

E.g. `AMPICILLIN-0,5 mM-6-1_B.JPG.iris`

Where a concentration is not relevant, put two dashes between condition and plate number:

E.g. `LB--1-2_A.JPG.iris`

If only one source plate and/or only one batch was produced, assign 1 for these:

E.g. `AMPICILLIN-0,5 mM-1-1_B.JPG.iris`

`platenumber` refers to the source plate number, i.e which mutants are on the experiment plate. This will match the plate information file number in later steps.


2- Enter a path to the folder you would like to save the output files to. Ensure you include a prefix which will be added to the start of all output file names e.g: ~/Projects/Test


3- Input the IRIS phenotype to analyse. The spelling should exactly match the IRIS file spelling and capitalisations. 


4- Press `Begin!`

#### Step_2_Threshold_Selector.py <a name="Step_2_Threshold_Selector"></a>

1- Type the threshold values into the corresponding boxes based on the bar plots below. 


2- Then select which statistical tests you would like to use to remove detrimental data. *Any combination can be selected.* 


3- Press `Begin Quality Control Tests!`


#### Step_3_S_Score_Calculator.py <a name="Step_3_S_Score_Calculator"></a>


1- Upload plate information files. 

2- Select if you want to score the original dataset, the curated dataset or both.

#### Step_4_Dataset_Comparison.py <a name="Step_4_Dataset_Comparison"></a>

1- Simply upload your cluster file. 

File should be in CSV format and consist of two columns; Cluster and Gene.


E.g: 

| Cluster | Gene |
|  ------------- | :-------------: |
| 1 | PA14_00050 |
| 2 | PA14_00060 |
| 2 | PA14_00070 |
| 4 | PA14_00080 |
| 5 | PA14_00090 |

### ChemGAPP Small

#### Step_1_chemgapp_small <a name="Step_1_ChemGAPP_Small"></a>

1- First upload all iris files you wish to include.
Ensure IRIS file names are in the format: `CONDITION-concentration-platenumber_replicate.JPG.iris`

E.g. `AMPICILLIN-50mM-6_B.JPG.iris`

Where concentrations have decimals, use a comma instead of a period:

E.g. `AMPICILLIN-0,5mM-6_B.JPG.iris`


2- Enter a path to the folder you would like to save the output files to. Ensure you include a prefix which will be added to the start of all output file names e.g: ~/Desktop/project/Test

Output files include:

`Test_Intial_dataset.csv`

`Test_Normalised_dataset.csv`

`Test_Scored_Dataset.csv`

`Test_Final_dataset.csv`


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


4- Enter the name of the wild type strain. 

This must match the name given in the plate information file.
E.g. `WT`


5- Optionally enter the name of a Seaborn colour palette to change the colour of the output plots. 

The default palette is `Spectral`.

6- Select how you want the bar plots to be group. Either by strain or by condition.

7- Press Begin!


8- To save bar plot images, right click on the image and select `Save as Image`. This will save the image as a png file.

### ChemGAPP GI

#### Step_1_Interaction_Scores.py <a name="Step_1_Interaction_Scores"></a>

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

2. Upload your IRIS Files.

Ensure IRIS file names are in the format: `SecondaryGeneName_replicate.JPG.iris`

E.g. `MexB_A.JPG.iris`



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
        


4. Press Begin. 

#### Step_2_Bar_Plot.py <a name="Step_2_Bar_Plot"></a>

1. Upload the desired interaction score files. 
These are the X_Interaction_Scores.csv files from the previous step.
Ensure the file names end with 'Interaction_Scores.csv'. 

```
E.g.

📂Interaction_Score_Files
 ┣ 📜 A_Interaction_Scores.csv
 ┣ 📜 B_Interaction_Scores.csv
 ┗ 📜 C_Interaction_Scores.csv
 
```


2. Input the name of the Primary Gene. This is the gene common to all the tested gene pairs.
    
E.g.

MexA::MexB
MexA::OmpM
MexA::MexY

MexA would be the Primary Gene



### Test files <a name="package"></a>

These files are located in the `ChemGAPP/Test_Files` folder and include:

1.  `ChemGAPP_Big/` Folder containing test IRIS files, plate information files and cluster data for ChemGAPP Big. 

2. `ChemGAPP_Small/`  Folder containing test IRIS files and plate information file ChemGAPP Small.

3. `ChemGAPP_GI/` Folder containing test IRIS files and plate information file ChemGAPP GI.

----
### Contact <a name="contact"></a>
For queries, please contact [Hannah Doherty](mailto:hxd476@student.bham.ac.uk?subject=[GitHub:ChemGAPP]), Institute of Microbiology and Infection, School of Biosciences, University of Birmingham. In collaboration with Center for Computational Biology, University of Birmingham.  

----