#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def get_options():
    parser = argparse.ArgumentParser(description="Produces a bar plot showing the counts of conditions with a certain number of plates lost at different thresholds of normality (z-score) and Mann-Whitney p-value.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--InputFile", help="The input file for the module. Uses the merged output file from Pass_Fail_Plates.py")
    parser.add_argument("-o", "--OutputPlot", help="Name of output file, a PDF of the bar chart")
    return parser.parse_args()

def main():
    options = get_options()
    pass_fail_file = options.InputFile
    otptplot = options.OutputPlot
    Pass_Fail = pd.read_csv(pass_fail_file)
    Pass_Fail['Plate'] = pd.to_numeric(Pass_Fail['Plate'])
    Pass_Fail = Pass_Fail.set_index(['Condition','Plate'])
    plt5 = {x[0:2] for x in Pass_Fail.index}
    bar = pd.DataFrame(columns = ['Condition','Plate','Batch','Replicates', 'Threshold','Percentage'])
    for p5 in sorted(plt5):
        dfpf14 = Pass_Fail.xs((p5), axis=0, drop_level=False)
        #splits the data by condition and source plate number so number of F's 
        # can be counted for each threshold across replicates to produce percentage
        # fails by dividing by overall count and multiplying by 100. This then allows
        # the plotting of bar charts with percentage of plates within a condition failing at certain thresholds
        for cols in dfpf14.columns:
            if cols[0] == 'M' or cols[0] == 'N':
                failed = ((len(dfpf14[dfpf14[cols] == 'F'])/dfpf14[cols].count())*100)
                name=(dfpf14.index[0][0],dfpf14.index[0][1],dfpf14['Batch'][0],dfpf14['Replicates within Condition/Plate'][0],cols,failed)
                columns = list(bar)
                data = []
                zipped = zip(columns, name)
                a_dictionary = dict(zipped)
                data.append(a_dictionary)
                bar = bar.append(data, True) 
        bar = bar.round(2)
    bar_list = bar['Percentage'].astype(float).values
    b3 = {x for x in bar_list}
    b3 = list(sorted(b3))
    bar['Percentage'] = pd.Categorical(bar['Percentage'], categories=b3,ordered=True)
    sns.set_theme(style="whitegrid")
    #stacks the percentages so you can see the how much each percentage contributes at different thresholds.
    b = sns.displot(bar, x='Threshold', hue='Percentage',
                    multiple='stack', legend = True, palette=sns.color_palette("Spectral_r",len(b3)))
    b.set_axis_labels("Threshold Types", "Count of Condition With % Plates Lost")
    b.set_xticklabels(rotation=90)
    plt.title("Count of Conditions with Different Percentages of Plates Lost at Different Thresholds")
    plt.savefig(otptplot, bbox_inches='tight')
    return dfpf14

if __name__ == "__main__":
    main()