#!/usr/bin/env python
# coding: utf-8

# In[ ]:
import argparse
parser = argparse.ArgumentParser(description="Produces a bar plot showing the counts of conditions with a certain number of plates lost at different thresholds of Variance and Mann-Whitney mean p-value variance.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--InputFile", help="The input file for the module. Uses the output files from Pass_Fail_Conditions.py, either Variance or Mann_Whitney.")
parser.add_argument("-o", "--OutputPlot", help="Name of output file, a PDF of the bar chart")
args = vars(parser.parse_args())
passfailfile = args["InputFile"]
outputplot = args["OutputPlot"]

def Bar_plot_Condition(pass_fail_file,otptplot):
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    Pass_Fail = pd.read_csv(pass_fail_file)
    Pass_Fail = Pass_Fail.set_index(['Condition','Batch'])
    bar = pd.DataFrame(columns = ['Condition','Batch', 'Threshold','Percentage'])
    for p5 in sorted(Pass_Fail.index):
        dfpf14 = Pass_Fail.xs((p5), axis=0, drop_level=False)
        dfpf14 = pd.DataFrame(dfpf14)
        dfpf14 = pd.DataFrame.transpose(dfpf14)
        for cols in dfpf14.columns:
            if cols[0] == 'M' or cols[0] == 'N' or cols[0] == 'A':
                failed = ((len(dfpf14[dfpf14[cols] == 'F'])/dfpf14[cols].count())*100)
                name=(dfpf14.index[0][0],dfpf14.index[0][1],cols,failed)
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
    b = sns.displot(bar, x='Threshold', hue='Percentage',
                    multiple='stack', legend = False, palette=("White","orange"))
    b.set_axis_labels("", "Count of Conditions Lost")
    b.set_xticklabels(rotation=90)
    plt.title("Count of Conditions lost at Different Thresholds")
    plt.savefig(otptplot, bbox_inches='tight')
    return bar

Bar_plot_Condition(passfailfile,outputplot)