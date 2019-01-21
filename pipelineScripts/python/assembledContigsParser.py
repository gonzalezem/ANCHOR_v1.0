#!/usr/bin/python
from datetime import datetime
import argparse
import pandas as pd
from pandas import *
from pandas.tools.plotting import table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns





def convert(x):
     try:
         return x.astype(int)
     except:
         return x



parser = argparse.ArgumentParser(description="Script creates a histogram of contig (i.e. merge reads) lengths and a table about average merge lengths. It uses as inpufile report file from Mothur s make.contigs (16S pipeline).",\
    epilog='2 files will be created: xx_Merge_Overview.pdf xx_ReadsMerge_Average_Stats.txt  .')
parser.add_argument('-r', '--report', help='[REQUIRED] Output report from make.contigs. Example: overlap_analysis.contigs.report. ', dest='input_report', action='store', required=True)

args = parser.parse_args()

report = args.input_report
prefix = report.split('.')[0]


print "\n-----------------------------------------"
print "Creating main file from the input files"
df = pd.read_csv(report, sep = '\t')

df['percentNs'] = df['Num_Ns'] / df['Overlap_Length'] * 100
df['percentMisMatches'] = df['MisMatches'] / df['Overlap_Length'] * 100

#Creating the output for all hits
df_all = df.copy()
df_all['Number_of_Seq'] = df_all.groupby(['Length'])['Length'].transform('count')
df_all = df_all.sort_values(by=['Number_of_Seq', 'Length'], ascending=[False, False])
print "done"                                        
print "-----------------------------------------\n"



print "\n-----------------------------------------"
print "Calculating Statistics for all data"
print "Note: only mean values will be reported in the output file"
df = df[['Length', 'Overlap_Length', 'MisMatches', 'Num_Ns', 'percentNs', 'percentMisMatches']]
#Number of elements with a particular overlap length
df['Number_of_Seq'] = df.groupby(['Length'])['Length'].transform('count')
df = df.sort_values(by=['Number_of_Seq', 'Length'], ascending=[False, False])
grouped = df.groupby(['Length'], sort=False)
stats = dict()
i=0
for name, group in grouped:
    i=i+1
    stats[i] = group.describe()
    #kepping only the mean statistics
    stats[i] = stats[i].drop(['std', 'count', 'min', 'max', '25%', '50%', '75%'])
    stats[i]
    if i == 1:
        #create the dataframe that we will output
        table_stats = pd.DataFrame.from_dict(stats[i])
    if i > 1:
        #incrementing the dataframe
        table_stats = pd.concat([stats[i], table_stats], axis=0)

table_stats = table_stats.sort_values(by=['Number_of_Seq', 'Length'], ascending=[False, False]).reset_index(drop = True)
#Converting table to integers:
table_stats = table_stats.apply(convert)
print "done"                                        
print "-----------------------------------------\n"



print "\n-----------------------------------------"
print "Creating overview of the read overlaps for all amplicons"
# Set up the matplotlib figure
f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(30, 10), sharex=False)
#Title
plt.suptitle('Overview of the read overlaps')

#plotting the amplicon distribution
sns.barplot(df['Length'], df['Number_of_Seq'], palette="Reds_r", ax=ax1)
ax1.set_ylabel("Number of amplicons")
ax1.set_xlabel("")
for item in ax1.get_xticklabels():
    item.set_rotation(90)

# Percentage of mismatches
sns.barplot(df['Length'], df['percentMisMatches'], palette="Blues_d", ax=ax2)
ax2.set_ylabel("Mismatches %")
ax2.set_xlabel("")
for item in ax2.get_xticklabels():
    item.set_rotation(90)

# Percentages of Ns
sns.barplot(df['Length'], df['percentNs'], palette="Greens_d", ax=ax3)
ax3.set_ylabel("Ambiguous bases %")
ax3.set_xlabel("Amplicon length (bp)")
for item in ax3.get_xticklabels():
    item.set_rotation(90)

# Finalize the plot
sns.despine(bottom=True)
#plt.setp(f.axes, yticks=[])
plt.tight_layout(h_pad=3)

#output
output_filename= prefix + "_Merge_Overview.pdf"
resolution = 600
plt.savefig(output_filename, bbox_inches='tight', dpi=resolution)
print "done"                                        
print "-----------------------------------------\n"


print "\nWriting files to disk"
outputfileAll = prefix + '_seq_all.txt'
#outputfileremove = prefix + '_seq_to_remove.txt'
print outputfileAll
df_all.to_csv(outputfileAll, sep = '\t', index = False)
#print outputfileremove
#dfremoved.to_csv(outputfileremove, sep = '\t', index = False)
outputfilestats = prefix + '_ReadsMerge_Average_Stats.txt'
print outputfilestats
table_stats.to_csv(outputfilestats, sep='\t', index=False, header=True, float_format='%.0f')
print "\n-----\nassembledContigsParser.py exiting normally!"














