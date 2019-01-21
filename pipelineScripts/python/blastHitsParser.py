#!/usr/bin/python
import numpy as np
import argparse
import pandas as pd
from pandas import *
import sys
####import matplotlib.gridspec as gridspec
####import matplotlib.cm as cm, matplotlib.font_manager as fm

print "Pandas version : " + pandas.__version__

startTime = datetime.now()


################################INPUT FILES#################################################
#Output file from blast from part6_parse_blastn_results.pbs.sh
#the columns are the folwing:
#queryid|subjectid|identity|alignmentlength|mismatches|gapopens|qstart|qend|sstart|send|evalue|bitscore|querylength|querydescription|coverageQ|subjectlength|coverageS|subjectdescription
#there is a header
################################################################################################

def HitsPerQuery(df):
    #calculating the number of sequences per query 
    init_number = len(df)
    df = df.drop_duplicates('queryid')
    final_number = len(df)
    return float(init_number) / float(final_number)


startTime = datetime.now()


parser = argparse.ArgumentParser(description="'output file : _output_py")
parser.add_argument('-i', '--inputfile', help='[REQUIRED] file from parseAnchorBlastOutput.sh script', dest='inputfile', action='store', required=True)

args = parser.parse_args()
inputfile = args.inputfile

output = "_output_py"



print "loading input file..."
data = pd.read_csv(inputfile, sep='\t')


#-----------------------------------------------------------------------------------------------------------------------"
#-----------------------------------------------------------------------------------------------------------------------"
#-----------------------------------------------------------------------------------------------------------------------"

#print data.head(n=2)

print '\n'
print "\nINPUT INFORMATION"
print "-------------------------------"
print ('Inputfile in tabular format : ' + str(inputfile))
print "Inputfile : ", "{:,}".format(data.shape[0]), " records"
print "Average number of hits per sequence : ", format(HitsPerQuery(data), "8.2f")
print "-------------------------------\n"



print "\n-----------------------------------------"
print "Keeping highest identity and coverage per sequence"
#I want to keep all % id cov at a given round number (e.g. 99%) and not only the best (99.78% could be the best), so I'll change id and cov to floor
data['Fidentity'] = np.floor(data['identity'])
data['FcoverageQ'] = np.floor(data['coverageQ'])
data['IdentityCov']=data['Fidentity'] + data['FcoverageQ']
#Now I'll be conservative with the output in order to have a maximum good hits: I'll reduce all values >199 down to 199 (our cutoff is 99% id and cov, i.e 198)
data.loc[data.IdentityCov > 199, 'IdentityCov'] = 199
data['MaxIdentityCov'] = data.groupby('queryid')['IdentityCov'].transform(lambda x: x.max())
#Keep the ones with the highest identity
data = data.loc[data['IdentityCov'] == data['MaxIdentityCov']]
#Remove duplicate hits (we'll use the subjectid based on the highest score identity ocoverage)
data = data.sort_values(by = ['queryid', 'IdentityCov']).reset_index(drop=True)
data = data.drop(['MaxIdentityCov', 'IdentityCov', 'Fidentity', 'FcoverageQ'],1)
print "Done! data have", "{:,}".format(data.shape[0]), "records left"
print "Average number of hits per sequence : ", format(HitsPerQuery(data), "8.2f")
print "-----------------------------------------\n"


data.to_csv(output, sep = '\t', index = False)

