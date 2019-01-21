#!/usr/bin/python
import numpy as np
from datetime import datetime
import argparse
import pandas as pd
from pandas import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import locale
import sys
print "Pandas version : " + pandas.__version__

startTime = datetime.now()

def convert(x):
     try:
         return x.astype(int)
     except:
         return x

#Custom Function 
RowCollapse = lambda x:";".join(x.astype(str)) 


def collapseAndSort(mainDF,colName):
    n = 10000  #chunk row size
    list_df = [mainDF[i:i+n] for i in range(0,mainDF.shape[0],n)]
    for j in range (0,len(list_df),1):
        secondaryAmbiSubOTU = list_df[j].groupby("subOTU").agg({'%s'%colName:RowCollapse})
        secondaryAmbiSubOTU = secondaryAmbiSubOTU.reset_index(drop=False)
        #1. Remove duplicates in AmbiguousLabels_2, subjectid_2 and Id/Cov_2
        #Now we will remove duplicated values within AmbiguousLabels
        s = secondaryAmbiSubOTU['%s'%colName].str.split(';').apply(Series, 1).stack()
        s.index = s.index.droplevel(-1) # to line up with data's index
        s.name = '%s'%colName# needs a name to join
        secondaryAmbiSubOTU = secondaryAmbiSubOTU.rename(columns = {'%s'%colName:'%s_OLD'%colName})
        secondaryAmbiSubOTU = secondaryAmbiSubOTU.join(s) #replace with new
        #remove dupplicated values
        secondaryAmbiSubOTU = secondaryAmbiSubOTU.drop_duplicates(["subOTU","%s"%colName])
        #get back to ; separated values in AmbiguousLabels
        secondaryAmbiSubOTU = secondaryAmbiSubOTU.sort_values(by = ['%s'%colName])
        secondaryAmbiSubOTU = secondaryAmbiSubOTU.groupby('subOTU')['%s'%colName].apply(';'.join).to_frame('%s'%colName)
        secondaryAmbiSubOTU = secondaryAmbiSubOTU.reset_index(drop=False)
        #merge the new column to the main data
        list_df[j] = list_df[j].drop(['%s'%colName],1)
        list_df[j] = list_df[j].merge(secondaryAmbiSubOTU, on=['subOTU'], how='inner')
    mainDF = pd.concat(list_df)
    return mainDF


################################INPUT FILES#################################################
#example : subOTUs_allDatabase_withTrueUnknowns_part4_parsed.txt
#inputfile columns:
##Index(['anchorID', 'subOTU', 'anchorCount', 'identity', 'coverageQ',
#       'subjectid_1', 'subjectdescription_1', 'taxonomy', 'taxon1',
#       'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7',
#       'lowestTaxLevel', 'AmbiguousHits', 'AmbiguousLabels_1', 'subjectid_2',
#       'Id/Cov_2', 'AmbiguousLabels_2'],
#      dtype='object')


#example: 98_afterBlast_mappingFile.txt
# highCounters              multipleCountLowCounters        singleCountLowCountSeq
# lab_aft_air_vent_10                   -                   lab_aft_air_vent_10
# lab_aft_air_vent_10                   -                   lab_aft_air_vent_10086
################################################################################################


parser = argparse.ArgumentParser(description="'inputfile: myfile.txt \noutput file : several count files")
parser.add_argument('-i', '--inputfile', help='[REQUIRED] anchor_table.txt (output from anchorParser.sh))', dest='inputfile', action='store', required=True)
parser.add_argument('-c', '--crib', help='[REQUIRED] contig_crib.txt. A crib file to map all sequences to an anchor ID. This inputfile is produced during OTUandCounts.sh', dest='crib', action='store', required=True)

args = parser.parse_args()
inputfile = args.inputfile
cribPath = args.crib

prefix = inputfile.split('_part4_parsed.txt')[0]



print "loading input file..."
data = pd.read_csv(inputfile, sep='\t', quoting=3, doublequote = False, error_bad_lines=False, encoding='utf-8')
crib = pd.read_csv(cribPath, sep='\t', quoting=3, doublequote = False, error_bad_lines=False, encoding='utf-8')

### NEED TO CHANGE UPSTREAM scripts LABEL !!!
crib = crib.rename(columns = {'anchors':'anchorID'})

#-----------------------------------------------------------------------------------------------------------------------"
#-----------------------------------------------------------------------------------------------------------------------"
#-----------------------------------------------------------------------------------------------------------------------"

#print data.head(n=2)

print '\n'
print "\nINPUT INFORMATION"
print "-------------------------------"
print ('Inputfile: ' + str(inputfile))
print ('Crib file: ' + str(cribPath))
print "Inputfile has  ", "{:,}".format(data.shape[0]), " records"
print "-------------------------------\n"



print "\n-----------------------------------------"
print "Parsing the crib file"
#create a column Sample in the crib
crib.loc[:,'sample']=[x.rsplit("_",1)[0] for x in crib['singleCountLowCountSeq']]
crib = crib[['sample','anchorID','singleCountLowCountSeq']]
#create a count column in the crib
crib['AnchorIDcount'] = crib.groupby(['anchorID','sample'])['anchorID'].transform('count')
crib = crib.drop_duplicates(['sample','AnchorIDcount','anchorID'])
crib = crib[['sample','AnchorIDcount','anchorID']]
print "Done!"
print "-----------------------------------------\n"

print "\n-----------------------------------------"
print "Mapping crib to inputfile"
df= pd.merge(data, crib, how='outer', on='anchorID', left_on=None, right_on=None,
                     left_index=False, right_index=False, sort=True,
                     suffixes=('_1', '_2'), indicator=False)

#Now we have the count per AnchorID, but what we want is the count per subOTU, so I'll sum all AnchorIDcount to get it
df['count'] = df.groupby(['sample','subOTU'])['AnchorIDcount'].transform(lambda x: x.sum())
taxa = df[['subOTU','anchorID','identity','coverageQ','count', 'taxonomy', 'sample', 'taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7','lowestTaxLevel','chimera','AmbiguousHits','AmbiguousLabels','subjectid','SubOTU_2','AmbiguousHits_2','AmbiguousLabels_2','subjectid_2','Id/Cov_2','anchorCollapsingDueToPrimer','sequence']].copy()

taxa = taxa.sort_values(by = ['subOTU','identity','coverageQ','taxonomy'], ascending = [False,False,False,True])
taxa = taxa.reset_index(drop=True) 

#now we want: 
#   1. collapse all ambiguous hits for a given subOTU
#   2. Evaluate number of Anchors per SubOTUs
#   2. Choose a representative anchor id for a given subOTU
#1. Collapse AmbiguousLabels_2
taxa = collapseAndSort(taxa,"AmbiguousLabels_2")
taxa = collapseAndSort(taxa,"subjectid_2")
taxa = collapseAndSort(taxa,"Id/Cov_2")
taxa = taxa.drop_duplicates()
#2. Evaluate number of Anchors per SubOTUs
temp = taxa.groupby('subOTU')['anchorID'].apply(lambda x: len(x.unique()))
temp = temp.reset_index(drop=False)
temp = temp.rename(columns = {'anchorID':'anchorPerSubOTU'})
taxa = taxa.merge(temp, on='subOTU', how='inner')
#3. extract representative anchor sequence
taxa = taxa.sort_values(by = ['count','identity','coverageQ'], ascending = [False,False,False])
taxa_subOTU = taxa.drop_duplicates(['subOTU','sample','count'])
#Now we only have one anchor per manOTU, so i'll change the name of the column
taxa_subOTU = taxa_subOTU.rename(columns = {'anchorID':'RefAnchorID'})
print "Done!"
print "-----------------------------------------\n"





print "\n-----------------------------------------"
print "OTU count + taxonomy table"
df=taxa_subOTU.copy()
#renamecolumns
df = df.rename(columns = {'taxon1':'Domain'})
df = df.rename(columns = {'taxon2':'Phylum'})
df = df.rename(columns = {'taxon3':'Class'})
df = df.rename(columns = {'taxon4':'Order'})
df = df.rename(columns = {'taxon5':'Family'})
df = df.rename(columns = {'taxon6':'Genus'})
df = df.rename(columns = {'taxon7':'Species'})
df = df.rename(columns = {'AmbiguousHits':'NumberOfAmbiguousHits'})


#Adding duplicates and removing them
df['Total'] = df.groupby(['sample', 'subOTU'])['count'].transform('sum')
del df['count']
df = df.drop_duplicates(['sample', 'subOTU'])
#Keeping just subOTU taxon
df2 = df[['Total','sample','subOTU']]
grouped = df2.groupby('sample')
temp = []
temp2 = []
i = 0
print "Processing samples:"
for name, group in grouped:
    print name
    if i == 0:
        temp = group.copy()
        temp = temp.rename(columns = {'Total':'%s'%name})
        del temp['sample']
        #print temp.head(n=2)
        temp = temp[['subOTU','%s'%name]]
        #print temp.head()
        temp2 = temp.copy()
        temp2['%s'%name] = temp2['%s'%name].astype(int)
        i = i + 1
        temp = []
    else:
        temp = group.copy()
        temp = temp.rename(columns = {'Total':'%s'%name})
        del temp['sample']
        temp2 = temp2.merge(temp, on='subOTU', how='outer')
        temp2 = temp2.fillna("0")
        temp2['%s'%name] = temp2['%s'%name].astype(int)
        temp = []
temp2 = temp2.apply(convert) #convert all values to integers because some may be strings
temp2['totalcounts'] = temp2.sum(axis=1)
temp2 = temp2.sort_values(by = 'totalcounts', ascending = False)
temp2 = temp2.drop_duplicates()
#Keeping All taxons
df_tax = df[['subOTU','lowestTaxLevel','NumberOfAmbiguousHits','anchorPerSubOTU','RefAnchorID','identity','coverageQ','chimera','AmbiguousLabels','subjectid','SubOTU_2','AmbiguousHits_2','AmbiguousLabels_2','subjectid_2','Id/Cov_2','Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species','anchorCollapsingDueToPrimer','sequence']].copy()
df_tax.loc[df_tax.lowestTaxLevel == 1, 'lowestTaxLevel'] = "D"
df_tax.loc[df_tax.lowestTaxLevel == 2, 'lowestTaxLevel'] = "P"
df_tax.loc[df_tax.lowestTaxLevel == 3, 'lowestTaxLevel'] = "C"
df_tax.loc[df_tax.lowestTaxLevel == 4, 'lowestTaxLevel'] = "O"
df_tax.loc[df_tax.lowestTaxLevel == 5, 'lowestTaxLevel'] = "F"
df_tax.loc[df_tax.lowestTaxLevel == 6, 'lowestTaxLevel'] = "G"
df_tax.loc[df_tax.lowestTaxLevel == 7, 'lowestTaxLevel'] = "S"

df_tax = df_tax.sort_values(by = ['subOTU','identity','coverageQ'], ascending = [True,False,False])
df_tax = df_tax.reset_index(drop=True) 
df_tax = df_tax.drop_duplicates(['subOTU','Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'])
all_taxon_countFile = df_tax.merge(temp2, on='subOTU', how='inner')
all_taxon_countFile = all_taxon_countFile.drop_duplicates()
all_taxon_countFile = all_taxon_countFile.sort_values(by = ['totalcounts'], ascending = [False])
all_taxon_countFile = all_taxon_countFile.reset_index(drop=True) 
all_taxon_countFile = all_taxon_countFile.rename(columns = {'subOTU':'OTU'})
all_taxon_countFile = all_taxon_countFile.rename(columns = {'SubOTU_2':'OTU_2'})
outputfile = 'OTU_table.txt'
all_taxon_countFile.to_csv(outputfile, sep='\t', index=False)
print "Done!"
print "-----------------------------------------\n"



print "\n     OUTPUT INFORMATION:     "
print "-------------------------------"
print ('Main subOTU count file : ' + str(outputfile))
print( '\nOTUandCounts.py was run in: ' + str(datetime.now()-startTime) + ' h:min:ss')

