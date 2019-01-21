#!/usr/bin/python
import numpy as np
from datetime import datetime
import argparse
import pandas as pd
from pandas import *
import locale
import sys
import os
import subprocess

startTime = datetime.now()



def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement


def batch_iterator(iterator, batch_size) :
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

def change_nucl(str, n):
    first_part = str[:n]
    last_part = str[n+1:]
    return first_part + "X" + last_part


def lowestTaxLevelReached(df_ref, taxonLevel_ref, taxonLevel):
    #check if taxon i-1 is the same as taxon i
    index = (df_ref["taxon%s"%taxonLevel_ref] == df_ref["taxon%s"%taxonLevel])
    df = df_ref.ix[index]
    df.is_copy = False
    #Assigning highest taxon level reached for a hit, for each chunk
    df['lowestTaxLevel'] = int(taxonLevel)
    return df

def fakeLastTaxLevelReached(df_ref, taxonLevel_ref, taxonLevel):
    #check if taxon i-1 is the same as taxon i
    index = (df_ref["fakeTaxon%i"%taxonLevel_ref] == df_ref["fakeTaxon%i"%taxonLevel])
    df = df_ref.ix[index]
    df.is_copy = False
    #Assigning highest taxon level reached for a hit, for each chunk
    df.loc[:,'lastFakeTaxon'] = taxonLevel
    return df


def listToString(mycolumn):
    mycolumn = mycolumn.astype(str).str.replace("\'","")
    mycolumn = mycolumn.astype(str).str.replace("'","")
    mycolumn = mycolumn.astype(str).str.replace("\[nan\]","-")
    mycolumn = mycolumn.astype(str).str.replace("\[ nan \]","-")
    mycolumn = mycolumn.astype(str).str.replace("\[ nan \]","-")
    mycolumn = mycolumn.astype(str).str.replace("\[ ","")
    mycolumn = mycolumn.astype(str).str.replace(" \]","")
    mycolumn = mycolumn.astype(str).str.replace(" \[","")
    mycolumn = mycolumn.astype(str).str.replace(" \]","")
    mycolumn = mycolumn.astype(str).str.replace("\[","")
    mycolumn = mycolumn.astype(str).str.replace("\]","")
    mycolumn = mycolumn.astype(str).str.replace(" ",";")
    mycolumn = mycolumn.astype(str).str.replace("nan;","")
    mycolumn = mycolumn.astype(str).str.replace("nan","")
    mycolumn = mycolumn.astype(str).str.replace('\n', '')
    mycolumn = mycolumn.astype(str).str.strip('"')
    return mycolumn 


def cleanTaxonName(data,mycolumn):
   data.loc[:,'%s'%mycolumn]=[x.rsplit("_",1)[0] for x in data.loc[:,'%s'%mycolumn]]
   return data

#Custome function to collapse several rows of a dataframe into 1
RowCollapse = lambda x:";".join(x.astype(str)) 

def sortCollapsedColumn(df,referenceCol,colName):
    df = df.reset_index(drop=True)
    #Sort values inside a collapsed column (ex: from 3;2;1 to 1;2;3) around a reference value (ex: queryid)
    s = df['%s'%colName].str.split(';').apply(Series, 1).stack()
    s.index = s.index.droplevel(-1) # to line up with data's index
    s.name = '%s'%colName# needs a name to join
    df = df.drop(["%s"%colName] ,1)
    df = df.join(s) #replace with new
    #sort 
    df = df.sort_values(by = ["%s"%referenceCol,"%s"%colName])
    #retack subjectid (; separator). That will create a 2 column dataframe (refereceCol as index)
    df2 = df.groupby("%s"%referenceCol)['%s'%colName].apply(';'.join).to_frame('%s'%colName)
    df2 = df2.reset_index(drop=False)
    #merge the new column to the main data
    df = df.drop(['%s'%colName],1)
    df = df.merge(df2, on=["%s"%referenceCol], how='inner')
    df = df.drop_duplicates()
    return df

def workOnOptions(data,optionName):
    data = data.rename(columns = {'taxon7':'subOTU'})

    print "\n-----------------------------------------"
    print "Change Multiple annotation to the last taxon level"
    #For example, for an OTU that is labelled as a order_MS (order being the highest common taxon between all annotation), we will look at the resolution of the last taxon. If it is species, then we will call it MS, if genus, then MG, if family, then MF etc.
    data_no_ambiguity = data.loc[data.AmbiguousHits == 0,]
    data_ambiguity = data.loc[data.AmbiguousHits > 0,]
    data_ambiguity_S = data_ambiguity.loc[data_ambiguity.lowestTaxLevel == 7,]
    #data.to_csv("DATA", sep = '\t', index = False)
    #data_ambiguity.to_csv("data_ambiguity", sep = '\t', index = False)
    #data_no_ambiguity.to_csv("data_no_ambiguity", sep = '\t', index = False)
    if len(data_ambiguity_S) !=0:
        data = data_ambiguity_S.copy()
    else:
        data = pd.DataFrame(data=None, columns=data_ambiguity.columns)
    data_ambiguity_G = data_ambiguity.loc[data_ambiguity.lowestTaxLevel == 6,]
    columnsToReplace=("subOTU", "taxonomy", "taxon1", "taxon2", "taxon3", "taxon4", "taxon5", "taxon6")
    if len(data_ambiguity_G) !=0:
        for column in columnsToReplace:
            data_ambiguity_G["%s"%column].replace(to_replace="_MS", value="_MG", inplace=True, regex=True)
        data = data.append(data_ambiguity_G)
    data_ambiguity_F = data_ambiguity.loc[data_ambiguity.lowestTaxLevel == 5,]
    if len(data_ambiguity_F) !=0:
        for column in columnsToReplace:
            data_ambiguity_F["%s"%column].replace(to_replace="_MS", value="_MF", inplace=True, regex=True)
        data = data.append(data_ambiguity_F)
    data_ambiguity_O = data_ambiguity.loc[data_ambiguity.lowestTaxLevel == 4,]
    if len(data_ambiguity_O) !=0:
        for column in columnsToReplace:
            data_ambiguity_O["%s"%column].replace(to_replace="_MS", value="_MO", inplace=True, regex=True)
        data = data.append(data_ambiguity_O)
    data_ambiguity_C = data_ambiguity.loc[data_ambiguity.lowestTaxLevel == 3,]
    if len(data_ambiguity_C) !=0:
        for column in columnsToReplace:
            data_ambiguity_C["%s"%column].replace(to_replace="_MS", value="_MC", inplace=True, regex=True)
        data = data.append(data_ambiguity_C)
    data_ambiguity_P = data_ambiguity.loc[data_ambiguity.lowestTaxLevel == 2,]
    if len(data_ambiguity_P) !=0:
        for column in columnsToReplace:
            data_ambiguity_P["%s"%column].replace(to_replace="_MS", value="_MP", inplace=True, regex=True)
        data = data.append(data_ambiguity_P)
    data_ambiguity_D = data_ambiguity.loc[data_ambiguity.lowestTaxLevel == 1,]
    if len(data_ambiguity_D) !=0:
        for column in columnsToReplace:
            data_ambiguity_D["%s"%column].replace(to_replace="_MS", value="_MD", inplace=True, regex=True)
        data = data.append(data_ambiguity_D)
    data = data.append(data_no_ambiguity)
    #data.to_csv("DATA2", sep = '\t', index = False)
    #clean
    data.loc[:,'taxon7']=[x.rsplit("_",1)[0] for x in data['subOTU']]
    for i in range (1,7,1):
        data = cleanTaxonName(data,'taxon%i'%i)

    outputfile = 'anchor_table.txt'
    data = data.rename(columns = {'queryid':'anchorID'})
    data = data.rename(columns = {'ContigCount':'anchorCount'})
    if args.primers == True:
        data = data.rename(columns = {'common_to':'anchorCollapsingDueToPrimer'})
    else:
        data['anchorCollapsingDueToPrimer'] = "-"
    data = data[['anchorID','subOTU','SubOTU_2','anchorCount','identity','coverageQ','subjectid','taxonomy','taxon1','taxon2','taxon3','taxon4','taxon5','taxon6','taxon7','lowestTaxLevel','AmbiguousHits','AmbiguousLabels','subjectid_2','Id/Cov_2','AmbiguousHits_2','AmbiguousLabels_2','anchorCollapsingDueToPrimer','sequence']]
    data = data.replace(np.nan, 'NA', regex=True)
    data.to_csv(outputfile, sep = '\t', index = False)
    print ('Output file: ' + str(outputfile))



def part4Work(data,prefixgroup):
    #Now that we have the main data ready, I am going to output the secondary annotation table for all sequences
    print "\n-----------------------------------------"
    print("Keep only the best hits for returns with the same accession ID (i.e. alignments in different part of genome)")
    data = data.sort_values(by = ['queryid', 'subjectid', 'bitscore'], ascending = [True,True,False]).reset_index(drop=True)
    data = data.drop_duplicates(['queryid','subjectid'])
    data = data.sort_values(by = ['queryid', 'subjectid']).reset_index(drop=True)
    print "Done! data has", "{:,}".format(data.shape[0]), "records."
    print "-----------------------------------------\n"

    print "\n-----------------------------------------"
    print "Separate the different annotations"

    for databs in dbList:
        print(databs)
        df_db = data.loc[data.db == databs,:]
        print "\n-----------------------------------------"
        print "Filtering factor: Filtering unknowns when more defined hits are as significant."
        pattern = 'Unknown|unknown'
        df_db.is_copy = False
        df_db['penalty'] = df_db.taxonomy.str.contains(pattern)
        #we'll change True to -1 and False to 0
        df_db['penalty'] = df_db.penalty.astype(int) *(-1)
        df_db['max_point'] = df_db.groupby('queryid')['penalty'].transform(lambda x: x.max()) 
        df_db = df_db[df_db['penalty'] == df_db['max_point']]
        #Removing duplicates bitscore for sequences that only have unknowns as hits (i.e. penalty==-1)
        gp=df_db.groupby('penalty')
        df1 = pd.DataFrame()
        df2 = pd.DataFrame()
        for name, group in gp: 
            if name ==-1:
                group.loc[:,'max_point'] = df_db.groupby('queryid')['bitscore'].transform(lambda x: x.max())
                group = group[group['bitscore'] == group['max_point']]
                group = group.drop_duplicates(['queryid','max_point'])
                df2 = group.copy()
            else:
                df1=group.copy()
        #rebuilding df_db from df1 and df2
        df_db = []
        if df2.empty:
            df_db=df1.copy()
        else:
            df_db = pd.concat([df1,df2])
        df1 = []
        df2 = []
        df_db = df_db.drop(['penalty', 'max_point'] ,1)
        df_db = df_db.reset_index(0, drop=True)
        print "Done! df_db has", "{:,}".format(df_db.shape[0]), "records left"
        print "-----------------------------------------\n"
        

        if numberOfNucleotidesOff != 0:
            print "\n-----------------------------------------"
            print "Keep hits with highest identity + coverage + 1 nucleotide OFF for ambiguous labels"
            #Ncbi report of coverage can exceed 100 when there are deletions, we'll put them back down to 100
            df_db.loc[df_db.coverageQ > 100.0, 'coverageQ'] = 100
            #Check highest identity/coverage reached per query
            df_db['IdentityCov'] = df_db['identity'] + df_db['coverageQ']
            #round down to one decimal
            df_db['IdentityCov'] = df_db['IdentityCov'] // 0.1 / 10
            df_db['MaxIdentityCov'] = df_db.groupby('queryid')['IdentityCov'].transform(lambda x: x.max())
            #allow for one nucleotide off from best bitscore:
            df_db['id1Off'] = (((df_db['identity'] * df_db['alignmentlength'] * 0.01) - int(numberOfNucleotidesOff)) / df_db['alignmentlength'] ) * 100
            df_db['IdentityCov1Off'] = df_db['id1Off'] + df_db['coverageQ']
            df_db['MaxIdentityCov1Off'] = df_db.groupby('queryid')['IdentityCov1Off'].transform(lambda x: x.max())
            df_db.loc[df_db.IdentityCov != df_db.MaxIdentityCov, 'IdentityCov1Off'] = 0
            #round down to one decimal
            df_db['MaxIdentityCov1Off'] = df_db['MaxIdentityCov1Off'] // 0.1 / 10
            #Keep the ones with the highest identity
            df_db = df_db.loc[df_db['IdentityCov'] >= df_db['MaxIdentityCov1Off']]
            #Tag 1OFF hits with a *
            df_db["%sOff"%numberOfNucleotidesOff]="N"
            df_db.loc[df_db.IdentityCov1Off == 0, "%sOff"%numberOfNucleotidesOff] = "Y"
            df_db = df_db.drop(['MaxIdentityCov', 'IdentityCov', 'id1Off', 'IdentityCov1Off', 'MaxIdentityCov1Off'],1)
            print "Done! df_db has", "{:,}".format(df_db.shape[0]), "records left"
            print "-----------------------------------------\n"
        else:
            print "\n-----------------------------------------"
            print "Keep hits with highest identity + coverage for ambiguous labels"
            #Ncbi report of coverage can exceed 100 when there are deletions, we'll put them back down to 100
            df_db.loc[df_db.coverageQ > 100.0, 'coverageQ'] = 100
            #Check highest identity/coverage reached per query
            df_db['IdentityCov'] = df_db['identity'] + df_db['coverageQ']
            #round down to one decimal
            df_db['IdentityCov'] = df_db['IdentityCov'] // 0.1 / 10
            df_db['MaxIdentityCov'] = df_db.groupby('queryid')['IdentityCov'].transform(lambda x: x.max())
            #Keep the ones with the highest identity
            df_db = df_db.loc[df_db['IdentityCov'] >= df_db['MaxIdentityCov']]
            df_db = df_db.drop(['MaxIdentityCov', 'IdentityCov'],1)
            print "Done! df_db has", "{:,}".format(df_db.shape[0]), "records left"
            print "-----------------------------------------\n"




        "\n-----------------------------------------"
        print "Assigning highest taxon reached for each blast hit"
        df_taxon = dict()
        df_taxon[7] = df_db.copy()
        df_taxon[7].loc[:,'lowestTaxLevel'] = "7"
        for i in range(6, 0,-1):
            df_taxon[i] = lowestTaxLevelReached(df_taxon[i+1], i+1, i)
        df_db = []
        #As I have several duplicated hits (df_taxon[1] is in all df_taxon[x]), I have to contactenate them in a certain order and keep the first value for duplication issue
        frames = [df_taxon[1], df_taxon[2], df_taxon[3], df_taxon[4], df_taxon[5], df_taxon[6], df_taxon[7]]

        df_db = pd.concat(frames).reset_index(0, drop=False)
        df_db = df_db.drop_duplicates(subset='index', keep='first')
        #Sort df_db again and remove the column 'index' now useless
        df_db = df_db.sort_values(by = ['queryid','bitscore','identity','taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7']).reset_index(drop=True)
        df_db = df_db.drop(['index'] ,1)
        print "\nDone! df_db has", "{:,}".format(df_db.shape[0]), "records left"
        print "-----------------------------------------\n"




        ################################# LAST COMMON MOST ABUNDANT TAXON LABEL:
        # Ex: 
        #taxon6             taxon7
        #Lactobacillus      L. A
        #Lactobacillus      L. B
        #Lactobacillus      L. C
        #would give:
        #taxon6             taxon7
        #Lactobacillus      [L. A, L. B, L. C]

        print "\n-----------------------------------------"
        print "Create columns for each most abundant taxonomic level for each contig"
        # empty values of taxon value higher than lowestTaxLevel value
        for i in range (1,8,1):
            for j in range (i+1,8,1):
                df_db.copy().loc[df_db.lowestTaxLevel == i, 'taxon%i'%j] = np.nan
        #Select the most abundant taxon per level and per queryid
        for i in range (1,8,1):
            #Number of particular taxon1 per query 
            df_db['t%i_total'%i] = df_db.groupby(['queryid','taxon%i'%i])['taxon%i'%i].transform('count')
            #Highest number of similar taxon1 per query
            df_db['maxt%i'%i] = df_db.groupby('queryid')['t%i_total'%i].transform(lambda x: x.max())
            #Keeping the most abundant taxon1 per query
            df_db["MostAbundantTax%i"%i] = df_db['taxon%i'%i].loc[df_db['t%i_total'%i] == df_db['maxt%i'%i]]
            df_db = df_db.drop(['maxt%i'%i, 't%i_total'%i],1)


        #Create a df_dbframe with all most abundant taxon at given levels (here could be several per queryid)
        for i in range (1,8,1):
            #create a series (could either be taxon separated by commas or just one taxon) for most abundant taxon per level and per query
            tempSeries=df_db.groupby('queryid')['MostAbundantTax%i'%i].apply(lambda x: x.unique())
            if i==1:
                MostAbundantTax=pd.DataFrame({'queryid':tempSeries.index, 'MostAbundantTax%i'%i:tempSeries.values})
                MostAbundantTax['MostAbundantTax%i'%i] = MostAbundantTax['MostAbundantTax%i'%i].astype(str).str.replace("nan", "-")
                MostAbundantTax['MostAbundantTax%i'%i]=listToString(MostAbundantTax['MostAbundantTax%i'%i])
            else:
                temp=pd.DataFrame({'queryid':tempSeries.index, 'MostAbundantTax%i'%i:tempSeries.values})
                MostAbundantTax = pd.merge(MostAbundantTax, temp, on='queryid')
                MostAbundantTax['MostAbundantTax%i'%i]=listToString(MostAbundantTax['MostAbundantTax%i'%i])
            tempSeries=[]
            temp=[]

        #Merge the df_dbframe just created to the main df_db, but before we remove the temporary MostAbundantTax%i columns we created above
        for i in range (1,8,1):
            df_db = df_db.drop(['MostAbundantTax%i'%i],1)
        df_db = pd.merge(df_db, MostAbundantTax, on='queryid')


        print "\n-----------------------------------------"
        print "Exporting ambiguous hits: multiple hits with same score for a same query"
        if numberOfNucleotidesOff != 0:
            ambiguousData=df_db.copy()
            ambiguousData.loc[ambiguousData["%sOff"%numberOfNucleotidesOff] == "Y", 'taxon7'] = ambiguousData.taxon7 + "*"
        #concatenate all subjectid for 'queryid' with the same 'taxonomy'
        secondary_annotation = df_db.groupby(["queryid","taxonomy"]).agg({"subjectid":RowCollapse})
        secondary_annotation['subjectid'].replace(to_replace=";", value=" ", inplace=True, regex=True) #i'll remove the ; so we can use the id numbers list directly on NCBI entrez
        secondary_annotation.columns=secondary_annotation.columns+"_2"
        secondary_annotation=secondary_annotation.reset_index(drop=False)
        #create a mabiguous df before removing df_db duplicates
        if numberOfNucleotidesOff != 0:
            ambiguous = ambiguousData[['queryid', 'subjectid', 'identity', 'coverageQ', 'taxonomy', "taxon7", "ContigCount"]]
        else:
            ambiguous = df_db[['queryid', 'subjectid', 'identity', 'coverageQ', 'taxonomy', "ContigCount"]]
        ambiguous = ambiguous.drop_duplicates(['queryid','taxonomy'])
        ambiguous = ambiguous.merge(secondary_annotation, on=['queryid','taxonomy'], how='left')
        ambiguous = ambiguous[['queryid', 'subjectid_2', 'identity', 'coverageQ','taxonomy', "ContigCount"]]
        print "Done! df_db has", "{:,}".format(df_db.shape[0]), "records left"
        print "-----------------------------------------\n"

        if numberOfNucleotidesOff != 0:
            print "\n-----------------------------------------"
            print "Remove all hits that do not reach the lowest taxonomic level : xOFF=N"
            df1=df_db[df_db["%sOff"%numberOfNucleotidesOff] == "N"]
            #Check highest taxonomic level reached per query
            df1['HitMaxTaxLevel'] = df1.groupby('queryid')['lowestTaxLevel'].transform(lambda x: x.max())
            #compare with lowestTaxLevel (which is the lowest taxonomic level reached by one hit) and keep the lowest per query
            df1 = df1.loc[df1['lowestTaxLevel'] == df1['HitMaxTaxLevel']]
            df1 = df1.drop(['HitMaxTaxLevel'],1)
            print "Done! df_db xOFF=N has", "{:,}".format(df1.shape[0]), "records left"
            print "-----------------------------------------\n"
    
    
            print "\n-----------------------------------------"
            print "Remove all hits that do not reach the lowest taxonomic level : xOFF=Y"
            if len(df_db[df_db["%sOff"%numberOfNucleotidesOff] == "Y"]) !=0:
                df2=df_db[df_db["%sOff"%numberOfNucleotidesOff] == "Y"]
                #Check highest taxonomic level reached per query
                df2['HitMaxTaxLevel'] = df2.groupby('queryid')['lowestTaxLevel'].transform(lambda x: x.max())
                #compare with lowestTaxLevel (which is the lowest taxonomic level reached by one hit) and keep the lowest per query
                df2 = df2.loc[df2['lowestTaxLevel'] == df2['HitMaxTaxLevel']]
                df2 = df2.drop(['HitMaxTaxLevel'],1)
                print "Done! df_db xOFF=Y has", "{:,}".format(df2.shape[0]), "records left"
            else:
                print "Not hit has xOFF=Y"
            print "-----------------------------------------\n"
    
            print "\n-----------------------------------------"
            print "Reconstruct main df_db"
            if len(df_db[df_db["%sOff"%numberOfNucleotidesOff] == "Y"]) !=0:
                #concatenate df_dbframes
                frames = [df1, df2]
                df_db = pd.concat(frames)
            else:
                df_db = df1.copy()
            df1=[]
            df2=[]
            #Remove xOFF that have lower taxonomic level than main hit
            df_db['HitMaxTaxLevel'] = df_db.groupby('queryid')['lowestTaxLevel'].transform(lambda x: x.max())
            df_db["Keep"]="Y"
            df_db.loc[(df_db["%sOff"%numberOfNucleotidesOff] == "Y") & (df_db.lowestTaxLevel < df_db.HitMaxTaxLevel), 'Keep'] = "NO"
            df_db = df_db.loc[df_db.Keep == "Y"]
            df_db = df_db.drop(['HitMaxTaxLevel','Keep'],1)
            #Reconstruct the taxon columns
            taxa1 = df_db['taxonomy'].apply(lambda x: pd.Series(x.split(';',)))
            taxa1.columns = ['taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7']
            df_db = df_db.drop(['taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7'],1)
            df_db = df_db.join(taxa1, how='inner')
    
            print "Done! df_db has", "{:,}".format(df_db.shape[0]), "records left"
            print "-----------------------------------------\n"

        else:
            print "\n-----------------------------------------"
            print "Remove all hits that do not reach the lowest taxonomic level"
            #Check highest taxonomic level reached per query
            df_db['HitMaxTaxLevel'] = df_db.groupby('queryid')['lowestTaxLevel'].transform(lambda x: x.max())
            #compare with lowestTaxLevel (which is the lowest taxonomic level reached by one hit) and keep the lowest per query
            df_db = df_db.loc[df_db['lowestTaxLevel'] == df_db['HitMaxTaxLevel']]
            df_db = df_db.drop(['HitMaxTaxLevel'],1)
            print "Done! df_db has", "{:,}".format(df_db.shape[0]), "records left"
            print "-----------------------------------------\n"
    
            print "\n-----------------------------------------"
            print "Reconstruct main df_db"
            taxa1 = df_db['taxonomy'].apply(lambda x: pd.Series(x.split(';',)))
            taxa1.columns = ['taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7']
            df_db = df_db.drop(['taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7'],1)
            df_db = df_db.join(taxa1, how='inner')
            print "Done! df_db has", "{:,}".format(df_db.shape[0]), "records left"
            print "-----------------------------------------\n"

        #subselecting output
        if numberOfNucleotidesOff != 0:
            df_db = df_db[['queryid', 'querylength', 'subjectid', 'subjectlength', "identity", "evalue", "coverageQ", "coverageS", 'subjectdescription', "bitscore", 'taxonomy', 'taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7', "MostAbundantTax1", "MostAbundantTax2", "MostAbundantTax3", "MostAbundantTax4", "MostAbundantTax5", "MostAbundantTax6", "MostAbundantTax7", "lowestTaxLevel", "ContigCount", "%sOff"%numberOfNucleotidesOff]]
        else:
            df_db = df_db[['queryid', 'querylength', 'subjectid', 'subjectlength', "identity", "evalue", "coverageQ", "coverageS", 'subjectdescription', "bitscore", 'taxonomy', 'taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7', "MostAbundantTax1", "MostAbundantTax2", "MostAbundantTax3", "MostAbundantTax4", "MostAbundantTax5", "MostAbundantTax6", "MostAbundantTax7", "lowestTaxLevel", "ContigCount"]]
        outputfile = prefixgroup + '_' + databs + '_Final.txt'
        df_db.to_csv(outputfile, sep = '\t', index = False)
        finalnumberofsequences = len(df_db)
        
        print "\n     OUTPUT INFORMATION:     "
        print "-------------------------------"
        print ('Parsed file : ' + str(outputfile))
        print "Parsed file has : ", "{:,}".format(df_db.shape[0]), " records\n-------------------------------\n"

def ambiguousLabelCreation(data):
    print "Evaluate the highest common taxon level for subOtu anchor sequences"
    df = data.groupby('queryid')['taxon7'].apply(lambda x: len(x.unique())).to_frame("DifferentT7")
    df_t6 = data.groupby('queryid')['taxon6'].apply(lambda x: len(x.unique())).to_frame("DifferentT6")
    df = df.join(df_t6, how='outer')
    df_t5 = data.groupby('queryid')['taxon5'].apply(lambda x: len(x.unique())).to_frame("DifferentT5")
    df = df.join(df_t5, how='outer')
    df_t4 = data.groupby('queryid')['taxon4'].apply(lambda x: len(x.unique())).to_frame("DifferentT4")
    df = df.join(df_t4, how='outer')
    df_t3 = data.groupby('queryid')['taxon3'].apply(lambda x: len(x.unique())).to_frame("DifferentT3")
    df = df.join(df_t3, how='outer')
    df_t2 = data.groupby('queryid')['taxon2'].apply(lambda x: len(x.unique())).to_frame("DifferentT2")
    df = df.join(df_t2, how='outer')
    df_t1 = data.groupby('queryid')['taxon1'].apply(lambda x: len(x.unique())).to_frame("DifferentT1")
    df = df.join(df_t1, how='outer')
    df["LastCommonTaxon"] = 0
    df.loc[(df['DifferentT1'] == 1), ['LastCommonTaxon']] = 1
    df.loc[(df['DifferentT2'] == 1), ['LastCommonTaxon']] = 2
    df.loc[(df['DifferentT3'] == 1), ['LastCommonTaxon']] = 3
    df.loc[(df['DifferentT4'] == 1), ['LastCommonTaxon']] = 4
    df.loc[(df['DifferentT5'] == 1), ['LastCommonTaxon']] = 5
    df.loc[(df['DifferentT6'] == 1), ['LastCommonTaxon']] = 6
    df.loc[(df['DifferentT7'] == 1), ['LastCommonTaxon']] = 7
    df = df.drop(['DifferentT1','DifferentT2','DifferentT3','DifferentT4','DifferentT5','DifferentT6','DifferentT7'],1)
    df = df.reset_index(0, drop=False)
    data = data.merge(df, on=['queryid'], how='left')
    return data

def createAnnotationName(data):
    #Create Ambiguous labels for all ambiguous hits
    pattern = ';'
    data.loc[(data['LastCommonTaxon'] == 0) & (data['AmbiguousHits'] > 1), ['taxon1','taxon2','taxon3','taxon4','taxon5','taxon6','taxon7']] = "SpeciesWithAmbiguousDomains"
    data.loc[(data['LastCommonTaxon'] == 1) & (data['AmbiguousHits'] > 1), ['taxon1','taxon2','taxon3','taxon4','taxon5','taxon6','taxon7']] = data.loc[(data['LastCommonTaxon'] == 1) & (data['AmbiguousHits'] > 1), 'taxon1'].map(str) + "_MS"
    data.loc[(data['LastCommonTaxon'] == 2) & (data['AmbiguousHits'] > 1), ['taxon2','taxon3','taxon4','taxon5','taxon6','taxon7']] = data.loc[(data['LastCommonTaxon'] == 2) & (data['AmbiguousHits'] > 1), 'taxon2'].map(str) + "_MS"
    data.loc[(data['LastCommonTaxon'] == 3) & (data['AmbiguousHits'] > 1), ['taxon3','taxon4','taxon5','taxon6','taxon7']] = data.loc[(data['LastCommonTaxon'] == 3) & (data['AmbiguousHits'] > 1), 'taxon3'].map(str) + "_MS"
    data.loc[(data['LastCommonTaxon'] == 4) & (data['AmbiguousHits'] > 1), ['taxon4','taxon5','taxon6','taxon7']] = data.loc[(data['LastCommonTaxon'] == 4) & (data['AmbiguousHits'] > 1), 'taxon4'].map(str) + "_MS"
    data.loc[(data['LastCommonTaxon'] == 5) & (data['AmbiguousHits'] > 1), ['taxon5','taxon6','taxon7']] = data.loc[(data['LastCommonTaxon'] == 5) & (data['AmbiguousHits'] > 1), 'taxon5'].map(str) + "_MS"
    data.loc[(data['LastCommonTaxon'] == 6) & (data['AmbiguousHits'] > 1), ['taxon6','taxon7']] = data.loc[(data['LastCommonTaxon'] == 6) & (data['AmbiguousHits'] > 1), 'taxon6'].map(str) + "_MS"
    #data.loc[(data['LastCommonTaxon'] == 7) & (data['AmbiguousHits'] > 1), 'taxon7'] = data.loc[(data['LastCommonTaxon'] == 7) & (data['AmbiguousHits'] > 1), 'taxon7'].map(str) + "_MS"
    data.loc[data['AmbiguousHits'] > 1, 'taxonomy'] = data['taxon1'] + ";" + data['taxon2'] + ";" + data['taxon3'] + ";" + data['taxon4'] + ";" + data['taxon5'] + ";" + data['taxon6'] + ";" + data['taxon7']
    
    
    #Some taxonomy have repeated upper levels (like Clostridiales in Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiales;Eubacterium;Eubacterium). Then with the operation just above, 
    # and the one just below, I come to have a problem where ambiguous species with a same name can be of different levels. For example:
    #Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiales;Eubacterium;Eubacterium    will produce    Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiales;Clostridiales_MS;Clostridiales_MS
    #while 
    #Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;Lachnospiraceae;Lachnospiraceae    will produce    Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiales_MS;Clostridiales_MS;Clostridiales_MS
    #They'll both have the same _MS name although their laast taxonomic level will be different.
    #Solution: I will lower the taxonomic level of the first so it will be identical to the last one
    # 1. Crete fake taxons (the same as the other ones, but without _MS)
    for i in range (1,8,1):
        data["Faketaxon%i"%i] = [x.rsplit("_MS",1)[0] for x in data['taxon%i'%i]]
    #2 . create a fake lowestTaxLevel column based on the faketaxons
    print "Parsing the df_db: assigning highest taxon reached for each blast hit"
    df_db_t7 = data.copy()
    #Create a column with the highest taxon level reached for a hit (fixing it at 7 initially) 
    df_db_t7['FakeCommonLastTax'] = 7
    #Parsing df_db for taxa that are identical between 2 levels
    index = (df_db_t7["Faketaxon7"] == df_db_t7["Faketaxon6"])
    df_db_t6 = df_db_t7.ix[index]
    index = (df_db_t6["Faketaxon6"] == df_db_t6["Faketaxon5"])
    df_db_t5  = df_db_t6.ix[index]
    index = (df_db_t5["Faketaxon5"] == df_db_t5["Faketaxon4"])
    df_db_t4  = df_db_t5.ix[index]
    index = (df_db_t4["Faketaxon4"] == df_db_t4["Faketaxon3"])
    df_db_t3  = df_db_t4.ix[index]
    index = (df_db_t3["Faketaxon3"] == df_db_t3["Faketaxon2"])
    df_db_t2  = df_db_t3.ix[index]
    index = (df_db_t2["Faketaxon2"] == df_db_t2["Faketaxon1"])
    df_db_t1  = df_db_t2.ix[index]
    #Assigning highest Faketaxon level reached for a hit, for each chunk
    df_db_t6.loc[:,'FakeCommonLastTax'] = 6
    df_db_t5.loc[:,'FakeCommonLastTax'] = 5
    df_db_t4.loc[:,'FakeCommonLastTax'] = 4
    df_db_t3.loc[:,'FakeCommonLastTax'] = 3
    df_db_t2.loc[:,'FakeCommonLastTax'] = 2
    df_db_t1.loc[:,'FakeCommonLastTax'] = 1
    df_db = []
    #As I have several duplicated hits (df_db_t1 is in all df_db_tx), I have to contactenate them in a certain order and keep the first value when duplication
    frames = [df_db_t1, df_db_t2, df_db_t3,df_db_t4, df_db_t5, df_db_t6, df_db_t7]
    df_db = pd.concat(frames).reset_index(0, drop=False)
    df_db = df_db.drop_duplicates(subset='index', keep='first')
    #Sort df_db again and remove the column 'index' now useless
    df_db = df_db.reset_index(drop=True)
    df_db = df_db.drop(['index'] ,1)
    data = df_db.copy()
    #3. I'll replace the taxon 1 to 7 according to the fake last tax
    pattern = ';'
    data.loc[(data['FakeCommonLastTax'] == 1) & (data['AmbiguousHits'] > 1), ['taxon1','taxon2','taxon3','taxon4','taxon5','taxon6','taxon7']] = data.loc[(data['FakeCommonLastTax'] == 1) & (data['AmbiguousHits'] > 1), 'Faketaxon1'].map(str) + "_MS"
    data.loc[(data['FakeCommonLastTax'] == 2) & (data['AmbiguousHits'] > 1), ['taxon3','taxon4','taxon5','taxon6','taxon7']] = data.loc[(data['FakeCommonLastTax'] == 2) & (data['AmbiguousHits'] > 1), 'Faketaxon2'].map(str) + "_MS"
    data.loc[(data['FakeCommonLastTax'] == 3) & (data['AmbiguousHits'] > 1), ['taxon4','taxon5','taxon6','taxon7']] = data.loc[(data['FakeCommonLastTax'] == 3) & (data['AmbiguousHits'] > 1), 'Faketaxon3'].map(str) + "_MS"
    data.loc[(data['FakeCommonLastTax'] == 4) & (data['AmbiguousHits'] > 1), ['taxon5','taxon6','taxon7']] = data.loc[(data['FakeCommonLastTax'] == 4) & (data['AmbiguousHits'] > 1), 'Faketaxon4'].map(str) + "_MS"
    data.loc[(data['FakeCommonLastTax'] == 5) & (data['AmbiguousHits'] > 1), ['taxon6','taxon7']] = data.loc[(data['FakeCommonLastTax'] == 5) & (data['AmbiguousHits'] > 1), 'Faketaxon5'].map(str) + "_MS"
    data.loc[(data['FakeCommonLastTax'] == 6) & (data['AmbiguousHits'] > 1), ['taxon7']] = data.loc[(data['FakeCommonLastTax'] == 6) & (data['AmbiguousHits'] > 1), 'Faketaxon6'].map(str) + "_MS"
    data.loc[data['AmbiguousHits'] > 1, 'taxonomy'] = data['taxon1'] + ";" + data['taxon2'] + ";" + data['taxon3'] + ";" + data['taxon4'] + ";" + data['taxon5'] + ";" + data['taxon6'] + ";" + data['taxon7']
    data = data.drop(['Faketaxon1','Faketaxon2','Faketaxon3','Faketaxon4','Faketaxon5','Faketaxon6','Faketaxon7'] ,1)
    return data

def mergeMultipleDatabaseAndParse(df):
    #The goal is to recreate columns based on concatenated multi-database entries
    #Columns to rebuild: AmbiguousLabels, LastCommonTaxon, subjectid, AmbiguousHits, taxonomy columns
    #Concatenate all the ambiguous labels of a given query
    df3 = df.groupby('queryid')['AmbiguousLabels'].apply(';'.join).to_frame("newAmbi")
    #Now we will remove duplicated values within the concatenated AmbiguousLabels
    s = df3['newAmbi'].str.split(';').apply(Series, 1).stack() #this creates an object (s) with a new line every separator (here ;)
    #s has 2 indexes. the first is original data's index, and the second is s's. We'll remove the latter
    s.index = s.index.droplevel(-1) # to line up with data's index
    s.name = 'AmbiguousLabels_1'# needs a name to join
    df3 = df3.join(s) #replace with new
    df3 = df3.reset_index(drop=False)
    df3 = df3.drop_duplicates(["queryid","AmbiguousLabels_1"])
    #get back to ; separated values in AmbiguousLabels
    df3 = df3.sort_values(by = ["queryid","AmbiguousLabels_1"])
    df3 = df3.groupby('queryid')['AmbiguousLabels_1'].apply(';'.join).to_frame("AmbiguousLabels")
    df3 = df3.reset_index(drop=False)
    #df3 has the good ambiguousLabels column, we'll need to remove it from df2and join df2 and df3
    df = df.drop(['AmbiguousLabels'] ,1)
    df = df.merge(df3, on=['queryid'], how='outer')
    #Now that Ambiguous labels have been re-created, we need to do the same for LastCommonTaxon column. It has to be adapted to the new AmbiguousLabels. I'll choose the lowest from a new LastCommonTaxon calculated from the new AmbiguousLabels and the previous LastCommonTaxon
    #Call LastCommonTaxonWithinDB the last common taxon within each database
    df = df.rename(columns = {'LastCommonTaxon':'LastCommonTaxonWithinDB'})
    #recreate last common taxon column between databases
    df = ambiguousLabelCreation(df)
    df = df.rename(columns = {'LastCommonTaxon':'LastCommonTaxonBetweenDB'})
    #Choose the lowest values between the 2 versions of LastCommonTaxon
    df["LastCommonTaxon"] = df[['LastCommonTaxonWithinDB','LastCommonTaxonBetweenDB']].min(axis=1)
    #Now use the minimum last common taxon for all hit of a same query
    df['LastCommonTaxon'] = df.groupby('queryid')['LastCommonTaxon'].transform(lambda x: x.min())
    df = df.drop(['LastCommonTaxonBetweenDB','LastCommonTaxonWithinDB'] ,1)
    #collapse subjectids values for a same query
    df3 = df.groupby("queryid").agg({"subjectid":RowCollapse})
    df3 = df3.reset_index(drop=False)
    df = df.drop(['subjectid'] ,1)
    #now join df to collapsed values in df3
    df = df.merge(df3, on=['queryid'], how='inner')
    #remove queryid duplicates    
    df = df.drop_duplicates(["queryid"])
    #Revaluate AmbiguousHits (taking into account between databases)
    #this creates a serie (s) with a new line every separator (here ;)
    s = df['AmbiguousLabels'].str.split(';').apply(Series, 1).stack() #
    #count how many lines per query id
    #s has 2 indexes. the first is original data's index, and the second is s's. We'll remove the latter
    s.index = s.index.droplevel(-1) # to line up with data's index
    #create another serie to count how many index have the same value: that will be the new number of ambiguous labels
    s2 = s.groupby(s.index).size()
    s2.name = 'AmbiguousHits'
    df = df.drop(['AmbiguousHits'] ,1)
    df = df.merge(s2.to_frame(), left_index=True, right_index=True)
    df.loc[df.AmbiguousHits == 1, 'AmbiguousHits'] = 0
    #re-evaluate taxon1 to taxon7 
    df = createAnnotationName(df)
    return df





parser = argparse.ArgumentParser(description="'inputfile: myfile_part3_parsed.txt \noutput file : anchors_table.txt. What programs are needed on the path:ktImportText(KRONA)")
parser.add_argument('-i', '--inputfile', help='[REQUIRED] Output file part3_bitscoreParser.py. It should have 29 columns: queryid|subjectid|identity|alignmentlength|mismatches|gapopens|qstart|qend|evalue|bitscore|querylength|querydescription|coverageQ|subjectlength|subjectdescription|coverageS|accession|accession_Full|taxonomy|taxon1|taxon2|taxon3|taxon4|taxon5|taxon6|taxon7', dest='inputfile', action='store', required=True)
parser.add_argument('-m', '--multiplier', help='[OPTIONAL] Multiplier file: counts of unique sequences. Such a file is the output of mothur (All.trim.contigs.good.count_table). It should have 2 colmns: "Representative_Sequence", "total"', dest='inputMultiplier', default=False, action='store', required=False)
parser.add_argument('-g', '--group_accession_id', help='[OPTIONAL] When all sequences have been parsed and that all remains is ambiguous calls, you can change the final accession ID call for a given sequence by choosing the one that is most common across all ambiguous sequences in the data.', dest='group_accession_id', default=False, action='store_true', required=False)
parser.add_argument('-o', '--output_option', help='[OPTIONAL] THe outut can be clustered by accession ID (SubOTU option, default), OTU labels (OTU) or no collapsing at all (denovo).', dest='output_option', default="SubOTU", action='store', required=False)
parser.add_argument('-n', '--numberOfNucleotidesOff', help='[OPTIONAL] This option creates a column with ambiguous hits at x nucleotide off from the best hit in each anchor sequence', dest='numberOfNucleotidesOff', default=0, action='store', required=False)
parser.add_argument('-p', '--primers', help='[OPTIONAL] Collapse anchors based on primer ambiguity nucleotides if exist. It will use the anchor fasta file produced in part5 (anchors_fasta.tsv) and the primer list used in part2 (metadata/primers.txt)', dest='primers', default=False, action='store_true', required=False)



args = parser.parse_args()
inputfile = args.inputfile
output_option = args.output_option
numberOfNucleotidesOff = args.numberOfNucleotidesOff

try:
    fasta = pd.read_csv("anchors_fasta.tsv", sep='\t', names = ["queryid", "sequence"])
except:
    print "\n---\nYou need a copy or link of anchors_fasta.tsv (check part5 folder) in this folder"
    sys.exit(1)



if args.primers == True:
    try:
        primers = pd.read_csv("primers.txt", sep='\t')
    except:
        print "\n---\nYou need a copy or link of primers.txt (check Metadata folder) in this folder"
        sys.exit(1)



print 'Preprocessing inputfile'
command = """/bin/bash -c "awk 'FNR>1' %s | sort -t $'\t' -k1,1 > __inputfile" """ %inputfile
command_run = subprocess.call(command, shell=True)
if command_run != 0:
    exit(1)
print 'done.'

print '\n'
print "\nINPUT INFORMATION"
print "-------------------------------"
print ('Inputfile in tabular format : ' + str(inputfile))
if args.inputMultiplier != False:
    multiplier = args.inputMultiplier
    print ('Count file : ' + str(multiplier))
if numberOfNucleotidesOff != 0:
    print ('Use ' + str(numberOfNucleotidesOff) + ' nucleotide(s) off as good hits')
print "-------------------------------\n"


print "\n-----------------------------------------"
print "Splitting inputfile\n"
itercount=0
f = open("__inputfile", "r")
for i, batch in enumerate(batch_iterator(f, 100000)) :
    filename = "./group_%i.txt" % (i+1)
    handle = open(filename, "w")
    handle.write("".join(batch))
    itercount =+ i
    handle.close()
    print "Wrote 100,000 records to %s" % (filename)
f.close()
print "Done."
print "Total number of chunk files : %i" % (itercount+1)


if itercount==0:
    print "one chunk, all good."
else:
    #Creating a bash code that will control that there are no queryid intersection between 2 chunks.
    #1: Take all the querid = to the last queryid
    #2: add it to the beginning of the next file
    #3: I'll eventually remove the duplicate (not in this script) queryid by keeping the last one. The reason I'm not doing it is time. If the file is really big, it will take unecessary time to be processed.
    controler_code = """#!/bin/bash
set -e
#setting all the group files in order and in a list
ls -1 group_*.txt | cut -d"." -f1 | cut -d"_" -f2 | grep [0-9] | sort -n > __groupList
#how many groups?
numberOfGroups=$(wc -l __groupList | cut -d" " -f1)
#preparing the script
regulator=1
while read i
do
    if [ ${regulator} != 1 ]; then 
        cat __intersect group_${i}.txt > __temp
        mv __temp group_${i}.txt
    fi
    #extracting the queryid of the last record in the group
    lastqueryid=$(tail -n1 group_${i}.txt | cut -f1)
    #extracting all this queryid records and putting them into a file (__intersect) so we can add it to the next group 
    grep "${lastqueryid}" group_${i}.txt > __intersect
    if [ ${regulator} != ${numberOfGroups} ]; then
        #removing this queryid records to this group
        grep -v "${lastqueryid}" group_${i}.txt > __temp
        mv __temp group_${i}.txt
    fi
    regulator=$((regulator+1))
done<__groupList
rm -f __intersect __groupList
    """
    myFile = open("__batch_controler.sh", "w")
    myFile.write(controler_code)
    myFile.close()
    print "Crontroling that a given queryid is not shared between 2 chunks"
    command = 'bash __batch_controler.sh'
    command_run = subprocess.call(command, shell=True)
    if command_run != 0:
        exit(1)

os.system('rm -f __groupList __inputfile __batch_controler.sh')
print "done!"
print "-----------------------------------------\n"




print "\n-----------------------------------------"
print "Loading and parsing each chunk"
#create a list where all different databses will be listed
alldatabaseList=[]
os.system('rm -f *_Final.txt')
for i in range(1,itercount+2):
    #print "here%i" % i
    inputfile = 'group_%i.txt' % (i)
    prefixgroup='group_%i' % (i)

    print "loading:", inputfile
    data = pd.read_csv(inputfile, sep='\t', header=None)


    try:
        len(data.columns) == 29
        data.columns = ["queryid", "subjectid", "identity", "alignmentlength", "mismatches", "gapopens", "qstart", "qend", "evalue", "bitscore", "querylength", "querydescription", "coverageQ", "subjectlength", "coverageS", "subjectdescription", "accession", "accession_Full", "taxonomy", "taxon1", "taxon2", "taxon3", "taxon4", "taxon5", "taxon6", "taxon7"]
        print "Done! data has", "{:,}".format(data.shape[0]), "records."
        print "-----------------------------------------\n"

    except:
        print "\n-----\nError detected!"
        print "Problem with :", inputfile
        print "Make sure that you have 29 columns : queryid|subjectid|identity|alignmentlength|mismatches|gapopens|qstart|qend|evalue|bitscore|querylength|querydescription|coverageQ|subjectlength|coverageS|subjectdescription|accession|accession_Full|taxonomy|taxon1|taxon2|taxon3|taxon4|taxon5|taxon6|taxon7"
        sys.exit(1)

    #create a column for the database used
    data.loc[:,'db']=[x.split("|",1)[0] for x in data['subjectid']]
    pattern = 'TrueUnknown_'
    data.loc[data.db.str.contains(pattern), 'db'] = "NA"
    dbList = list(data.db.unique())
    alldatabaseList = alldatabaseList + dbList



    if args.inputMultiplier == False:
        data['ContigCount']=1
    else:
        multiplier = args.inputMultiplier
        multiFile=pd.read_csv(multiplier, sep='\t')
        multiFile=multiFile[["Representative_Sequence", "total"]]
        multiFile = multiFile.rename(columns = {'Representative_Sequence':'queryid'})
        multiFile = multiFile.rename(columns = {'total':'ContigCount'})
        data_check = len(data)
        df_merge = pd.merge(data, multiFile, how='inner', on='queryid')
        data = df_merge.copy()
        if len(data) != data_check:
            print "Not all sequences in your input file were found in the multiplier file"
            missing = pd.merge(data, multiFile, how='left', on='queryid', left_on=None, right_on=None,
                 left_index=False, right_index=False, sort=True,
                 suffixes=('_x', '_y'), copy=True, indicator=True)
            missing = missing.loc[missing['_merge'] != 'both']
            missing.to_csv("_MISSING_DATA.txt", sep = '\t', index = False)
            print "Check _MISSING_DATA.txt in the output folder"
            sys.exit("Error")

    
    part4Work(data,prefixgroup)



print "\n-----------------------------------------"
print "All Chunks have been processed. Gathering them into one file"
#first clean the alldatabaseList into unique values
alldatabaseList = unique(alldatabaseList)
mainData = []
counter = 0
for databs in alldatabaseList:
    counter = counter + 1
    print(databs)
    #final main output
    os.system('cat group_*_%s_Final.txt > %s_almostThere_final.txt' %(databs,databs))
    os.system('sed -i \'0,/queryid/! {/queryid/d}\' %s_almostThere_final.txt' %databs)

    last_inputfile = databs + "_almostThere_final.txt"
    
    #Main data
    data = pd.read_csv(last_inputfile, sep='\t')
    os.system('rm -f %s_almostThere_final.txt' %databs)
    print "Done! data has", "{:,}".format(data.shape[0]), "records"
    print "-----------------------------------------\n"
    

    print "\n-----------------------------------------"
    print "Keep hits with highest identity + coverage values"
    #Check highest identity/coverage reached per query
    data['IdentityCov'] = data['identity'] + data['coverageQ']
    data['MaxIdentityCov'] = data.groupby('queryid')['IdentityCov'].transform(lambda x: x.max())
    #Keep the ones with the highest identity
    data = data.loc[data['IdentityCov'] == data['MaxIdentityCov']]
    data = data.drop(['MaxIdentityCov', 'IdentityCov'],1)
    print "Done! data has", "{:,}".format(data.shape[0]), "records left"
    print "-----------------------------------------\n"


    print "\n-----------------------------------------"
    print "Evaluate the highest common taxon level for subOtu anchor sequences"
    df = data.groupby('queryid')['taxon7'].apply(lambda x: len(x.unique())).to_frame("DifferentT7")
    for i in range(6,0,-1):
        df_temp = data.groupby('queryid')['taxon%i'%i].apply(lambda x: len(x.unique())).to_frame("DifferentT%i"%i)
        df = df.join(df_temp, how='outer')
        df_temp = []
    df["LastCommonTaxon"] = 0
    for i in range(1,8,1):
        df.loc[(df['DifferentT%i'%i] == 1), ['LastCommonTaxon']] = i
    df = df.drop(['DifferentT1','DifferentT2','DifferentT3','DifferentT4','DifferentT5','DifferentT6','DifferentT7'],1)
    df = df.reset_index(0, drop=False)
    data = data.merge(df, on=['queryid'], how='left')
    print "\nDone! df_db has", "{:,}".format(data.shape[0]), "records left"
    print "-----------------------------------------\n"
  
    
        
    print "\n-----------------------------------------"
    print "Create a column listing all the ambiguous hits"
    df = data.drop_duplicates(['queryid','taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7'])
    df.is_copy = False
    #count number of ambiguous hits
    df['AmbiguousHits'] = df.groupby(['queryid'])['queryid'].transform('count')
    df=df[['queryid', 'AmbiguousHits']]
    #differentiate the 1OFF if any
    data["tax7Ambi"] = data["taxon7"]
    if numberOfNucleotidesOff != 0:
        data.loc[data["%sOff"%numberOfNucleotidesOff] == "Y", 'tax7Ambi'] = data.taxon7 + "*"
    #Merge df and data (i.e add column with the number of ambiguous hits)
    data = data.merge(df, on=['queryid'], how='left')
    #Change 1 to 0 in the ambiguous column 
    data.loc[data.AmbiguousHits == 1, 'AmbiguousHits'] = 0
    df=[]
    #Prepare for AmbiguousLabels column
    df = data.drop_duplicates(['queryid','taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7'])
    ambiguousAnnotation = df.groupby("queryid").agg({'tax7Ambi':RowCollapse})
    ambiguousAnnotation.columns=ambiguousAnnotation.columns+"_2"
    ambiguousAnnotation = ambiguousAnnotation.rename(columns = {'tax7Ambi_2':'AmbiguousLabels'})
    ambiguousAnnotation=ambiguousAnnotation.reset_index(drop=False)
    data = data.drop(['tax7Ambi'],1)
    print "Done! data has", "{:,}".format(data.shape[0]), "records"
    print "-----------------------------------------\n"
    
    if args.group_accession_id == True:
        print "\n-----------------------------------------"
        print "Now try to find a common home (i.e. subjectid) for all ambiguous results that we cannot separate score-wise"
        #count number of occurence of a given ambiguous subjectID accross all samples
        data['num_totals'] = data.groupby(['subjectid'])['subjectid'].transform('count')
        #Select the most abundant ambiguous subjectID accross all samples if it exists
        data['MaxCommonSubjectID'] = data.groupby('queryid')['num_totals'].transform(lambda x: x.max())
        #Keep the ones with the highest identity
        data = data.loc[data['num_totals'] == data['MaxCommonSubjectID']]
        data = data.drop(['num_totals', 'MaxCommonSubjectID'],1)
        print "Done! data has", "{:,}".format(data.shape[0]), "records"
        print "-----------------------------------------\n"
    
    
    print "\n-----------------------------------------"
    print "Remove remaining ambiguity alphabetically"
    data = data.sort_values(by = ['queryid', 'taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7']).reset_index(drop=True)
    data = data.drop_duplicates(['queryid', 'taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7'])
    data = data.sort_values(by = ['queryid', 'subjectid']).reset_index(drop=True)
    data = data.drop_duplicates(['queryid', 'subjectid'])
    #create a new accession id based on best ambiguous hits (just to be able to differentiate hits when grouping on same accession ID)
    #limit the number of perfect ambiguous hits to 10 max
    data['subindex'] = data.groupby(['queryid']).cumcount()
    data = data.loc[data['subindex'] < 11 ]
    data = data.drop(['subindex'],1)
    ambiguousSubj = data.groupby("queryid").agg({'subjectid':RowCollapse})
    ambiguousSubj = ambiguousSubj.reset_index(drop=False)
    data = data.drop(['subjectid'],1)
    data = data.merge(ambiguousSubj, on=['queryid'], how='outer')
    data = data.drop_duplicates(['queryid'])
    data = data.sort_values(by = ['queryid']).reset_index(drop=True)
    #create amabiguous df before removing data duplicates
    data = data.merge(ambiguousAnnotation, on=['queryid'], how='outer')
    #Join that dataframe to the main data
    if counter == 1:
        mainData = data.copy()
    else:
        mainData = data.append(mainData)

    print "Done! data has", "{:,}".format(data.shape[0]), "records"
    print "-----------------------------------------\n"



#Between databases annotation selection rules:
#   1. Lowest annotation taxonomic level is selected
#   2. Best identity and coverage hits are selected
#   3. On equal identity and coverage stats, 16SMicrobial is selected as primary annotation (df1) and others are selected as secondary annotation (df2)
#   4. If 16SMicrobial is not the best identity and coverage then the (other) best identity and coverage is selected as primary (df1)

#create a column for the database used
mainData.loc[:,'db']=[x.split("|",1)[0] for x in mainData['subjectid']]
pattern = 'TrueUnknown_'
mainData.loc[mainData.db.str.contains(pattern), 'db'] = "NA"
#Between databases annotation selection rules:
#   1. Highest annotation taxonomic level is selected
#   2. Best identity and coverage hits are selected
#   3. On equal identity and coverage stats, 16SMicrobial is selected as primary annotation (df1) and others are selected as secondary annotation (df2)


#1
mainData['HitMaxTaxLevel'] = mainData.groupby('queryid')['lowestTaxLevel'].transform(lambda x: x.max())
mainData = mainData.loc[mainData['lowestTaxLevel'] == mainData['HitMaxTaxLevel']]
mainData = mainData.drop(['HitMaxTaxLevel'],1)


#2
mainData['IdentityCov'] = mainData['identity'] + mainData['coverageQ']
mainData['MaxIdentityCov'] = mainData.groupby('queryid')['IdentityCov'].transform(lambda x: x.max())
#Keep the ones with the highest identity/Coverage
mainData = mainData.loc[mainData['IdentityCov'] == mainData['MaxIdentityCov']]
mainData = mainData.drop(['MaxIdentityCov', 'IdentityCov'],1)
mainData = mainData.reset_index(0, drop= True)

#3Primary annotation
df1 = mainData.copy()
df1["dbScores"]=0
df1.loc[df1.db == '16SMicrobial', 'dbScores'] = 1
df1['max_point'] = df1.groupby('queryid')['dbScores'].transform(lambda x: x.max())
df1 = df1[df1['dbScores'] == df1['max_point']]
df1 = df1.drop(['dbScores', 'max_point'],1)
#select index so I can choose secondary annotation
indexToRemove = df1.index
#reduce the number of columns
df1 = df1[["queryid","subjectid","identity","coverageQ","taxonomy","taxon1","taxon2","taxon3","taxon4","taxon5","taxon6","taxon7","lowestTaxLevel","LastCommonTaxon","ContigCount","AmbiguousHits","AmbiguousLabels"]]
df1.loc[df1.AmbiguousHits == 1, 'AmbiguousHits'] = 0
#Here I have a raw table (potentitally multi returns per query), I'll parse it with the function above

df1 = mergeMultipleDatabaseAndParse(df1)


#I'll keep the second best hit annotations. 
df2=mainData.copy()
#remove the best hits index from df2
df2.drop(indexToRemove, inplace=True)
if len(df2) != 0:
    #Reduce the number of columns before joining to best hits df
    df2["Id/Cov"]= df2["identity"].map(str) + "/" + df2["coverageQ"].map(str)
    df2 = df2[["queryid","subjectid","Id/Cov","taxonomy","taxon1","taxon2","taxon3","taxon4","taxon5","taxon6","taxon7","lowestTaxLevel","LastCommonTaxon","ContigCount","AmbiguousHits","AmbiguousLabels"]]
    #Here I have a raw table (potentitally multi returns per query), I'll parse it with the function above
    df2 = mergeMultipleDatabaseAndParse(df2)
    #Select columns and add _2 to specify the secondary annotation
    df2 = df2[['queryid', 'Id/Cov','subjectid', 'taxon7', 'AmbiguousHits','AmbiguousLabels']]
    df2 = df2.rename(columns = {'taxon7':'SubOTU_2'})
    df2 = df2.rename(columns = {'Id/Cov':'Id/Cov_2'})
    df2 = df2.rename(columns = {'subjectid':'subjectid_2'})
    df2 = df2.rename(columns = {'AmbiguousHits':'AmbiguousHits_2'})
    df2 = df2.rename(columns = {'AmbiguousLabels':'AmbiguousLabels_2'})
    data = pd.merge(df1, df2, how='outer', on='queryid', left_on=None, right_on=None,
                 left_index=False, right_index=False, sort=True,
                 suffixes=('_x', '_y'), indicator=False)
    data.fillna('-', inplace=True)

else:
    data = df1.copy()
    data["SubOTU_2"] = "-"
    data["subjectid_2"] = "-"
    data["Id/Cov_2"] = "-"
    data["AmbiguousLabels_2"] = "-"
    data["AmbiguousHits_2"] = "-"


checkDataLen=len(data)

#Last part: removing ambiguity from primers (if exists) and collapse anchors on anchor with most abundance
if args.primers == True:
    print "\n-----\nCollapsing anchor sequences based on primer ambiguity"
    primers = pd.read_csv("primers.txt", sep='\t')
    pattern = 'R|Y|S|W|K|M|B|V|D|H|N'
    checkprimers_F = primers.forward.str.contains(pattern).sum()
    checkprimers_R = primers.reverse.str.contains(pattern).sum()
    if ((checkprimers_F > 0) | (checkprimers_R > 0)):
        data_tempFinal=[]
        for i in range (0,len(primers)):
            print "Primer pair number: ", i+1
            fwd = primers.at[i, 'forward']
            #Subset the sequences containing exactly the forward primer
            fasta = pd.read_csv("anchors_fasta.tsv", sep='\t', names = ["queryid", "sequence"])
            print  "Number of anchor sequences:", len(fasta)
            print "Forward primer: ", fwd
            fasta['originalSeq'] = fasta['sequence']
            #FIRST PART: REPLACE ANY AMBIGUOUS CHARACTER POSITION IN THE FASTA BY AN "X"
            #Since the option was chosen, we assume that an ambiguous primer is in the data
            #list of the ambiguous nucleotide
            ambiChar=["R","Y","S","W","K","M","B","V","D","H","N"]
            #fwd = primer sequence (1st line; we asume there is only one forward primer)
            fwd = primers.at[i, 'forward']
            #loop over all ambiguous characters
            for char in ambiChar:
                #check the first position (there could be several) of the ambiguous character (checkAmbiPos = -1 if no ambiguous character)
                checkAmbiPos = fwd.find('%s'%char,1)
                while checkAmbiPos > -1:
                    #replace the ambiguous character position in the fasta with an X. We'll change it also in the variable to satisfy the while loop.
                    fasta['sequence'] = fasta['sequence'].apply(lambda x: change_nucl(x, checkAmbiPos))
                    fwd = change_nucl(fwd, checkAmbiPos)
                    checkAmbiPos = fwd.find('%s'%char,1)
            fastaRejected = fasta[fasta.sequence.str.startswith(fwd) == False]
            fasta = fasta[fasta.sequence.str.startswith(fwd) == True]
            #Check if some sequences were reversed in which case we have to reverse them so they all match the same strand
            if len(fastaRejected) !=0:
                fastaRejected['sequence'] = fastaRejected['originalSeq'].apply(lambda x: reverse_complement(x))
                fastaRejected['originalSeq'] = fastaRejected['sequence']
                fwd = primers.at[i, 'forward']
                for char in ambiChar:
                    #check the first position (there could be several) of the ambiguous character (checkAmbiPos = -1 if no ambiguous character)
                    checkAmbiPos = fwd.find('%s'%char,1)
                    while checkAmbiPos > -1:
                        #replace the ambiguous character position in the fasta with an X. We'll change it also in the variable to satisfy the while loop.
                        fastaRejected['sequence'] = fastaRejected['sequence'].apply(lambda x: change_nucl(x, checkAmbiPos))
                        fwd = change_nucl(fwd, checkAmbiPos)
                        checkAmbiPos = fwd.find('%s'%char,1)
                fastaToBeAdded = fastaRejected[fastaRejected.sequence.str.startswith(fwd) == True]
                if len(fastaToBeAdded) !=0:
                    fasta = fasta.append(fastaToBeAdded)
                    fastaToBeAdded = []
                    fastaRejected = []
                else:
                    fastaRejected = []
    
            print  "Number of anchor sequences containing the exact forward primer:", len(fasta)
            #it's more annoying with reverse as I have to inverse the sequences. No need to reverse the primer as it is already reversed ompared to the amplicons
            reverse = primers.at[i, 'reverse']
            print "Reverse primer: ", reverse
            #I'll also reverse the amplicon sequences
            fasta['sequenceRev'] = fasta['sequence'].apply(lambda x: x[::-1])
            for char in ambiChar:
                #check the first position (there could be several) of the ambiguous character (checkAmbiPos = -1 if no ambiguous character)
                checkAmbiPos = reverse.find('%s'%char,1)
                while checkAmbiPos > -1:
                    #replace the ambiguous character position in the fasta with an X. We'll change it also in the variable to satisfy the while loop.
                    fasta['sequenceRev'] = fasta['sequenceRev'].apply(lambda x: change_nucl(x, checkAmbiPos))
                    reverse = change_nucl(reverse, checkAmbiPos)
                    checkAmbiPos = reverse.find('%s'%char,1)
            #we'll reverse the sequences back
            fasta['sequence'] = fasta['sequenceRev'].apply(lambda x: x[::-1])
            #subset to sequences that contain the reverse primer
            reverse2=reverse_complement(reverse)[::-1]
            fasta = fasta[fasta.sequenceRev.str.startswith(reverse2) == True]
            fasta = fasta.drop(['sequenceRev'],1)
            print  "Number of anchor sequences containing the exact reverse primer:", len(fasta)
    
            #SECOND PART: REMOVE DUPLICATES AND ADJUST COUNTS
            #count the number of identical sequences
            fasta['UniqSeqCount'] = fasta.groupby(['sequence'])['queryid'].transform('count')
            #create a  column "common to" to see which anchors have the same sequence
            fasta2 = fasta.groupby('sequence')['queryid'].apply(','.join).to_frame("common_to")
            #join fasta2 to fasta 
            fasta2 = fasta2.reset_index(drop=False)
            fasta = fasta.merge(fasta2, on=['sequence'], how='outer')
            #remove the 'sequence' column and change name of orignal sequences
            fasta = fasta.drop(['sequence'],1)
            fasta = fasta.rename(columns = {'originalSeq':'sequence'})
            #join to the main output
            data_temp = data.merge(fasta, on=['queryid'], how='inner')
            if i==0:
                data_tempCheck=data_temp
            else:
                data_tempCheck=data_tempCheck.append(data_temp)
    
            #Cumulate counts from identical primers
            data_temp['CollapsedCounts'] = data_temp.groupby(['common_to'])['ContigCount'].transform('sum')
            #Now we'll check what is the most abundant and best alignment stats per new groups of anchors
            #On duplicates, we'll use the best annotation results
            #First I need to downgrade TrueUnknowns which have an identity and coverage of 100%
            pattern = 'Unknown|unknown'
            data_temp.is_copy = False
            data_temp['penalty'] = data_temp.taxonomy.str.contains(pattern)
            #we'll change True to -10 and False to 0
            data_temp['penalty'] = data_temp.penalty.astype(int) *(-10)
            data_temp['idCov'] = data_temp['identity'] + data_temp['coverageQ'] + data_temp['penalty']
            data_temp['max_idCov'] = data_temp.groupby('common_to')['idCov'].transform(lambda x: x.max())
            data_temp = data_temp[data_temp['idCov'] == data_temp['max_idCov']]
            #If there are still duplicates, we'll use the best annotation results
            data_temp['max_counts'] = data_temp.groupby('common_to')['ContigCount'].transform(lambda x: x.max())
            #Only keep the most abundant anchor as representative sequence
            data_temp = data_temp[data_temp['ContigCount'] == data_temp['max_counts']]
            #IF there are still duplicated common_to values, then we'choose the first one
            data_temp = data_temp.drop_duplicates('common_to')
            data_temp = data_temp.drop(['ContigCount','idCov','max_idCov','max_counts','penalty'],1)
            data_temp = data_temp.rename(columns = {'CollapsedCounts':'ContigCount'})
            #Extract a crib (reference anchor / similarAnchor_Primer_ )
            crib = data_temp[['queryid','common_to']]
            s = crib['common_to'].str.split(',').apply(Series, 1).stack() #this creates an object (s) with a new line every separator (here ;)
            #s has 2 indexes. the first is original data's index, and the second is s's. We'll remove the latter
            s.index = s.index.droplevel(-1) # to line up with data's index
            s.name = 'common_to'# needs a name to join
            del crib['common_to'] #remove old column
            crib = crib.join(s) #replace with new
            crib = crib.rename(columns = {'queryid':'referenceAnchor'})
            crib = crib.rename(columns = {'common_to':'anchorDemotedByPrimerAmbiguity'})
            #data_temp.to_csv("data_primerAmbiguity_%i.txt"%i,sep="\t", index=False)
            if i==0:
                data_tempFinal = data_temp
                crib_final = crib
            else:
                data_tempFinal=data_tempFinal.append(data_temp)
                crib_final=crib_final.append(crib)
    
        #Gather subsets of data
        print "Original data size: ", checkDataLen
        print "New data size: ", len(data_tempFinal)
        if checkDataLen == len(data_tempCheck):
            data = data_tempFinal.copy()
            crib_final.to_csv("crib_primerAmbiguity.txt",sep="\t", index=False)
            data_temp = []
            data_tempCheck = []
            data_tempFinal = []
            crib_final = []
            print "Primer ambiguity collapsing: done!\n"
        else:
            print "\n---\nThe number of sequences don't agree with what would be expected from the primers. Normally all sequences have been screened by the exact primers, now when screened again, it seems we have less (or more) sequences:"
            print "Original data size: ", checkDataLen
            print "Data rescreened size: ", len(data_tempFinal)
            sys.exit(1)
    else:
        fasta = pd.read_csv("anchors_fasta.tsv", sep='\t', names = ["queryid", "sequence"])
        data = data.merge(fasta, on=['queryid'], how='inner')
        data['anchorCollapsingDueToPrimer'] = "-"
else:
    fasta = pd.read_csv("anchors_fasta.tsv", sep='\t', names = ["queryid", "sequence"])
    data = data.merge(fasta, on=['queryid'], how='inner')
    data['anchorCollapsingDueToPrimer'] = "-"


if output_option == "SubOTU":
    print "\n-----------------------------------------"
    print "OPTION1: We'll separate OTUs into sub-subOTU based on common NCBI accessionID"
    #sort collapsed subjectid columns to be able to rmove all duplicates:
    data = sortCollapsedColumn(data, "queryid", "subjectid")
    #remove all subjectid duplicates
    data_uniq = data[['subjectid', 'FakeCommonLastTax', 'taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7']]
    data_uniq = data_uniq.drop_duplicates('subjectid')
    #remove the unknowns
    data_uniq = data_uniq.loc[data_uniq.taxon7 != "TrueUnknown",:]
    data_uniq = data_uniq.sort_values(by = ['subjectid']).reset_index(drop=True)
    
    #rename taxon7 otu with  a sub-level name based on NCBI accession ID
    for i in range (min(data_uniq.FakeCommonLastTax),8,1):
        if len(data_uniq.loc[data_uniq.FakeCommonLastTax == i])>0:
            data_uniq.loc[data_uniq.FakeCommonLastTax == i, 'Newtaxon'] = data_uniq.loc[data_uniq.FakeCommonLastTax == i].groupby(['taxon%i'%i])['taxon%i'%i].cumcount()+1
            data_uniq["Newtaxon"] = data_uniq["Newtaxon"].astype(str).str.replace("\.0", "")
            if i<7:
                for j in range (i+1,8,1):
                    data_uniq.loc[data_uniq.FakeCommonLastTax == i, 'taxon%i'%j] = data_uniq.loc[data_uniq.FakeCommonLastTax == i, 'taxon%i'%j]+'_'+data_uniq.loc[data_uniq.FakeCommonLastTax == i, 'Newtaxon'].astype(str)
            else:
                data_uniq.loc[data_uniq.FakeCommonLastTax == i, 'taxon%i'%i] = data_uniq.loc[data_uniq.FakeCommonLastTax == i, 'taxon%i'%i]+'_'+data_uniq.loc[data_uniq.FakeCommonLastTax == 7, 'Newtaxon'].astype(str)#Now that we have the right taxon names, we'll merge the 2 dataframe, but before we'll remove common columns names from the first one
    
    #Deal with the unknowns (they all have the same taxonomy but we don't what to collapse them into 1 OTU):
    data_unknowns = data.loc[data.taxon7 == "TrueUnknown",:]
    data_unknowns = data_unknowns.reset_index(drop=True)
    #I'll add index value +1 to all taxon7 of TrueUnknowns so I'll differentiate them
    data_unknowns.is_copy = False
    data_unknowns["newIndex"] = data_unknowns.index + 1
    data_unknowns["taxon7"] = data_unknowns["taxon7"] + "_" + data_unknowns["newIndex"].map(str)

    #remove unknowns from data (I'll add them back below)
    data = data.loc[data.taxon7 != "TrueUnknown",:]

    data = data.drop(['taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7'], 1)
    data_uniq = data_uniq.drop(['FakeCommonLastTax','Newtaxon'], 1)
    


    data_SUBOTUs= pd.merge(data, data_uniq, how='left', on='subjectid', left_on=None, right_on=None,
                 left_index=False, right_index=False, sort=True,
                 suffixes=('_x', '_y'), indicator=False)
    #Now concatenate the true unknowns with data
    frames = [data_SUBOTUs, data_unknowns] 

    data = pd.concat(frames)
    
    data_temp=[]
    data_uniq=[]
    print "Done."
    print "-----------------------------------------\n"

    
if output_option == "OTU":
    print "\n-----------------------------------------"
    print "OPTION2: We'll organize the data into common taxon7 label (i.e. OTUs)"
    #remove all subjectid duplicates
    data_uniq = data[['taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7']]
    data_uniq = data_uniq.drop_duplicates('taxon7')
    data_uniq = data_uniq.sort_values(by = ['taxon7']).reset_index(drop=True)
    #Now that we have the right taxon names, we'll merge the 2 dataframe, but before we'll remove common columns names from the first one
    data_temp = data.drop(['taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6'], 1)
    data_OTUs= pd.merge(data_temp, data_uniq, how='left', on='taxon7', left_on=None, right_on=None,
                 left_index=False, right_index=False, sort=True,
                 suffixes=('_x', '_y'), indicator=False)
    #data_OTUs=data_OTUs[['queryid', 'querylength', 'subjectid', 'subjectlength', 'identity', 'evalue', 'coverageQ', 'coverageS', 'subjectdescription', 'bitscore', 'taxonomy', 'taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7', 'MostAbundantTax1', 'MostAbundantTax2', 'MostAbundantTax3', 'MostAbundantTax4', 'MostAbundantTax5', 'MostAbundantTax6', 'MostAbundantTax7', 'lowestTaxLevel', 'ContigCount', 'AmbiguousHits', 'AmbiguousLabels']]
    data_temp=[]
    data_uniq=[]
    data = data_SUBOTUs.copy()
    print "Done."
    print "-----------------------------------------\n"


if output_option == "denovo":
    print "\n-----------------------------------------"
    print "OPTION3: de novo treatment: no collapsing at all"
    #remove all subjectid duplicates
    
    data_denovo = data.copy()
    #rename taxon7 otu with  a sub-level name based on queryid
    for i in range (min(data_denovo.lowestTaxLevel),8,1):
        if len(data_denovo.loc[data_denovo.lowestTaxLevel == i])>0:
            data_denovo.loc[data_denovo.lowestTaxLevel == i, 'Newtaxon'] = data_denovo.loc[data_denovo.lowestTaxLevel == i].groupby(['taxon%i'%i])['taxon%i'%i].cumcount()+1
            data_denovo["Newtaxon"] = data_denovo["Newtaxon"].astype(str).str.replace("\.0", "")
            if i<7:
                for j in range (i+1,8,1):
                    data_denovo.loc[data_denovo.lowestTaxLevel == i, 'taxon%i'%j] = data_denovo.loc[data_denovo.lowestTaxLevel == i, 'taxon%i'%j]+'_'+data_denovo.loc[data_denovo.lowestTaxLevel == i, 'Newtaxon'].astype(str)
            else:
                data_denovo.loc[data_denovo.lowestTaxLevel == i, 'taxon%i'%i] = data_denovo.loc[data_denovo.lowestTaxLevel == i, 'taxon%i'%i]+'_'+data_denovo.loc[data_denovo.lowestTaxLevel == 7, 'Newtaxon'].astype(str)
    #Now that we have the right taxon names, we'll merge the 2 dataframe, but before we'll remove common columns names from the first one
    data_denovo=data_denovo[['queryid', 'querylength', 'subjectid', 'subjectlength', 'identity', 'evalue', 'coverageQ', 'coverageS', 'subjectdescription', 'bitscore', 'taxonomy', 'taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7', 'MostAbundantTax1', 'MostAbundantTax2', 'MostAbundantTax3', 'MostAbundantTax4', 'MostAbundantTax5', 'MostAbundantTax6', 'MostAbundantTax7', 'lowestTaxLevel', 'ContigCount', 'AmbiguousHits', 'AmbiguousLabels']]
    data = data_SUBOTUs.copy()

    print "Done."
    print "-----------------------------------------\n"


if output_option == "SubOTU":
    workOnOptions(data,"subOTUs")
if output_option == "OTU":
    workOnOptions(data,"OTUs")
if output_option == "denovo":
    workOnOptions(data,"deNovo")


os.system('rm -f NA_* group_*')

print( '\nWhole process took : ' + str(datetime.now()-startTime) + ' h:min:ss')
print "anchorParser.py is exiting normally!"

