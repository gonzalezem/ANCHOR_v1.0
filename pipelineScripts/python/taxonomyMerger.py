#!/usr/bin/python
import sys
from datetime import datetime
import argparse
import pandas as pd
from pandas import *
import subprocess

def HitsPerQuery(df):
    #calculating the number of sequences per query 
    init_number = len(df)
    df = df.drop_duplicates('queryid')
    final_number = len(df)
    return float(init_number) / float(final_number)



def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
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





startTime = datetime.now()

parser = argparse.ArgumentParser(description='inputfile: myfile_part2_parsed.txt \noutput file : myfile_part3_parsed.txt')
parser.add_argument('-i', '--inputfile', help='[REQUIRED] Output file part2_bitscoreParser.py. It should have 20 columns: queryid|subjectid|identity|alignmentlength|mismatches|gapopens|qstart|qend|sstart|send|evalue|bitscore|querylength|querydescription|Q_AlignmentLength|coverageQ|subjectlength|subjectdescription|S_AlignmentLength|coverageS', dest='inputfile', action='store', required=True)
parser.add_argument('-t', '--taxonomyfile', help='Taxonomy cured file. It should contain 3 columns: AcessionID(without extension) | cured taxonomy. No header is assumed',dest='taxonomyfile', action='store', required=True)
parser.add_argument('-id', '--identity', help='[OPTIONAL] Filter input by identity threshold. Removes any record that has an identity lower than the threshold. Example : -id 99.', dest='identity', default = 0.1, action='store', type=float)
parser.add_argument('-c', '--coverage', help='[OPTIONAL] Filter input by alignment coverage on queryid. Removes any record that has an identity lower than the threshold. Example : -c 99.', dest='coverage', default = 0.1, action='store', type=float)


args = parser.parse_args()
inputfile = args.inputfile
taxonomyfile = args.taxonomyfile
identity = args.identity
coverage = args.coverage


print "\nloading taxonomy files."

curatedTax = pd.read_csv(taxonomyfile, sep='\t', header=None)
try:
    len(curatedTax.columns) == 2
    curatedTax.columns = ['accession', 'taxonomy']
    curatedTaxCheck = curatedTax.drop_duplicates('accession')
except:
    print "\n***************"
    print "Error detected!"
    print "Problem with :", taxonomyfile
    print "Make sure this file has 2 columns : accession id and cured taxonomy (no header)"
    print "***************"
    sys.exit(1)


print "\n-----------------------------------------"
print "Split full taxonomy into 7 taxons"
check_number_of_hits = curatedTax.shape[0]
taxa1 = curatedTax['taxonomy'].apply(lambda x: pd.Series(x.split(';',)))
taxa1.columns = ['taxon1', 'taxon2', 'taxon3', 'taxon4', 'taxon5', 'taxon6', 'taxon7']
curatedTax_test = curatedTax.join(taxa1, how='inner')
if curatedTax_test.shape[0] ==  check_number_of_hits:
    curatedTax = curatedTax_test.copy()
else:
    print "\n***************"
    print "Error detected!"
    print "{:,}".format(curatedTax_test.shape[0]), "records found a taxonomy hit in your file, instead of", "{:,}".format(check_number_of_hits), "."
    print "Check file: _curatedTax_problem.txt"
    curatedTax = curatedTax.reset_index(drop=False)
    curatedTax_problem = pd.merge(curatedTax, curatedTax, on='taxref', how='left', indicator=True)#, lsuffix='_a', rsuffix='_b')    print "\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    curatedTax_problem = curatedTax_problem.loc[curatedTax_problem['_merge'] != 'both']
    print curatedTax_problem.head(n=2)
    curatedTax_problem.to_csv("_curatedTax_problem.txt", sep = '\t', index = False)
    sys.exit(1)


print "\n-----------------------------------------"
print "Loading BLASTn file and split it\n"
itercount=0
f = open(inputfile, "r")
for i, batch in enumerate(batch_iterator(f, 500000)) :
    filename = "./group_%i.txt" % (i+1)
    handle = open(filename, "w")
    handle.write("".join(batch))
    itercount =+ i
    handle.close()
    print "Wrote 500,000 records to %s" % (filename)
f.close()
print "Done."
print "Total number of chunk files : %i" % (itercount+1)
print "-----------------------------------------\n"



print "\n-----------------------------------------"
print "Creating taxonomy for each hit"
for i in range(1,itercount+2):
    inputfile = 'group_%i.txt' % (i)
    outputfile = '_chunk_parsed_%i.temp' % (i)

    print "loading:", inputfile
    if i==1:
        data = pd.read_csv(inputfile, sep='\t', skiprows=1, header=None)
    else:
        data = pd.read_csv(inputfile, sep='\t', header=None)


    try:
        len(data.columns) == 18
        data.columns = ["queryid", "subjectid", "identity", "alignmentlength", "mismatches", "gapopens", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "querylength", "querydescription", "coverageQ", "subjectlength", "coverageS", "subjectdescription"]
    except:
        print "\n***************"
        print "Error detected!"
        print "Problem with :", inputfile
        print "Make sure that you have 18 columns : queryid|subjectid|identity|alignmentlength|mismatches|gapopens|qstart|qend|sstart|send|evalue|bitscore|querylength|querydescription|coverageQ|subjectlength|coverageS|subjectdescription"
        print "***************"
        sys.exit(1)




    print '\n'
    print "\nCHUNK INFORMATION"
    print "-------------------------------"
    print ('Chunk name: ' + str(inputfile))
    print ('Database parsed fasta file: ' + str(taxonomyfile))
    print "Chunk has: ", "{:,}".format(data.shape[0]), " records"
    print "Chunk has: ", "{:,}".format(data.drop_duplicates("queryid").shape[0]), "unique amplicons"
    print "Taxonomy file has: ", "{:,}".format(curatedTax.shape[0]), " records"
    print "Average number of hits per sequence : ", format(HitsPerQuery(data), "8.2f")
    print "-------------------------------\n"


    if identity == 0.1:
        pass
    else:
        print "\n-----------------------------------------"
        print "Removing hits with identity <", identity
        data = data[data['identity']>=identity] #Keep rows that have 'identity' value lower or identical to 'identity' threshold
        data = data.reset_index(drop=True) #reset index
        print "Chunk now has: ", "{:,}".format(data.shape[0]), " records"
        print "Taxonomy file has: ", "{:,}".format(curatedTax.shape[0]), " records"
        print "Average number of hits per sequence : ", format(HitsPerQuery(data), "8.2f")
        print "-----------------------------------------\n"
    
    if coverage == 0.1:
        pass
    else:
        print "\n-----------------------------------------"
        print "Removing hits with alignment coverage of amplicon <", coverage
        data = data[data['coverageQ']>=coverage] #Keep rows that have 'coverage' value lower or identical to 'coverage' threshold
        data = data.reset_index(drop=True) #reset index
        print "Chunk now has: ", "{:,}".format(data.shape[0]), " records"
        print "Taxonomy file has: ", "{:,}".format(curatedTax.shape[0]), " records"
        print "Average number of hits per sequence : ", format(HitsPerQuery(data), "8.2f")
        print "-----------------------------------------\n"
    



    print "\n-----------------------------------------"
    print "Creating accession id column"
    #FROM FJ968741.1
    #TO FJ968741
    data['accession'] = [x.split(".",1)[0] for x in data['subjectid']]
    data['accession_Full'] = data['subjectid']
    curatedTax['accession'] = [x.split(".",1)[0] for x in curatedTax['accession']]
    print "Done!"
    print "-----------------------------------------\n"



    print "\n-----------------------------------------"
    print "merge blastn output with taxonomy file"
    df1=merge(data, curatedTax, left_on='accession', right_on='accession', how='inner')
    dataCheck = data.drop_duplicates("queryid")
    df1Check = df1.drop_duplicates("queryid")
    if len(dataCheck) ==  len(df1Check):
        data = df1.copy()
        #in case of duplicates
        data = data.drop_duplicates()
        print "Chunk now has: ", "{:,}".format(data.shape[0]), " records"
        print "-----------------------------------------\n"
    else:
        print "\n***************"
        print "Error detected!"
        print "Check the columns accession on:\n", inputfile, "(", "{:,}".format(data.shape[0]), "records)\nand\n",taxonomyfile,"(", "{:,}".format(curatedTax.shape[0]), "records ).\nMerging these two created a table with", "{:,}".format(len(df1)), "records."
        print "--------------------------------------\n\nIF THE NUMBER IS LOWER: look at the _merge column in _problematic_merging_lower.txt table."
        df1=merge(data, curatedTax, left_on='accession', right_on='accession', how='outer', indicator=True)
        df1 = df1[df1["_merge"] == "left_only"]
        df1.to_csv("_problematic_merging_lower.txt", sep='\t', index=False)
        print df1.tail(n=1)
        print "--------------------------------------\n\nIF THE NUMBER IS HIGHER: look at the duplicates columns in _problematic_merging_higher.txt table."
        df2=pd.concat(g for _, g in data.groupby("queryid") if len(g) > 1)
        df2.to_csv("_problematic_merging_higher.txt", sep='\t', index=False)
        print df2.head(n=5)
        sys.exit(1)

    data = data.sort_values(by = ['queryid', 'bitscore', 'identity', 'coverageQ', 'evalue'], ascending = [True, False, False, False, True]).reset_index(drop=True)


    print "\n-----------------------------------------"
    print "Removing sstart and send columns"
    data = data.drop(['sstart', 'send'],1)
    data = data.drop_duplicates()
    print "Chunk now has: ", "{:,}".format(data.shape[0]), " records"
    print "Taxonomy file has: ", "{:,}".format(curatedTax.shape[0]), " records"
    print "Average number of hits per sequence : ", format(HitsPerQuery(data), "8.2f")
    print "-----------------------------------------\n"

    data.to_csv(outputfile, sep='\t', index=False)
    print "Next chunk!\n\n"





print "\n     ASSEMBLING FINAL OUTPUT:     "
print "-------------------------------"

outputfile = 'blastnTaxMerged.txt'


command = 'head -n1 _chunk_parsed_1.temp > __header && cat _chunk_parsed_*.temp | sort | uniq > __thereThere && cat __header __thereThere>' + outputfile + ' && sed -i \'0,/queryid/! {/queryid/d}\' '+ outputfile + ' && rm -f _chunk_parsed_*.temp group_*.txt __thereThere __header'
command_run = subprocess.call(command, shell=True)
if command_run != 0:
     exit(1)



print ('Output file: ' + str(outputfile))
print( '\nWhole process took : ' + str(datetime.now()-startTime) + ' h:min:ss')
print "taxonomyMerger.py exiting normally"

