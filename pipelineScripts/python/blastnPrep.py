#!/usr/bin/python

from Bio import SeqIO
#from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils import GC
import os
from datetime import datetime
import argparse
import pandas as pd
from pandas import *
import io

startTime = datetime.now()

parser = argparse.ArgumentParser(description="Returns a tab separated table with sequences identification, length, and description of a fasta file : prefix_ParsedFasta.txt. If Trinity input 2 files will be generated, see below. Note : don't use this script to parse nr protein database fasta file. Use the script fastaparser_for-nr-fasta-file.py instead.")
parser.add_argument('-i', '--ifile', help='[REQUIRED] Fasta file', dest='inputfile', action='store', required=True)
parser.add_argument('-p', '--path', help='[REQUIRED] Path where you want the folders and files to be installed', dest='path', action='store', required=True)
parser.add_argument('-dp', '--databasePath', help='[REQUIRED] path to database', dest='dirdb', action='store', required=True)
parser.add_argument('-dn', '--databaseName', help='[REQUIRED] Database name', dest='dbName', action='store', required=True)
parser.add_argument('-cc', '--chunkCounts', help='[OPTIONAL] Number of sequences inside each chunks. Default is 500. Note: We have hardcoded Swissprot to have 2000 sequences chunks as it is fast.).', dest='chunkCounts', default = 500, action='store', type=float)

args = parser.parse_args()

inputfile = args.inputfile
mypath = args.path
dbName = args.dbName
chunkCounts = args.chunkCounts
dirdb = args.dirdb


def SimpleFastaParser(handle):
    """Generator function to iterator over Fasta records (as string tuples).

For each record a tuple of two strings is returned, the FASTA title
line (without the leading '>' character), and the sequence (with any
whitespace removed). The title line is not divided up into an
identifier (the first word) and comment or description.

>>> for values in SimpleFastaParser(open("Fasta/dups.fasta")):
... print(values)
('alpha', 'ACGTA')
('beta', 'CGTC')
('gamma', 'CCGCC')
('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
('delta', 'CGCGC')

"""
    #Skip any text before the first record (e.g. blank lines, comments)
    while True:
        line = handle.readline()
        if line == "":
            return # Premature end of file, or just empty?
        if line[0] == ">":
            break

    while True:
        if line[0] != ">":
            raise ValueError(
                "Records in Fasta files should start with '>' character")
        title = line[1:].rstrip()
        lines = []
        line = handle.readline()
        while True:
            if not line:
                break
            if line[0] == ">":
                break
            lines.append(line.rstrip())
            line = handle.readline()

        #Remove trailing whitespace, and any internal spaces
        #(and any embedded \r which are possible in mangled files
        #when not opened in universal read lines mode)
        yield title, "".join(lines).replace(" ", "").replace("\r", "")

        if not line:
            return # StopIteration

    assert False, "Should not reach this line"

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

chunkname = inputfile.split('/')[-1]
chunkname = chunkname.split('.fasta')[0]





######################################################FASTA CHUNKS#################################################################
os.system('rm -rf %s/fastachunks %s/scripts %s/blastResults_%s' %(mypath,mypath,mypath,dbName))
os.system('mkdir -p %s/fastachunks %s/scripts %s/blastResults_%s' %(mypath,mypath,mypath,dbName))



record_iter = SeqIO.parse(open(inputfile),"fasta")


itercount = 0

print "Total number of fasta chunks chosen by user: %i" % (chunkCounts)


print "Split %s for %s databases" % (inputfile,dbName)
print record_iter
for i, batch in enumerate(batch_iterator(record_iter, chunkCounts)) :
    filename = "%s_%i.fasta" % (chunkname,i+1)
    handle = open(filename, "w")
    count = SeqIO.write(batch, handle, "fasta")
    os.system('mv %s_%i.fasta %s/fastachunks'% (chunkname,i+1,mypath))
    itercount =+ i
    handle.close()
print "Total number of fasta chunks : %i" % (itercount+1)
    


###############################################shell_scripts maker : nr, trembl, 





itercount = itercount +1
    


print "-----\nCreate launch scripts for BLASTn"
with io.FileIO("blastn_vs_%s.sh"%dbName, "w") as file:
    file.write("#!/bin/bash\nset -e\n")
    file.write("dirdb=\"%s\"\n" %(dirdb))
    file.write("for i in {1..%i..1}\ndo\n\techo \"Blasting ${i}_fasta against %s database\"\n\tdate\n\t" %(itercount,dbName))
    file.write("blastn -query %s/fastachunks/%s_${i}.fasta -db ${dirdb}/%s_index/%s -out %s/blastResults_%s/output_%s_${i}.txt -evalue 1e-10 -num_threads 12 -word_size 30 -perc_identity 80 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle qseq sseq'\n\tdate\n" %(mypath,chunkname,dbName,dbName,mypath,dbName,dbName))
    file.write("done\n")


os.system('mv blastn_vs_%s.sh %s/scripts' %(dbName,mypath))







print "\nOUTPUT INFORMATION"
print "---------------------"
print( '\nProcess took : ' + str(datetime.now()-startTime) + ' h:min:ss')
