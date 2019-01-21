#!/usr/bin/python

from Bio import SeqIO
#from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqUtils import GC
import os
from datetime import datetime
import argparse
import pandas as pd
from pandas import *

startTime = datetime.now()

parser = argparse.ArgumentParser(description="Returns a tab separated table with sequences identification, length, and description of a fasta file : prefix_ParsedFasta.txt. If Trinity input 2 files will be generated, see below. Note : don't use this script to parse nr proteine database fasta file. Use the script fastaparser_for-nr-fasta-file.py instead.")
parser.add_argument('-i', '--ifile', help='[REQUIRED] Fasta file', dest='inputfile', action='store', required=True)
parser.add_argument('-p', '--prefix', help='[REQUIRED] Prefix is a string that will inserted at the beginning of the outputfile name. Input file example : -p Swissprot',dest='prefix', action='store', required=True)
parser.add_argument('-t', '--trinity', help='[OPTIONAL] Subject sequences are Trinity subject sequences (e.g comp1_c0_seq1). If -t is applied, two files will be generated: one descriptive (prefix_ContigDescrp.txt) with 4 columns (contigid, ContigLength, description, and GCContent) and another one that will map isoforms to compounds/genes (prefix_Contig_Coumpound_List.txt).', dest='trinity', default=False, action='store_true')

args = parser.parse_args()

inputfile = args.inputfile
prefix = args.prefix
Trinity = args.trinity
#inputfile = '/home/gonzaleze/Programs/ncbi-blast-2.2.27+/db/nr/nr.fasta'
#inputfile = '/home/gonzaleze/Desktop/Dropbox/Linux_Scripts/blastalignmentproportion_workspace/test_nr.fasta'

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

if Trinity == True:
    outputfile1 = ( prefix + '_ContigDescrp_with_GC.txt')
    outputfile2 = ( prefix + '_Contig_Coumpound_List.txt')
    outputfile3 = ( prefix + '_ContigDescrp.txt')
    with open(inputfile) as fasta_file:  # Will close handle cleanly
        identifier = []
        length = []
        description = []
        gccontent = []
        for title, sequence in SimpleFastaParser(fasta_file):
            identifier.append(title.split(None, 1)[0])  # First word is ID
            length.append(len(sequence))
            gccontent.append(GC(sequence))
            description.append("No Description") # Description is "No Description"
    #ContigDescrp = DataFrame(dict(subjectid = Series(identifier, name = 'subjectid'), subjectlength = Series(length, name = 'subjectlength'))).set_index(['subjectid'])
    ContigDescrp = DataFrame(dict(Contigid = Series(identifier, name = 'Contigid'), ContigLength = Series(length, name = 'ContigLength'), GCContent = Series(gccontent, name = 'GCContent'), Description =Series(description, name = 'Description') )).set_index(['Contigid'])
    #print ContigDescrp
    ContigDescrp=ContigDescrp[["ContigLength", "Description", "GCContent"]]
    ContigDescrp.to_csv(outputfile1, sep='\t', index=True)
    ContigDescrp = ContigDescrp.drop('GCContent',1)
    ContigDescrp.to_csv(outputfile3, sep='\t', index=True)
    #Getting another column : gene id
    ContigDescrp['compoundid'] = [x.split("_i")[0] for x in ContigDescrp.index]
    Contig_to_Compounds = ContigDescrp.compoundid

    #Contig_to_Compounds.reset_index().set_index(['compoundid'])


    Contig_to_Compounds.to_csv("temp.txt", sep='\t', index=True)
    os.system('awk -F"\t" \'{print $2 "\t" $1}\' temp.txt > '+ outputfile2)
    os.system('rm temp.txt')
    print "\nOUTPUT INFORMATION"
    print "---------------------"
    print ('Output description file (with GC content): ' + str(outputfile1))
    print ('Output description file : ' + str(outputfile3))
    print ('Output isoform/compound maping file : ' + str(outputfile2))
    print('\nProcess took : ' + str(datetime.now()-startTime) + ' h:min:ss')


else:
    outputfile = ( prefix + '_ParsedFasta.txt')
    record_iter = SeqIO.parse(open(inputfile),"fasta")
    print record_iter
    for i, batch in enumerate(batch_iterator(record_iter, 500000)) :
        filename = "group_%i.fasta" % (i+1)
        handle = open(filename, "w")
        count = SeqIO.write(batch, handle, "fasta")
        itercount =+ i
        handle.close()
        print "Wrote %i records to %s" % (count, filename)
    print "Total number of new fasta files : %i" % (itercount+1)
    
    for i in range(1,itercount+2):
        #print "here%i" % i
        inputfile = './group_%i.fasta' % (i)
        with open(inputfile) as fasta_file:  # Will close handle cleanly
            identifier = []
            length = []
            description = []
            for title, sequence in SimpleFastaParser(fasta_file):
                identifier.append(title.split(None, 1)[0])  # First word is ID
                length.append(len(sequence))
                test = title.split(None, 1)[-1] == title.split(None, 1)[0]
                #print test
                if test is True:
                    description.append("No Description") # Description is "No Description"
                else:
                    description.append(title.split(None, 1)[1])  # First word is ID
        Parsedfasta = DataFrame(dict(subjectid = Series(identifier, name = 'subjectid'), subjectlength = Series(length, name = 'subjectlength'), subjectdescription =Series(description, name = 'subjectdescription') )).set_index(['subjectid'])
        Parsedfasta=Parsedfasta[["subjectlength", "subjectdescription"]]
        Parsedfasta.to_csv(('ParsedFasta_%i.txt' % i), sep='\t', index=True)
    os.system('cat ParsedFasta_*.txt >'+ outputfile)
    os.system('sed -i \'0,/ubjectid/! {/ubjectid/d}\' ' + outputfile)
    os.system('rm ParsedFasta_*.txt group_*.fasta')

    print "\nOUTPUT INFORMATION"
    print "---------------------"
    print ('Output file : ' + str(outputfile))
    print( '\nProcess took : ' + str(datetime.now()-startTime) + ' h:min:ss')
