#!/usr/bin/python
import argparse
import pandas as pd
from pandas import *
import os
import os.path
import sys
from datetime import datetime
import subprocess


outputfile='nt_taxonomy_parsed.dat'
startTime = datetime.now()

parser = argparse.ArgumentParser(description="From a list of accession number (without version number: e.g. BA000005), it will create a file with full taxonomy.\nNote: this program will download large files, it could take some time.", epilog="Ouputfile will be: nt_taxonomy_parsed.dat.")
parser.add_argument('-i', '--inputfile', help='[REQUIRED] list of accession numbers, no headers', dest='inputfile', action='store', required=True)
parser.add_argument('-t', '--taxfolder', help='[OPTIONAL] Folder where accession2taxid.txt.gz, nodes.dmp.gz and names.dmp.gz are located. This scripts goes much faster if you provide these, it will download thenm otherwise', dest='taxfolder', action='store', required=False)
args = parser.parse_args()
inputfile = args.inputfile
taxfolder = args.taxfolder

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


# Definition of the classe Node
class Node:
    """Noeud"""
    def __init__(self):
        self.tax_id = 0       # Number of the tax id.
        self.parent = 0       # Number of the parent of this node
        self.children = []    # List of the children of this node
        self.tip = 0          # Tip=1 if it's a terminal node, 0 if not.
        self.name = ""        # Name of the node: taxa if it's a terminal node, numero if not.       
    def genealogy(self):      # Trace genealogy from root to leaf
        ancestors = []        # Initialise the list of all nodes from root to leaf.
        tax_id = self.tax_id  # Define leaf
        while 1:
            if name_object.has_key(tax_id):
                ancestors.append(tax_id)
                tax_id = name_object[tax_id].parent
            else:
                break
            if tax_id == "1":
                # If it is the root, we reached the end.
                # Add it to the list and break the loop
                ancestors.append(tax_id)
                break
        return ancestors # Return the list

# Function to find common ancestor between two nodes or more
def common_ancestor(node_list):
    global name_object
    list1 = name_object[node_list[0]].genealogy()  # Define the whole genealogy of the first node
    for node in node_list:
        list2 = name_object[node].genealogy()      # Define the whole genealogy of the second node
        ancestral_list = []                             
        for i in list1:
            if i in list2:                         # Identify common nodes between the two genealogy
                ancestral_list.append(i)                 
        list1 = ancestral_list                     # Reassing ancestral_list to list 1.
    common_ancestor = ancestral_list[0]            # Finally, the first node of the ancestra_list is the common ancestor of all nodes.
    return common_ancestor                         # Return a node







print "loading input file..."
data = pd.read_csv(inputfile, sep='\t')


print "\n-----------------------------------------"
if (taxfolder != ""):
    print "\nUnzipping taxonomy files"
    os.system('gunzip <' + taxfolder + '/accession2taxid.txt.gz > accession2taxid.txt')
    os.system('gunzip <' + taxfolder + '/names.dmp.gz > names.dmp')
    os.system('gunzip <' + taxfolder + '/nodes.dmp.gz > nodes.dmp')
else:
    if os.path.isfile("accession2taxid.txt"):
        print "No Download needed: accession2taxid.txt is already inside the folder.\n"
    else:
        print "Downloading nucl_gb.accession2taxid.gz"
        os.system('wget \"ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz\"')
        print "\nDownloading dead_nucl.accession2taxid.gz"
        os.system('wget \"ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/dead_nucl.accession2taxid.gz\"')
        print "\nUnzipping files"
        os.system('gunzip <nucl_gb.accession2taxid.gz > nucl_gb.accession2taxid')
        os.system('gunzip <dead_nucl.accession2taxid.gz > dead_nucl.accession2taxid')
        os.system('cat nucl_gb.accession2taxid dead_nucl.accession2taxid > accession2taxid.txt')
        os.system('rm nucl_gb.accession2taxid dead_nucl.accession2taxid')
    
    if os.path.isfile("names.dmp"):
        print "No Download needed: names.dmp is already inside the folder.\n"
    else:
        print "\nDownloading taxdump.tar.gz"
        os.system('wget \"ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz\"')
        print "\nUnzipping file"
        os.system('tar -zxvf taxdump.tar.gz')
    


print "Extracting your hits from accession2taxid files"
subprocess.call(['bash', '-c', 'join -1 1 -2 1 -a2 -o 2.1 1.3 -e 0 -t $"\t" <(sort -t $"\t" -k1,1 accession2taxid.txt) <(sort '+ inputfile +') > _mytaxonomy'])
print "-----------------------------------------\n"

#############################
#                           #
#   Read taxonomy files     #
#                           #
#############################

######################
# 
print "\n-----------------------------------------"
print "Load names defintion"

name_dict = {}          # Initialise dictionary with TAX_ID:NAME
name_dict_reverse = {}  # Initialise dictionary with NAME:TAX_ID

# Load  NCBI names file ("names.dmp")
name_file =  open("names.dmp","r")
while 1:
    line = name_file.readline()
    if line == "":
        break
    line = line.rstrip()
    line = line.replace("\t","")
    tab = line.split("|")
    if tab[3] == "scientific name":
        tax_id, name = tab[0], tab[1]     # Assign tax_id and name ...
        name_dict[tax_id] = name          # ... and load them
        name_dict_reverse[name] = tax_id  # ... into dictionaries
name_file.close()
print "Done."
print "-----------------------------------------\n"

######################
# 
print "\n-----------------------------------------"
print "Load taxonomy"

# Define taxonomy variable
global name_object
name_object = {}


# Load taxonomy NCBI file ("nodes.dmp")
taxonomy_file = open("nodes.dmp","r")
while 1:
    line = taxonomy_file.readline()
    if line == "":
        break
    #print line
    line = line.replace("\t","")
    tab = line.split("|")
    
    tax_id = str(tab[0])
    tax_id_parent = str(tab[1])
    division = str(tab[4])

    # Define name of the taxid
    name = "unknown"
    if tax_id in name_dict:
        name = name_dict[tax_id]
    
    if not name_object.has_key(tax_id):
        name_object[tax_id] = Node()
    name_object[tax_id].tax_id   = tax_id        # Assign tax_id
    name_object[tax_id].parent   = tax_id_parent # Assign tax_id parent
    name_object[tax_id].name     = name          # Assign name
    
    if  tax_id_parent in name_object:
        children = name_object[tax_id].children  # If parent is is already in the object
        children.append(tax_id)                  # ...we found its children.
        name_object[tax_id].children = children  # ... so add them to the parent
taxonomy_file.close()

print "Done."
print "-----------------------------------------\n"

print "\n-----------------------------------------"
if os.path.isfile("group_1.txt"):
    print "Chunk files already produced. Skipping"
else:
    print "\n-----------------------------------------"
    print "Creating chunk files"
        #record_iter = SeqIO.parse(open(inputfile),"fasta")
    itercount=0
    record_iter =  open("_mytaxonomy","r")
    for i, batch in enumerate(batch_iterator(record_iter, 10000)) :
        filename = "./group_%i.txt" % (i+1)
        handle = open(filename, "w")
        handle.write(", ".join(batch))
        #, '184\t9913\n', '187\t9913\n',
        #handle.write(str(batch))
        itercount =+ i
        handle.close()
        print "Wrote 10,000 records to %s" % (filename)
    print "Done."
    print "Total number of accession files : %i" % (itercount+1)
print "-----------------------------------------\n"


print "\n-----------------------------------------"
print "Creating taxonomy for each hit"
for i in range(1,itercount+2):
    #print "here%i" % i
    inputfile = './group_%i.txt' % (i)
    outputfile = './_nt_taxonomy_parsed_%i.dmp' % (i)
    with open(inputfile) as fasta_file:  # Will close handle cleanly
           accession_number=pd.read_csv(inputfile, sep='\t', header=None)
           accession_number.columns = ['accessionNumber', 'TaxID']
           accession_number['Taxonomy']='nan'
           accession_number['accessionNumber'] = accession_number['accessionNumber'].str.replace(', ', '')
           length = len(accession_number)
           for j in range(0,length,1):
               if j > 0 and accession_number.ix[j,1] == accession_number.ix[j-1,1]:
                   accession_number.ix[j,2] =accession_number.ix[j-1,2]
               else:
                   try:
                    tax_id_wanted = accession_number.ix[j,1].astype(str)
                    wanted_genealogy = name_object[tax_id_wanted].genealogy()
                    for tax_id in wanted_genealogy:
                        accession_number.ix[j,2] = name_object[tax_id].name + ";" + accession_number.ix[j,2]  
                   except KeyError:
                        continue

    accession_number['Taxonomy'] = accession_number['Taxonomy'].str.replace(';nan', '')
    accession_number['Taxonomy'] = accession_number['Taxonomy'].str.replace('root;', '')

    #accession_number = accession_number[['accessionNumber','Taxonomy']]
    #print accession_number.head(n=3)
    accession_number.to_csv(outputfile, sep='\t', index=False)        

os.system('cat _nt_taxonomy_parsed_*.dmp >'+ "nt_taxonomy_parsed.dat")
os.system('sed -i \'0,/accessionNumber/! {/accessionNumber/d}\' nt_taxonomy_parsed.dat')
os.system('rm -f _nt_taxonomy_parsed*.dmp group_*.txt _mytaxonomy citations.dmp delnodes.dmp division.dmp gc.prt gencode.dmp merged.dmp names.dmp nodes.dmp readme.txt taxdump.tar.gz dead_nucl.accession2taxid.gz nucl_gb.accession2taxid.gz accession2taxid.txt')

print "Well done dude!"
print( '\nProcess took : ' + str(datetime.now()-startTime) + ' h:min:ss')