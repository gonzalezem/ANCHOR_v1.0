#!/bin/bash
set -e

#fastaToCsv.sh inputfile outputfile

awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' $1 > $2
sed -i "s/^>//" $2