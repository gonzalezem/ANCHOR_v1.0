#!/bin/bash
set -e

:<<'END'
Uses as inpufile Mothur's unique.seqs() output : fastaname.trim.contigs.good.names 
The goal is to remove the comma separated list in column two and have only one sequence per linel in column 2 instead of a list.

Outputfile: inputfileName_parsed.txt")
bash mothur_uniq_report_parser fastaname.trim.contigs.good.names outputname
END

echo -e "\nmothur_uniq_report_parser.sh"
echo -e "-------------------------------"
echo -e "Inputfile: $1"
echo -e "Outputfile: $2"
echo -e "-------------------------------"


sed "s/,/\n\t/g" $1 | awk 'NF==1{print p "\t" $1; next} {p=$1} 1' > $2

echo -e "mothur_uniq_report_parser.sh is exiting normally"

