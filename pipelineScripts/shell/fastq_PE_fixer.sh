#!/bin/bash
set -e

: <<'END'
This script order the reads by name and remove missing paired ends. Read names needs to end by /1 and /2
Example:
fastq_PE_fixer.sh Reads_Forward_R1.fastq  Reads_Forward_R2.fastq

The output are 2 files:
Reads_Forward_R1_ordered.fastq
Reads_Forward_R2_ordered.fastq

END

outputname1=$(basename "$1" | cut -d. -f1)
outputname2=$(basename "$2" | cut -d. -f1)

#unzip inputfile in case it is gzipped
if [[ $1 = *.gz ]]; then
	inputfile1=$(basename $1 .gz)
	gunzip -c $1 > $inputfile1
else
	inputfile1=$1
fi

if [[ $2 = *.gz ]]; then
	inputfile2=$(basename $2 .gz)
	gunzip -c $2 > $inputfile2
else
	inputfile2=$2
fi


#transform the file to csv
cat $inputfile1 | awk  -F'\t' '{ printf "%s%s", $0, (NR%4==1 || NR%4==2 || NR%4==3 ? FS : RS)}' > __read1
cat $inputfile2 | awk  -F'\t' '{ printf "%s%s", $0, (NR%4==1 || NR%4==2 || NR%4==3 ? FS : RS)}' > __read2

cut -f1 __read1 | sed "s/ *\/1//" | sort | uniq | sort -t $' ' -k 1,1 > __read1Names
cut -f1 __read2 | sed "s/ *\/2//" | sort | uniq | sort -t $' ' -k 1,1 > __read2Names

join __read1Names __read2Names > __goodReads

sed "s/$/ \/1/" __goodReads | sort -t $' ' -k 1,1 > __goodR1
sed "s/$/ \/2/" __goodReads | sort -t $' ' -k 1,1 > __goodR2

#The order is tricky so the last sort is on the space (assuming the reads are like @6401_16S_10002 \1) and not the tab
join -1 1 -2 1 -t $'\t' <(sort -t $'\t' -k 1,1 __read1) <(sort -k 1 __goodR1) | sort -t $' ' -k 1,1 | sed "s/\t/\n/g" > ${outputname1}_ordered.fastq
join -1 1 -2 1 -t $'\t' <(sort -t $'\t' -k 1,1 __read2) <(sort -k 1 __goodR2) | sort -t $' ' -k 1,1 | sed "s/\t/\n/g" > ${outputname2}_ordered.fastq

rm -f __read[12] __read[12]Names __goodReads __goodR[12]
echo -e "${outputname1}.fastq and ${outputname1}.fastq are now cleaned"



