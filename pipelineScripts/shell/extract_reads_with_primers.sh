#!/bin/bash
set -e

: <<'END'
This script detects the presence of a sequence in a fastq file (or fastq.gz) and remove either the content of it or a part of it. 
It needs 2 or 3 parameters:
remove_primers_from_fastq.sh inputfile.fastq primer_sequences number_of_bases_to_trim
or
remove_primers_from_fastq.sh inputfile.fastq primer_sequences

If the 3rd parameter is empty, the full primer will be removed from the sequence.

The output are 2 files:
inputfile_noPrimerFound.fastq (sequences were no primer were found)
inputfile_withPrimers.fastq (sequences were the primer was found)

END

#Extract prefix from the input file for the output file
outputname=$(basename "$1" | cut -d. -f1)

#unzip input file
if [[ $1 = *.gz ]]; then
	inputfile=$(basename $1 .gz)
	gunzip -c $1 > $inputfile
else
	inputfile=$1
fi

#check whether parameter 3 is empty
if [[ -z "$3" ]]; then
	lengthToRemove=${#2}
else
	lengthToRemove=$3
fi

#Modify IUPAC ambiguity codes
primerSeq=$(echo $2 | sed -e "s/H/\[ACT\]/g" | sed -e "s/V/\[ACG\]/g" | sed -e "s/N/\[ACGT\]/g" | sed -e "s/W/\[AT\]/g" | sed -e "s/R/\[AG\]/g" | sed -e "s/Y/\[CT\]/g" | sed -e "s/S/\[GC\]/g" | sed -e "s/K/\[GT\]/g" | sed -e "s/M/\[AC\]/g" | sed -e "s/B/\[CGT\]/g" | sed -e "s/D/\[AGT\]/g")


#transform the file to tsv
cat $inputfile | awk  -F'\t' '{ printf "%s%s", $0, (NR%4==1 || NR%4==2 || NR%4==3 ? FS : RS)}' > __full_file
inputLength=$(wc -l __full_file | cut -d" " -f1)

#check if the primer is in the sequence
grep $'\t'"$primerSeq" __full_file > __with_primers | true  


if [[ -s __with_primers ]];
then
	primedLength=$(wc -l __with_primers | cut -d" " -f1)
	#separate the seqs without primers in another file
	join  -1 1 -2 1 -t $'\t' -v1 <(sort -t $'\t' -k1,1 __full_file) <(sort -t $'\t' -k1,1 __with_primers) > __without_primers
	untouchedLength=$(wc -l __without_primers | cut -d" " -f1)
	#transform it back to fastq
	sed "s/\t/\n/g" __without_primers > ${outputname}_noPrimerFound.fastq
	rm -f __full_file __without_primers
	#remove the primer
	awk -F'\t' -v var="$lengthToRemove" '{$2=substr($2,var+1);$4=substr($4,var+1)}1' OFS='\t' __with_primers  | sed "s/\t/\n/g" > ${outputname}_withPrimers.fastq
	rm -f __with_primers
	#remove the unzipped file in case the inputfile was gzipped 
	if [[ $1 = *.gz ]]; then
		rm -f $inputfile
	fi
else
	sed "s/\t/\n/g" __full_file > ${outputname}_noPrimerFound.fastq
	rm -f __full_file __with_primers
fi


