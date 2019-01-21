#!/bin/bash
set -e

: <<'END'
1. trim all reads to first 50 nucl
2. look for primers (ambiguous nuclotides accepted) in the 50 nucl
3.When found, remove all nucl before the primer in original reads

Note1: pandas (python) is used

Note2: Expects primers.txt in the folder and expects a header. Example:
forward	reverse
ATTYCGGCGRCTGCTGG	CTACGGGAGGCAGCAG
CTACGGGAGGCAGCAG	ATTACCGCGGCTGCTGG



Example:
If primer is: ATTACCGCGGCTGCTGG

from 
@P01_post_CRSsNP_NR_1 /1
ATAGATTACCGCGGCTGCTGGCACGGAATTAGCCGGTCCTTATTCATGCGGTACCTGCAATGAAGGACACGTCCCTCACTTTATCCCCGCATAAAAGCAGTTTACAACC
+
AAA1AFFFFFADGGGGGGGGGGHH?0/FFHFHHH?A/AAFHHHFA2BFGGCGGEFFGFHHHEHHHHGFGH?>/>EEGHAHHBE@FFGG?@EAEHHGFFH0BBEGFGFFC

to
@P01_post_CRSsNP_NR_10 /1
ATTACCGCGGCTGCTGGCACGGAATTAGCCGGTCCTTATTCATGCGGTACCTGCAATGAAGGACACGTCCCTCACTTTATCCCCGCATAAAAGCAGTTTACAACC
+
AFFFFFADGGGGGGGGGGHH?0/FFHFHHH?A/AAFHHHFA2BFGGCGGEFFGFHHHEHHHHGFGH?>/>EEGHAHHBE@FFGG?@EAEHHGFFH0BBEGFGFFC

command line:
bash uncover_primers_from_reads.sh forward_R1.fastq.gz (or fastq)
END


dir0="YourMainPath"
dirpipe="${dir0}/pipelineScripts/shell"
dirMeta="${dir0}/metadata"
dirreads="${dir0}/raw_reads"
successDir="${dir0}/run/successfulRuns"
primers="${dirMeta}/primers.txt"
primerSelectionBypass=$(grep "^primerSelectionBypass" $1 | cut -d"=" -f2)





echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Main directory: $dir0"
echo -e "Pipeline scripts directory: $dirpipe"
echo -e "Metadata directory: $dirMeta"
echo -e "Raw reads directory: ${dirreads}"
echo -e "Bypass primer selection: $primerSelectionBypass"
echo -e "-----------------------------------\n\n"


echo -e "Checking input before running the script"


FILE="${dirpipe}/fastq_PE_fixer.sh"
if [ ! -e "$FILE" ];
then
   echo -e "Shell script $FILE does not exist. You need it to run. Find it and place it here: ${dirpipe}/" >&2
   exit 1
fi


if [ ! -d "$dirreads" ]; then
  echo -e "Couldn't find the reads directory (${dirreads}). Did you run part1_remove_phix.pbs.sh before running this script?"
fi

cd ${dirreads}
#unzip raw reads if needed
checkFastqs=$(ls -1 *.* | rev | cut -d"." -f1 | rev | grep "fastq" | wc -l)
if [ ${checkFastqs} -eq 0 ]; then
	for i in *.gz
	do
		removeFastq="YES"
		newfile=$(basename $i .fastq.gz)
		gunzip -c $i > $newfile.fastq
	done
fi

mkdir -p ${dir0}/run
cd ${dir0}/run

if [ "${primerSelectionBypass}"  == "YES" ]; then
	echo -e "By-passing primer removal step."
	mkdir -p ${dir0}/run/fastqsWithPrimers
	cd ${dirreads}
	#remove fastq in raw reads directory if they were zipped
	if [ "${removeFastq}"  == "YES" ]; then
		mv ${dirreads}/*.fastq ${dir0}/run/fastqsWithPrimers
	else
		cd ${dirreads}
		cp *.fastq ${dir0}/run/fastqsWithPrimers
		for i in *.fastq; do echo -e "zipping --> ${i}";gzip $i; done
	fi
	mkdir -p ${successDir}
	mkdir -p ${dir0}/run/Summary
	touch ${successDir}/primerSelection.ok
	echo -e "\n---\nNo primer were removed at the user's choice. Script is exiting normally!"
	exit 0
fi

mkdir -p ${dir0}/fastqsWithPrimers ${dir0}/fastqsWithoutExactPrimers
cd ${dir0}/fastqsWithPrimers


FILE="${primers}"
if [ ! -e "$FILE" ];
then
   echo -e "Primer file $FILE does not exist. You need it to run. Find it and place it here: ${primers}" >&2
   exit 1
fi


cd ${dir0}/fastqsWithPrimers
ls -1 ${dirreads}/*.fastq | sort > __readsList
echo -e "Read_Name\tTotal_reads\tprimers_found\tNo_primer_found\tpercentage_kept_reads" > primerSelection_summary.log
rm -f _final_good_1 _final_good_2
while read myfastq
do
	inputfile=${myfastq}

	#Extract prefix from the input file for the output file
	outputname=$(basename "${myfastq}" | cut -d. -f1)

	#check if file has already been run
	FILE="${dir0}/fastqsWithPrimers/${outputname}_withPrimers.fastq"
	if [ ! -e "$FILE" ]; then	
		#Modify IUPAC ambiguity codes
		awk 'FNR>1' ${primers} | sed -e "s/H/\[ACT\]/g" | sed -e "s/V/\[ACG\]/g" | sed -e "s/N/\[ACGT\]/g" | sed -e "s/W/\[AT\]/g" | sed -e "s/R/\[AG\]/g" | sed -e "s/Y/\[CT\]/g" | sed -e "s/S/\[GC\]/g" | sed -e "s/K/\[GT\]/g" | sed -e "s/M/\[AC\]/g" | sed -e "s/B/\[CGT\]/g" | sed -e "s/D/\[AGT\]/g" | sed "s/\t/\n/g" | sort | uniq > __primers
		
		
		#transform the file to tsv
		cat $inputfile | awk  -F'\t' '{ printf "%s%s", $0, (NR%4==1 || NR%4==2 || NR%4==3 ? FS : RS)}' > __full_file
		inputLength=$(wc -l __full_file | cut -d" " -f1)
		
		cp __full_file _not_found
		while read myPrimer
		do
			echo -e "Primer : ${myPrimer}\nRead: ${myfastq}"
			#We'll cut the 50 first bases and look for the primer there
			cut -f2 _not_found | sed "s/\(^...................................................\)\(..*\)/\1/" > _temp1
			paste _temp1 <(cut -f1 _not_found) > _temp2
			#Now we have 2 columns in _temp2: truncated sequence and sequence ID
			#Let's extract the sequences starting with the primer
			grep $'\t'"${myPrimer}" _not_found >> _final_good_1 | true
			
			#Now we'll look for the primer anywhere in the truncated sequence
			cat _temp2 | grep -v "^${myPrimer}" | grep "${myPrimer}" > _Potgood2 | true
			if [[ -e "_Potgood2" && -s "_Potgood2" ]];
			then 
				sed "s/\(^[ACGT][ACGT]*\)\(${myPrimer}..*\)/\2/" _Potgood2 > _good2
				#if _Potgood2 is not empty, we'll extract all those reads from the initial file (to have the original length)
				join -1 1 -2 2 -t $'\t' -o 1.2 1.1 <(sort -t $'\t' -k1,1 _not_found) <(sort -t $'\t' -k2,2 _good2) > _Potgood3
				#do it again, but keep the full sequence now 
				sed "s/\(^[ACGT][ACGT]*\)\(${myPrimer}..*\)/\2/" _Potgood3 > _Potgood4
				join -1 1 -2 2 -t $'\t' -o 1.1 1.2 1.3 1.4 2.1 <(sort -t $'\t' -k1,1 _not_found) <(sort -t $'\t' -k2,2 _Potgood4) >> _final_good_2
			fi
			#Finally we'll extract all the sequences where the primer has not been found
			cat _temp2 | grep -v "^${myPrimer}" | grep -v "${myPrimer}" > _rest
			join -1 1 -2 2 -t $'\t' -o 1.1 1.2 1.3 1.4 <(sort -t $'\t' -k1,1 _not_found) <(sort -t $'\t' -k2,2 _rest) > _temp3 |true
			mv _temp3 _not_found
		done<__primers
		
		#Let's use python to cut the sequences and quality (speed advantage over shell)
		cat <<EOF > _pythonScript.py
#!/usr/bin/python
import pandas as pd
from pandas import *
print "Pandas version : " + pandas.__version__
inputfile="_final_good_2"
data = pd.read_csv(inputfile, sep='\t', names=("id","seqO","plus","qualO","seqF"))

data["seqO_len"] = data.seqO.str.len()
data["seqF_len"] = data.seqF.str.len()
data["qualO_len"] = data.seqO.str.len()
data["removed"] = data["seqO_len"]-data["seqF_len"]
data = data.drop('seqO',1)
data.head(n=2)
#data['Report qualO'].str[data["removed"]:]

data2 = data.apply(lambda x: x['qualO'][x['removed']:], axis=1)
data3 = data.join(data2.to_frame(name="qualF"), how='outer')

data["qualF_len"] = data3.qualF.str.len()
data = data3[["id","seqF","plus","qualF",]]
data.to_csv("_final_good_2_pythoned", sep="\t", index=False, header=False)
EOF
		fastqName=$(basename ${myfastq} .fastq)
		initialReads=$(awk '{s++}END{print s/4}' ${myfastq})
		if [[ -e "_final_good_2" && -s "_final_good_2" ]];
		then 
			python _pythonScript.py
			cat _final_good_1 _final_good_2_pythoned | sed "s/\t/\n/g" > ${outputname}_withPrimers.fastq
			yesPrimers=$(awk '{s++}END{print s/4}'  ${outputname}_withPrimers.fastq)
			percentageKept=$(echo "scale=3; ${yesPrimers} / ${initialReads} * 100" | bc)
		else
			if [[ -e "_final_good_1" && -s "_final_good_1" ]]; then
				sed "s/\t/\n/g" _final_good_1 > ${outputname}_withPrimers.fastq
				yesPrimers=$(awk '{s++}END{print s/4}'  ${outputname}_withPrimers.fastq)
				percentageKept=$(echo "scale=3; ${yesPrimers} / ${initialReads} * 100" | bc)
			else
				touch ${outputname}_withPrimers.fastq
				yesPrimers=0
				percentageKept=0
			fi
		fi
		
		
		if [[ -e "_not_found" && -s "_not_found" ]]; then 
			cat _not_found | sed "s/\t/\n/g" > ${outputname}_noExactPrimerFound.fastq
			noPrimers=$(awk '{s++}END{print s/4}'  ${outputname}_noExactPrimerFound.fastq)
			mv ${outputname}_noExactPrimerFound.fastq ${dir0}/fastqsWithoutExactPrimers
		else
			touch ${dir0}/fastqsWithoutExactPrimers/${outputname}_noExactPrimerFound.fastq
			noPrimers=0
		fi

		rm -f _temp1 _temp2 _good1 _good2 _Potgood2 _rest _Potgood3 _Potgood4 _pythonScript.py _not_found _final_good_2_pythoned _final_good_1 __full_file __primers _final_good_2
		#stats
		echo -e "${fastqName}\t${initialReads}\t${yesPrimers}\t${noPrimers}\t${percentageKept}" >> primerSelection_summary.log
	fi
done<__readsList





cd ${dir0}/fastqsWithPrimers
echo -e "\n-----\nCleaning the Paired-end fastqs"
for i in *R1_withPrimers.fastq; 
do 
	core=$(basename $i R1_withPrimers.fastq)
	bash ${dirpipe}/fastq_PE_fixer.sh ${core}R1_withPrimers.fastq ${core}R2_withPrimers.fastq
done

cd ${dir0}/fastqsWithoutExactPrimers
echo -e "\n-----\nCleaning the Paired-end fastqs"
for i in *R1_noExactPrimerFound.fastq; 
do 
	core=$(basename $i R1_noExactPrimerFound.fastq)
	#echo -e "${i}     |    ${core}"
	bash ${dirpipe}/fastq_PE_fixer.sh ${core}R1_noExactPrimerFound.fastq ${core}R2_noExactPrimerFound.fastq
	mv ${core}R1_noExactPrimerFound_ordered.fastq ${core}R1.fastq
	mv ${core}R2_noExactPrimerFound_ordered.fastq ${core}R2.fastq
	rm -f ${core}R1_noExactPrimerFound.fastq ${core}R2_noExactPrimerFound.fastq
	gzip ${core}R1.fastq
	gzip ${core}R2.fastq
done


cd ${dir0}/fastqsWithPrimers
rm -f *_withPrimers.fastq

echo -e "Add statistics of processed files"
sed -i "1s/$/\tPEFixing\tPercentageReadsKept/" primerSelection_summary.log
i=2
while read readfile;
do
	fullname=$(basename "$readfile")
	filename="${fullname##*/}"
	outputname=$(basename "$filename" | cut -d. -f1)
	NumbSeqs=$(awk '{s++}END{print s/4}' ${outputname}_withPrimers_ordered.fastq)
	#now the percentage of reads kept
	initialSeqNumber=$(awk '{s++}END{print s/4}' ${dirreads}/${outputname}.fastq)
	percentageKept=$(echo "scale=3; ${NumbSeqs} / ${initialSeqNumber} * 100" | bc)
	mv ${outputname}_withPrimers_ordered.fastq ${outputname}.fastq
	sed -i "${i}s/$/\t${NumbSeqs}\t${percentageKept}/" primerSelection_summary.log
	unset NumbSeqs
	((i = i + 1))
done<__readsList



rm -f __readsList


#CLEAN WORKDIR

mv ${dir0}/fastqsWithPrimers ${dir0}/run/
mv ${dir0}/fastqsWithoutExactPrimers ${dir0}/run/
mkdir -p ${successDir}
mkdir -p ${dir0}/run/Summary
cp ${dir0}/run/fastqsWithPrimers/primerSelection_summary.log ${dir0}/run/Summary/primer_selection.txt
touch ${successDir}/primerSelection.ok
#remove fastq in raw reads directory if they were zipped
if [ "${removeFastq}"  == "YES" ]; then
	rm -f ${dirreads}/*.fastq
else
	cd ${dirreads}
	for i in *.fastq; do echo -e "zipping --> ${i}";gzip $i; done
fi


echo -e "\n---\nprimerSelection.sh exiting normally.\nWell done!"




rm -f __readsList