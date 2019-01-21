#!/bin/bash
set -e

: <<'END'
What's needed: mothur software. Place it in ANCHOR/pipelineScripts/mothur/
What is it doing? Merge reads with Mothur and run a python script to parse and produce merge graphs
END

dir0="YourMainPath"
iniFile="${dir0}/metadata/pipe.ini"
dirpipe="${dir0}/pipelineScripts"
dirmothur="${dir0}/pipelineScripts/mothur"
dirphix=$(grep "^dirphix" ${iniFile} | cut -d"=" -f2)
expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"
dirreads="${runDir}/fastqsWithPrimers"
procNumber=$(grep "^procNumber" ${iniFile} | cut -d"=" -f2)


FILE="${dirmothur}/mothur"
if [ ! -e "$FILE" ];
then
   echo -e "You need a mothur file within ${dirmothur}/mothur folder. Download mothur and come back" >&2
   exit 1
fi

FILE="${dirpipe}/python/assembledContigsParser.py"
if [ ! -e "$FILE" ];
then
   echo -e "Python script $FILE does not exist. You need it to run. Find it and place it here: ${dirpipe}/python/" >&2
   exit 1
fi



echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "pipeline folder: $dirpipe"
echo -e "Mothur folder: $dirmothur"
echo -e "Reads directory: $dirreads"
echo -e "-----------------------------------\n\n"



cd ${dir0}

mkdir -p ${dir0}/readsMerge
cd ${dir0}/readsMerge
rm -f reads_filename_list.txt
ls -1 ${dirreads}/*.fastq | sort >> reads_filename_list.txt
sed "s/^..*\///g" reads_filename_list.txt  | sed "s/_R[12].fastq//g" | sort | uniq > samplelist.txt
echo -e "Create input file for mothur's make.contigs function"
rm -f ${expName}.files
while read sample
do
	set +e
	fwd=$(grep "/${sample}_R1.fastq$" reads_filename_list.txt)
	rvrs=$(grep "/${sample}_R2.fastq$" reads_filename_list.txt)
	set -e
	echo -e "${sample}    |    ${fwd}    |    ${rvrs}"
	if [[ -n "${fwd}" && -n "${rvrs}" ]]; then
		echo -e "${sample}\t${fwd}\t${rvrs}" >> ${expName}.files
	else
		echo -e "Couldn't find either\n${fwd}\nor\n${rvrs}\nfiles in ${dir0}/readsMerge folder check these and run this script again."
		exit 1
	fi
done<samplelist.txt


echo -e "\n\nCreate Mothur file"
echo -e "set.logfile(name=readsMerge.log)" > readsMerge.mothur
echo -e "make.contigs(file=${expName}.files, processors=${procNumber}, trimoverlap=F)" >> readsMerge.mothur
echo -e "summary.seqs(fasta=${expName}.trim.contigs.fasta, processors=${procNumber})" >> readsMerge.mothur

echo -e "\n\nRun mothur"
${dirmothur}/mothur readsMerge.mothur
#remove a tab in the fasta file (it is coming from a space within the fatsqs)
sed -i "s/\t$//" ${expName}.trim.contigs.fasta
cd ${dir0}/readsMerge
rm *.qual
cd ${dir0}/readsMerge
python ${dirpipe}/python/assembledContigsParser.py -r ${expName}.contigs.report

#CLEAN WORKDIR
#Zip reads 
cd ${dirreads}
for i in *.fastq; do echo -e "zipping --> ${i}";gzip $i; done
#moving folders
mv ${dir0}/readsMerge ${dir0}/run/
mkdir -p ${successDir}
mkdir -p ${dir0}/run/Summary
cp ${dir0}/run/readsMerge/${expName}_ReadsMerge_Average_Stats.txt ${dir0}/run/Summary/readsMerge.txt
cp ${dir0}/run/readsMerge/${expName}_Merge_Overview.pdf ${dir0}/run/Summary/
touch ${successDir}/readsMerge.ok

echo -e "---\nreadsMerge.pbs.sh exiting normally.\nWell done!"
