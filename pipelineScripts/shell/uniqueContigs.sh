#!/bin/bash
set -e

dir0="YourMainPath"
iniFile="${dir0}/metadata/pipe.ini"
dirpipe="${dir0}/pipelineScripts"
dirPipePython="${dirpipe}/python"
dirmothur="${dir0}/pipelineScripts/mothur"
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"

expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
procNumber=$(grep "^procNumber" ${iniFile} | cut -d"=" -f2)
#get the following 2 values from uniqueContigs folder (open the pdf and look how many amplicons you want to use)
amplicon_min_length=$(grep "^amplicon_min_length" ${iniFile} | cut -d"=" -f2)
amplicon_max_length=$(grep "^amplicon_max_length" ${iniFile} | cut -d"=" -f2)




echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "pipeline folder: $dirpipe"
echo -e "Mothur folder: $dirmothur"
echo -e "Python scripts folder: ${dirPipePython}"
echo -e "sample_type: ${sample_type}"
echo -e "Mothur program loction: $dirmothur"
echo -e "Number of processor used for Mothur: $procNumber"
echo -e "Minimum amplicon length: $amplicon_min_length"
echo -e "Maxmimum amplicon length: $amplicon_max_length"
echo -e "-----------------------------------\n\n"



mkdir -p ${dir0}/uniqueContigs
cd ${dir0}/uniqueContigs
ln -nsf ${dir0}/run/readsMerge/${expName}.trim.contigs.fasta
ln -nsf ${dir0}/run/readsMerge/${expName}.contigs.groups
ln -nsf ${dir0}/run/readsMerge/${expName}.trim.contigs.summary
echo -e "\n\nCreate Mothur file"
cat <<EOF > uniqueContigs.mothur
set.logfile(name=${expName}_uniqueContigs.log)
screen.seqs(fasta=${expName}.trim.contigs.fasta, group=${expName}.contigs.groups, summary=${expName}.trim.contigs.summary, minlength=${amplicon_min_length}, maxlength=${amplicon_max_length}, processors=${procNumber})
summary.seqs(fasta=${expName}.trim.contigs.good.fasta)
unique.seqs(fasta=${expName}.trim.contigs.good.fasta)
summary.seqs(fasta=${expName}.trim.contigs.good.unique.fasta, processors=${procNumber})
count.seqs(name=${expName}.trim.contigs.good.names, group=${expName}.contigs.good.groups)
summary.seqs(fasta=${expName}.trim.contigs.good.unique.fasta, count=${expName}.trim.contigs.good.count_table, processors=${procNumber})
EOF
echo -e "\n\nRun mothur"
${dirmothur}/mothur uniqueContigs.mothur

#Use a python script to extract length of contig sequences
python ${dirpipe}/python/fastaparser.py -i ${expName}.trim.contigs.good.unique.fasta -p ${expName}.trim.contigs.good.unique

echo -e "\n---\nParsing the name file: creating a mapping file between the uniq sequenece and all the sequences"
cd ${dir0}/uniqueContigs
bash ${dirpipe}/shell/mothur_uniq_report_parser.sh ${expName}.trim.contigs.good.names ${expName}_mapping_file.txt

#getting some number from mothur log:
grep "sequences\.$" ${expName}_uniqueContigs.log | rev | cut -d" " -f2 | rev | sort | uniq | sort -nr > _temp
originalSeq=$(head -n1 _temp)
filteredSeq=$(head -n2 _temp | tail -n1)
uniqueSeq=$(tail -n1 _temp)
rm -f _temp

#CLEAN WORKDIR
#moving folders
mv ${dir0}/uniqueContigs ${runDir}/
mkdir -p ${successDir}
mkdir -p ${runDir}/Summary
echo -e "Assembled Contigs\t${originalSeq}\nLength Filtered Contigs\t${filteredSeq}\nUnique Contigs\t${uniqueSeq}" > ${runDir}/Summary/uniqueContigs.txt
touch ${successDir}/uniqueContigs.ok

echo -e "\n---\nEnd of uniqueContigs.sh exiting normally\nWell done!"