#!/bin/bash
set -e


: <<'END'
What you need:
1. USEARCH9 executable file in the folder /pipelineScripts/usearch9. Change the variable below if you prefer to use another path
what it does:
1. Flag any chimera with USEARCH9
END

dir0="YourMainPath"
iniFile="${dir0}/metadata/pipe.ini"
dirpipe="${dir0}/pipelineScripts"
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"
expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
cutoff=$(grep "^cutoff" ${iniFile} | cut -d"=" -f2)
AnchorMinBlastIdentity=$(grep "^AnchorMinBlastIdentity" ${iniFile} | cut -d"=" -f2)
lowCountSeqThreshold=$(grep "^lowCountSeqThreshold" ${iniFile} | cut -d"=" -f2)
usearch9="${dir0}/pipelineScripts/usearch9/usearch9"


echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Experiment name: ${expName}"
echo -e "Main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "Pipeline folder: $dirpipe"
echo -e "Database list: $databaseList"
echo -e "Anchor count threshold: $cutoff"
echo -e "Anchor identity threshold for BLASTn: $AnchorMinBlastIdentity"
echo -e "Identity and coverage thresholds for BLASTn: $lowCountSeqThreshold"
if [[ -e "${usearch9}" && -s "${usearch9}" ]];
then 
	echo -e "USEARCH version 9 path: ${usearch9}"
else 
	echo -e "\n-----\nError Detected:\nYou need USEARCH version 9 installed: download it here: https://www.drive5.com"
	exit 1
fi
echo -e "-----------------------------------\n\n"




echo -e "\n\n\n------------------------------------------------------------------\nCHIMERA DETECTION: UCHIME"
mkdir -p ${dir0}/chimeraFlag
cd ${dir0}/chimeraFlag
ln -nsf ${runDir}/anchorSelection/anchors_fasta.tsv
ln -nsf ${runDir}/uniqueContigs/${expName}.trim.contigs.good.count_table
#Add abundance to fasta file to comply with uchime requirements
join -1 1 -2 1 -t $'\t' -a 1 -o 1.1 2.2 1.2 -e EMPTY <(sort -t $'\t' -k1,1 anchors_fasta.tsv) <(cut -f-2 ${expName}.trim.contigs.good.count_table | sort -t $'\t' -k1,1) | sort -t $'\t' -k2,2 -nr > uchime_input
sed -i "s/^/>/" uchime_input
sed -i "s#\t#;size=#" uchime_input
sed -i "s#\t#;\n#" uchime_input



echo -e " Running UCHIME2"
${usearch9} -uchime2_denovo uchime_input -uchimeout results.uchime -abskew 16
checkPositives=$(grep "perfect_chimera" results.uchime | wc -l)
if [ ${checkPositives} -gt 0 ]; then
	#results are in the last column (Y or N) and the contig name in the second column. 
	grep "perfect_chimera" results.uchime | cut -d";" -f1 > positives.uchime
	sed -i "s/$/\tY/" positives.uchime
	grep -v "perfect_chimera" results.uchime | cut -d";" -f1 > negatives.uchime
	sed -i "s/$/\tN/" negatives.uchime
	cat positives.uchime negatives.uchime > parsedResults.uchime
	sed -i "1s/^/anchorID\tchimera\n/" parsedResults.uchime
	rm -f positives.uchime negatives.uchime
fi
rm -f uchime_input
#add column chimera to anchor table
ln -nsf ${runDir}/anchorParser/anchor_table.txt anchor_table_preChimera.txt
checkChimera="${dir0}/chimeraFlag/parsedResults.uchime"
if [[ -e "${checkChimera}" && -s "${checkChimera}" ]];
then 
	join -1 1 -2 1 -t $'\t' -a 1 <(awk 'FNR>1' anchor_table_preChimera.txt | sort -t $'\t' -k1,1) <(awk 'FNR>1' ${checkChimera} | sort -t $'\t' -k1,1) > _temp
	sed -i "s/\t$/\tN/" _temp
	cat <(head -n1 anchor_table_preChimera.txt | sed "1s/$/\tchimera/") _temp > anchor_table.txt
	rm -f _temp
else
	sed "s/$/\tN/" anchor_table_preChimera.txt > _temp
	sed -i "1s/\tN$/\tchimera/" _temp
	mv _temp anchor_table.txt
	rm -f _temp
fi


#CLEAN WORKDIR
cd ${dir0}/chimeraFlag
mv ${dir0}/chimeraFlag ${dir0}/run/
mkdir -p ${successDir}
touch ${successDir}/chimeraFlag.ok


echo -e "\n---\nchimeraFlag.sh is exiting normally.\nWell done!"




