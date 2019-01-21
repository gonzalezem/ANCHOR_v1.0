#!/bin/bash
set -e

: <<'END'
Goal: Parse newlyt merged BLASTn + taxonomy files after icluding hits with no blast return (TrueUnknowns) and anchor counts (unique sequence count > anchor threshold)
END



dir0="YourMainPath"
iniFile="${dir0}/metadata/pipe.ini"
dirpipe="${dir0}/pipelineScripts"
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"
expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
databaseList=$(grep "^databaseList" ${iniFile} | cut -d"=" -f2)
AnchorMinBlastIdentity=$(grep "^AnchorMinBlastIdentity" ${iniFile} | cut -d"=" -f2)
lowCountSeqThreshold=$(grep "^lowCountSeqThreshold" ${iniFile} | cut -d"=" -f2)
primerSelectionBypass=$(grep "^primerSelectionBypass" ${iniFile} | cut -d"=" -f2)



echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Experiment name: ${expName}"
echo -e "Main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "Pipeline folder: $dirpipe"
echo -e "Database list: $databaseList"
echo -e "Primer selection: $primerSelectionBypass"
echo -e "Identity threshold for BLASTn: $AnchorMinBlastIdentity"
echo -e "Identity and coverage thresholds for BLASTn: $lowCountSeqThreshold"
echo -e "-----------------------------------\n\n"



##################################### CREATE A NEW SCRIPT FROM HERE ############################################




#Include hits with no blast return: TrueUnknowns
mkdir -p ${dir0}/anchorParser
cd ${dir0}/anchorParser
FILE="blastnTaxMerged_plusUnknowns.txt"
if [[ -e "${FILE}" && -s "${FILE}" ]];
then 
	echo -e "${FILE} exists and not empty. We'll skip that part!"
else
	ln -nsf ${runDir}/uniqueContigs/${expName}.trim.contigs.good.count_table
	ln -nsf ${runDir}/anchorSelection/anchors_seqList.txt
	ln -nsf ${runDir}/mergeBlastnTaxonomy/blastnTaxMerged.txt
	echo -e "Adding True unknown i.e. NoBlastHits to main data"
	join  -1 1 -2 1 -t $'\t' -v1 <(sort -t $'\t' -k1,1 anchors_seqList.txt) <(cut -f1 blastnTaxMerged.txt|sort|uniq) > _TrueUnknowns
	rm -f unknown_fake_data
	touch unknown_fake_data #in case there is no unknown
	while read unknown
	do
		unknownLength=$(grep "^${unknown}"$'\t' ${runDir}/uniqueContigs/${expName}.trim.contigs.good.unique_ParsedFasta.txt | cut -f2)
		echo -e "${unknown}\tTrueUnknown\t100\t${unknownLength}\t0\t0\t0\t0\t0\t0\t${unknownLength}\tTrue Unknown\t100\t0\t0\tNoBlastHit\tTrueUnknown\tTrueUnknown\tTrueUnknown;TrueUnknown;TrueUnknown;TrueUnknown;TrueUnknown;TrueUnknown;TrueUnknown\tTrueUnknown\tTrueUnknown\tTrueUnknown\tTrueUnknown\tTrueUnknown\tTrueUnknown\tTrueUnknown" >>unknown_fake_data
	done<_TrueUnknowns
	cat blastnTaxMerged.txt unknown_fake_data > blastnTaxMerged_plusUnknowns.txt
	rm -f _TrueUnknowns unknown_fake_data blastnTaxMerged.txt
	checkFile=$(grep $'\t'$'\t' blastnTaxMerged_plusUnknowns.txt | wc -l )
	if [ ${checkFile} -ne 0 ];then
		echo -e "There is a problem with:\n${dir0}/anchorParser/blastnTaxMerged_plusUnknowns.txt\n There are some double tabs which shouldn't happen"
		exit 1
	fi
fi

#Inlude Anchor counts: number of unique sequences above a given threshold
FILE="anchorsCount.txt"
if [[ -e "${FILE}" && -s "${FILE}" ]];
then 
	echo -e "${FILE} exists and not empty. Skipping multiplier counts calculations."
else
	ln -nsf ${runDir}/anchorSelection/anchors_seqList.txt
	cut -f1,2 ${expName}.trim.contigs.good.count_table | awk 'FNR>1'> _multiplier
	join  -1 1 -2 1 -t $'\t' -a 2 -o 2.1 1.2 -e ERROR <(sort -t $'\t' -k1,1 _multiplier) <(sort -t $'\t' -k1,1 anchors_seqList.txt) > _multiplier_anchors
	head -n1 ${expName}.trim.contigs.good.count_table | cut -f1,2 > _headers
	set +e
	grep "ERROR" _multiplier_anchors > ERROR_multiplier_anchors.txt
	if [[ -e "ERROR_multiplier_anchors.txt" && -s "ERROR_multiplier_anchors.txt" ]]; then
		echo -e "Check the errors in the following file:\n${dir0}/anchorParser/ERROR_multiplier_anchors.txt\nAll anchors (anchors_seqList.txt) should have a count within the file (${runDir}/uniqueContigs/${expName}.trim.contigs.good.count_table)"
		exit 1
	else
		rm -f ERROR_multiplier_anchors.txt
	fi
	set -e
	cat _headers _multiplier_anchors > anchorsCount.txt
	rm -f _headers _multiplier_anchors anchors_seqList.txt ${expName}.trim.contigs.good.count_table _multiplier
fi


#Parse all hits
FILE="deNovo_allDatabase_withTrueUnknowns_part4_parsed.txt"
if [[ -e "${FILE}" && -s "${FILE}" ]];
then 
	echo -e "${FILE} exists and not empty. Skipping python suite part4."
else
	echo -e "\nPython script suite part 4"
	ln -nsf ${runDir}/anchorSelection/anchors_fasta.tsv
	if [ "${primerSelectionBypass}"  == "YES" ]; then
		python ${dirpipe}/python/anchorParser.py -i blastnTaxMerged_plusUnknowns.txt -m anchorsCount.txt
	else
		ln -nsf ${dir0}/metadata/primers.txt
		python ${dirpipe}/python/anchorParser.py -i blastnTaxMerged_plusUnknowns.txt -m anchorsCount.txt -p
	rm -f _*
	fi
fi


rm -f ${dir0}/anchorParser/__lowCountersCountList primers.txt anchors_fasta.tsv

#CLEAN WORKDIR
cd ${dir0}/anchorParser
mv ${dir0}/anchorParser ${dir0}/run/
mkdir -p ${successDir}
touch ${successDir}/anchorParser.ok


echo -e "\n---\nanchorParser.sh is exiting normally.\nWell done!"


