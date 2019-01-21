#!/bin/bash
set -e


: <<'END'
What you need:
1. anchor_table.txt from OTUandCounts.sh

what it does:
1. creates the OTU count + taxonomy file


END




dir0="YourMainPath"
iniFile="${dir0}/metadata/pipe.ini"
dirpipe="${dir0}/pipelineScripts"
runDir="${dir0}/run"
dirMeta="${dir0}/metadata"
successDir="${runDir}/successfulRuns"
expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
cutoff=$(grep "^cutoff" ${iniFile} | cut -d"=" -f2)
AnchorMinBlastIdentity=$(grep "^AnchorMinBlastIdentity" ${iniFile} | cut -d"=" -f2)
lowCountSeqThreshold=$(grep "^lowCountSeqThreshold" ${iniFile} | cut -d"=" -f2)
primerSelectionBypass=$(grep "^primerSelectionBypass" ${iniFile} | cut -d"=" -f2)


echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Experiment name: ${expName}"
echo -e "Main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "Pipeline folder: $dirpipe"
echo -e "Database list: $databaseList"
echo -e "Design file: $Conditions"
echo -e "Anchor count threshold: $cutoff"
echo -e "Anchor identity threshold for BLASTn: $AnchorMinBlastIdentity"
echo -e "Primer selection: $primerSelectionBypass"
echo -e "Identity and coverage thresholds for BLASTn: $lowCountSeqThreshold"
echo -e "-----------------------------------\n\n"





cd ${dir0}

mkdir -p ${dir0}/OTUandCounts
cd ${dir0}/OTUandCounts

echo -e "Calculating counts per samples"
ln -nsf ${runDir}/chimeraFlag/anchor_table.txt
ln -nsf ${runDir}/countSummary/afterBlast_mappingFile.txt
if [ ! "${primerSelectionBypass}"  == "YES" ]; then
	numberAmbiguousNucl=$(awk 'FNR>1' ${dirMeta}/primers.txt | grep -e "R" -e "Y" -e "S" -e "W" -e "K" -e "M" -e "B" -e "V" -e "D" -e "H" -e "N" | wc -l)
	echo -e "\n- Primers have ${numberAmbiguousNucl} ambiguous nucleotides"
fi

if [ "${primerSelectionBypass}"  == "YES" ]; then
	cp afterBlast_mappingFile.txt contig_crib.txt
	rm -f afterBlast_mappingFile.txt
elif [ ${numberAmbiguousNucl} -eq 0 ]; then
	cp afterBlast_mappingFile.txt contig_crib.txt
	rm -f afterBlast_mappingFile.txt
else
	#Now I need to replace the anchors that were demoted (as a result to ambiguous nuclotides in the primer sequence, see part4_bitscoreParser.py in part8_merging_all_database.pbs.sh) by a reference anchor.
	#In part8_merging_all_database.pbs.sh, we generated a crib for them(2 columns: referenceAnchor	anchorDemotedByPrimerAmbiguity). I will just replace the demoted by the reference in afterBlast_mappingFile.txt
	ln -nsf ${runDir}/anchorParser/crib_primerAmbiguity.txt
	originalAnchorsCumCount=$(awk 'FNR>1' afterBlast_mappingFile.txt | wc -l | cut -d" " -f1)
	originalAnchorNumber=$(cut -f1 afterBlast_mappingFile.txt | awk 'FNR>1' | sort | uniq | wc -l)
	join -1 2 -2 1 -t $'\t' -o 1.1 2.2 2.3 <(sort -t $'\t' -k2,2 crib_primerAmbiguity.txt) <(sort -t $'\t' -k1,1 afterBlast_mappingFile.txt) > _temp
	cat <(head -n1 afterBlast_mappingFile.txt) _temp > contig_crib.txt
	rm -f _temp afterBlast_mappingFile.txt
	trueAnchorsCumCount=$(awk 'FNR>1' contig_crib.txt | wc -l | cut -d" " -f1)
	trueAnchorsNumber=$(cut -f1 anchor_table.txt | awk 'FNR>1' | sort | uniq | wc -l)
	echo -e "Original number of Anchors\t${originalAnchorNumber}\noriginal Anchors Cumulative Count\t${originalAnchorsCumCount}\nPrimer ambiguity filtered anchors\t${trueAnchorsNumber}\nPrimer ambiguity filtered anchors total count\t${trueAnchorsCumCount}" > ${runDir}/Summary/anchor_collapsing_due_to_primers.txt
fi
#running a python script to extract final counts
python ${dirpipe}/python/OTUandCounts.py -i anchor_table.txt -c contig_crib.txt

rm -f ConditionList crib_primerAmbiguity.txt _temp anchor_table.txt



#CLEAN WORKDIR
cd ${dir0}/OTUandCounts
mv ${dir0}/OTUandCounts ${dir0}/run/
mkdir -p ${successDir}
touch ${successDir}/OTUandCounts.ok


echo -e "\n---\nOTUandCounts.sh is exiting normally.\nWell done!"


