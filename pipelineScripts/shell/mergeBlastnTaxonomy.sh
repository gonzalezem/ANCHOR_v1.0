#!/bin/bash
set -e

: <<'END'
Goal: Merge parsed BLASTn files with taxonomy parsed files.
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


echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Experiment name: ${expName}"
echo -e "Main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "Pipeline folder: $dirpipe"
echo -e "Database list: $databaseList"
echo -e "Identity threshold for BLASTn: $AnchorMinBlastIdentity"
echo -e "Identity and coverage thresholds for BLASTn: $lowCountSeqThreshold"
echo -e "-----------------------------------\n\n"




mkdir -p ${dir0}/mergeBlastnTaxonomy
allBDFile="${dir0}/mergeBlastnTaxonomy/blastnRes_anchors_vs_allDatabases.txt"
if [[ -e "${allBDFile}" && -s "${allBDFile}" ]];then 
	echo -e "${allBDFile} exists and not empty. Skipping this part."
else
	rm -f ${dir0}/mergeBlastnTaxonomy/allDatabase_parsed_tax.txt
	rm -f ${dir0}/mergeBlastnTaxonomy/blastnRes_anchors_vs_allDatabases.txt
	echo "Concatenate taxonomy and BLASTn output files"
	for db in ${databaseList}
	do
		cd ${dir0}/mergeBlastnTaxonomy
		#txonomy files
		cat ${runDir}/AnchorTaxonomyParser/${db}/${db}_curated_parsed.txt >> allDatabase_parsed_tax.txt
		#blastn output
		cat ${runDir}/parseAnchorBlastOutput/${db}/blastnRes_anchors_vs_${db}.txt >> blastnRes_anchors_vs_allDatabases.txt
	done
	sed -i '0,/^queryid/! {/^queryid/d}' blastnRes_anchors_vs_allDatabases.txt
fi


echo -e "Merge taxonomy and BLASTn output files"
if python ${dirpipe}/python/taxonomyMerger.py -i blastnRes_anchors_vs_allDatabases.txt -t allDatabase_parsed_tax.txt
then
	gzip < blastnTaxMerged.txt > blastnTaxMerged.txt.gz 
else
	echo -e "\n-----\ntaxonomyMerger.py exited with an error, check the log."
	exit 1
fi


#CLEAN WORKDIR
cd ${dir0}/mergeBlastnTaxonomy
mv ${dir0}/mergeBlastnTaxonomy ${dir0}/run/
mkdir -p ${successDir}
touch ${successDir}/mergeBlastnTaxonomy.ok



echo -e "\n---\nmergeBlastnTaxonomy.sh is exiting normally.\nWell done!"
