#!/bin/bash
set -e

: <<'END'
Goal: Parse the anchor sequences vs. public repositories blast output: extract hits with the percentage identity and the coverage (alignment length / contig length) above a given threshold.
END





dir0="YourMainPath"
iniFile="${dir0}/metadata/pipe.ini"
dirpipe="${dir0}/pipelineScripts"
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"
expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
databaseList=$(grep "^databaseList" ${iniFile} | cut -d"=" -f2)
AnchorMinBlastIdentity=$(grep "^AnchorMinBlastIdentity" ${iniFile} | cut -d"=" -f2)
primerSelectionBypass=$(grep "^primerSelectionBypass" $1 | cut -d"=" -f2)



echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Experiment name: ${expName}"
echo -e "main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "Pipeline folder: $dirpipe"
echo -e "Anchor identity and coverage threshold: $AnchorMinBlastIdentity"
echo -e "Database list: $databaseList"
if [ "${primerSelectionBypass}"  == "YES" ]; then
	blastnDir="anchorSelection"
else
	blastnDir="anchorSelection_primerIntegration"
fi
echo -e "Using blastn output from: ${dir0}/run/${blastnDir}"
echo -e "-----------------------------------\n\n"



cd ${dir0}


for db in ${databaseList}
do
	echo -e "\n---\nParsing ${db} BLAST output"
	mkdir -p ${dir0}/parseAnchorBlastOutput/${db}
	echo -e "Calculate coverage on query"
	cd ${dir0}/parseAnchorBlastOutput/${db}
	resFile="${dir0}/parseAnchorBlastOutput/${db}/${db}_speciesID_list.txt"
	if [[ -e "${resFile}" && -s "${resFile}" ]];then
		echo -e "${resFile} exists and not empty. Skipping blastn parsing"
	else
		rm -f blastnRes_anchors_vs_${db}.txt
		for output in ${dir0}/run/${blastnDir}/BLASTn_${db}/blastResults_${db}/output_*.txt
		do
			#Keep what will be needed from blast output
			cut -f-15 ${output} > _temp0
			#I'll remove the subject description (15th column) for now as some hits have so much inside that field that it reaches the limits of awk. I'll include  them later
			cut -f-14 _temp0 > _temp1
			#extract subjectid and subject description
			cut -f2,15 _temp0 | sort | uniq > _temp2
			#limit very very large annotation definition and add database information within the subject id
			cut -d, -f-2 _temp2 | sed "s/^/${db}|/" > _temp3
			#Calculate coverage on query
			awk -v OFS='\t' -v var="${db}" '{ print $1, var"|"$2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, "NoDescription", $4/$13*100, $14, $4/$14*100}' _temp1 > _temp4
			#I may have a coverage >100 when there are gaps within the alignment. In this case I'll lower the coverages to 100% to keep the gap penalty reflected upon the identity value only.
			awk -v OFS='\t' '{if ($15>100.0){$15=100.0} print}' _temp4 > _temp5
			awk -v OFS='\t' '{if ($17>100.0){$17=100.0} print}' _temp5 > _temp4
			#joining with the shortened subject descriptions
			join  -1 2 -2 1 -t $'\t' -a1 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 2.2 -e ERROR <(sort -t $'\t' -k2,2 _temp4) <(sort -t $'\t' -k1,1 _temp3) > _temp6 
			#For 16SMicrobial, change 16SMicrobial|gi|1277396189|ref|NR_151900.1| to 16SMicrobial|NR_151900.1
			if [ "${db}" == "16SMicrobial" ]; then
				paste <(cut -f1 _temp6) <(cut -f2 _temp6 | sed "s/|$//" | rev | cut -d"|" -f1 | rev | sed "s/^/16SMicrobial|/") <(cut -f3- _temp6) > _temp7
				mv _temp7 _temp6
			fi
			echo "keep only best hits (potential hits threshold is the rounded floor of the higest id & cov) to ${threshold}% identity and coverage hits"
			awk -F"\t" -v var="${AnchorMinBlastIdentity}" '($3>=var && $15>=var)' _temp6 > _temp7
			sed -i "1s/^/queryid\tsubjectid\tidentity\talignmentlength\tmismatches\tgapopens\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tquerylength\tquerydescription\tcoverageQ\tsubjectlength\tcoverageS\tsubjectdescription\n/" _temp7
			if [[ -e "_temp7" && -s "_temp7" ]];then
				if python ${dirpipe}/python/blastHitsParser.py -i _temp7
				then
					cat _output_py >> ${dir0}/parseAnchorBlastOutput/${db}/blastnRes_anchors_vs_${db}.txt
				else
					echo -e "blastHitsParser_for_part6_shell_scripts.py exited with an error, check the log"
					exit 1
				fi
			fi
			rm -f _*
		done

		#remove all duplicate headings
		awk -F"\t" '!seen[$0]++' blastnRes_anchors_vs_${db}.txt > _temp
		mv _temp blastnRes_anchors_vs_${db}.txt

		#check if there are some errors from the join above
		set +e
		checkErr=$(grep -c "ERROR"  blastnRes_anchors_vs_${db}.txt)
		set -e
		if [ ${checkErr} -gt 0 ]
		then
			grep "ERROR"  blastnRes_anchors_vs_${db}.txt > ERRORS.txt 
			echo "Not all the taxonomy have a definition which is not normal. Check ${dir0}/parseAnchorBlastOutput/${db}/ERRORS.txt"
			exit 1
		fi

		echo -e "Output ${db} hit list"
		if [ "${db}" == "nt" ]; then
			#remove the version of the hit in nt database (ex: KR708860.1 will become KR708860)
			cut -f2 blastnRes_anchors_vs_${db}.txt | awk 'FNR>1' | cut -d"." -f1 | sort | uniq > ${db}_speciesID_list.txt
		else
			cut -f2 blastnRes_anchors_vs_${db}.txt | awk 'FNR>1' | sort | uniq > ${db}_speciesID_list.txt
		fi
	fi
done

#CLEAN WORKDIR
cd ${dir0}/parseAnchorBlastOutput/
mv ${dir0}/parseAnchorBlastOutput ${dir0}/run/
mkdir -p ${successDir}
touch ${successDir}/parseAnchorBlastOutput.ok

echo -e "\n---\nparseAnchorBlastOutput exiting normally.\nWell done!"