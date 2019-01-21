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



##################################### CREATE A NEW SCRIPT FROM HERE ############################################



#creating multiplier file
echo -e "Inflate unique counts with the counts coming from low counters"
set +e
dirLowCounters="${dir0}/part5_HC_blastn_LC_binning/low_counters"
set -e
if [ ! -d "${dirLowCounters}" ]; then
	echo -e "\n---\n${dirLowCounters} doesn't exist although it should. Did you run the part5 bash scripts?"
	exit 1
fi
set +e
#ls -1 ${dirLowCounters}/low_counters_inflation_at*.txt > __lowCountersCountList
rm -f __lowCountersCountList
for idLC in ${lowCountSeqThreshold}
do
	echo -e "${dirLowCounters}/low_counters_inflation_at_${idLC}cutoff.txt" >> __lowCountersCountList
done
set -e
if [[ -e "__lowCountersCountList" && -s "__lowCountersCountList" ]];
then
    while read lowCounters
    do
    	shorterName=$(basename ${lowCounters} .txt)
    	final=$(echo -e "${shorterName}" | rev | cut -d"/" -f1 | rev)
    	mkdir -p ${dir0}/mergeBlastnResults/${final}
    	cd ${dir0}/mergeBlastnResults/${final}
    	FILE="allDatabase_withTrueUnknowns_part3_parsed.txt"
		if [[ -e "${FILE}" && -s "${FILE}" ]];
		then 
			echo -e "${FILE} exists and not empty. We'll skip that part!"
		else
			ln -nsf ${dir0}/part4_unique_sequences/${expName}.trim.contigs.good.count_table
			ln -nsf ${dir0}/part5_HC_blastn_LC_binning/high_counters/high_counters_seqList.txt
  			ln -nsf ${lowCounters} ${final}.txt
  			ln -nsf ${dir0}/mergeBlastnResults/blastnTaxMerged.txt
			echo -e "Adding True unknown i.e. NoBlastHits to main data"
			join  -1 1 -2 1 -t $'\t' -v1 <(sort -t $'\t' -k1,1 high_counters_seqList.txt) <(cut -f1 blastnTaxMerged.txt|sort|uniq) > _TrueUnknowns
			i=1
			rm -f unknown_fake_data
			touch unknown_fake_data #in case there is no unknown
			while read unknown
			do
				unknownLength=$(grep "^${unknown}"$'\t' ${dir0}/part4_unique_sequences/${expName}.trim.contigs.good.unique_ParsedFasta.txt | cut -f2)
				echo -e "${unknown}\tTrueUnknown_${i}\t100\t${unknownLength}\t0\t0\t0\t0\t0\t0\t${unknownLength}\tTrue Unknown\t100\t0\t0\tNoBlastHit\tTrueUnknown_${i}\tTrueUnknown_${i}\tTrueUnknown;TrueUnknown${i};TrueUnknown${i};TrueUnknown${i};TrueUnknown${i};TrueUnknown${i};TrueUnknown${i}\tTrueUnknown\tTrueUnknown_${i}\tTrueUnknown_${i}\tTrueUnknown_${i}\tTrueUnknown_${i}\tTrueUnknown_${i}\tTrueUnknown_${i}" >>unknown_fake_data
				((i = i + 1))
			done<_TrueUnknowns
			cat blastnTaxMerged.txt unknown_fake_data > allDatabase_withTrueUnknowns_part3_parsed.txt
			checkFile=$(grep $'\t'$'\t' allDatabase_withTrueUnknowns_part3_parsed.txt | wc -l )
			if [ ${checkFile} -ne 0 ];then
				echo -e "There is a problem with:\n${dir0}/mergeBlastnResults/${final}/allDatabase_withTrueUnknowns_part3_parsed.txt\n There are some double tabs which shouldn't happen"
				exit 1
			fi
		fi


		FILE="multiplier_plus_${final}.txt"
		if [[ -e "${FILE}" && -s "${FILE}" ]];
		then 
			echo -e "${FILE} exists and not empty. Skipping multiplier counts calculations."
		else
			ln -nsf ${dir0}/part5_HC_blastn_LC_binning/high_counters/high_counters_seqList.txt
			cut -f1,2 ${expName}.trim.contigs.good.count_table | awk 'FNR>1'> _multiplier
			join  -1 1 -2 1 -t $'\t' -a 2 -o 2.1 1.2 -e ERROR <(sort -t $'\t' -k1,1 _multiplier) <(sort -t $'\t' -k1,1 high_counters_seqList.txt) > _multiplier_highCounters
			head -n1 ${expName}.trim.contigs.good.count_table | cut -f1,2 > _headers
			set +e
			grep "ERROR" _multiplier_highCounters > ERROR_multiplier_highCounters.txt
			if [[ -e "ERROR_multiplier_highCounters.txt" && -s "ERROR_multiplier_highCounters.txt" ]]; then
				echo -e "Check the errors in the following file:\n${dir0}/mergeBlastnResults/${final}/ERROR_multiplier_highCounters.txt\nAll high counters (high_counters_seqList.txt) should have a count within the file (${dir0}/part4_unique_sequences/${expName}.trim.contigs.good.count_table)"
				exit 1
			else
				rm -f ERROR_multiplier_highCounters.txt unknown_fake_data
			fi
			set -e
			cat _headers _multiplier_highCounters > anchorsCount.txt
			rm -f _headers _multiplier_highCounters
		fi
		FILE="deNovo_allDatabase_withTrueUnknowns_part4_parsed.txt"
		if [[ -e "${FILE}" && -s "${FILE}" ]];
		then 
			echo -e "${FILE} exists and not empty. Skipping python suite part4."
		else
			echo -e "\nPython script suite part 4"
			check_file="${dir0}/part5_HC_blastn_LC_binning/high_counters/high_counters_fasta.tsv"
			if ! [[ -e "${FILE}" && -s "${FILE}" ]];
			then 
				cd ${dir0}/part5_HC_blastn_LC_binning/high_counters/
				${dirpipe}/shell/fastaToTab.sh high_counters.fasta high_counters_fasta.tsv
				sed -i "s/\t\t/\t/" high_counters_fasta.tsv
				cd ${dir0}/mergeBlastnResults/${final}
			fi
			ln -nsf ${dir0}/part5_HC_blastn_LC_binning/high_counters/high_counters_fasta.tsv
			if [ "${part2Bypass}"  == "YES" ]; then
				python ${dirpipe}/python/part4_bitscoreParser.py -i allDatabase_withTrueUnknowns_part3_parsed.txt -m anchorsCount.txt
			else
				ln -nsf ${dir0}/Metadata/adapters_and_primers.txt
				python ${dirpipe}/python/part4_bitscoreParser.py -i allDatabase_withTrueUnknowns_part3_parsed.txt -m anchorsCount.txt -p
			rm -f _*
			fi
		fi
    done<__lowCountersCountList
else 
    echo -e "${dirLowCounters}/low_counters_inflation_at*.txt does not exist or is empty. Did the part5 bash script ran without any error? Check and come back"
    exit 1 
fi

rm -f ${dir0}/mergeBlastnResults/__lowCountersCountList
echo -e "\n---\npart8_merging_all_database.pbs.sh is exiting normally.\nWell done dude!"


