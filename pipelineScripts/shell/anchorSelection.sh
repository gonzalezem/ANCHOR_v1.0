#!/bin/bash
set -e

: <<'END'
What does it need:
	1. NCBI BLASTn.
	2. databses index should in the db folder within the main directory (${dir0}/db). They should have databse_index name (ex: ${dir0}/db/nt_index)
What is it doing? blastn anchor sequences vs. Public repositories
END



dir0="YourMainPath"
iniFile="${dir0}/metadata/pipe.ini"
dirpipe="${dir0}/pipelineScripts"
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"
procNumber=$(grep "^procNumber" ${iniFile} | cut -d"=" -f2)
cutoff=$(grep "^cutoff" ${iniFile} | cut -d"=" -f2)
expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
databaseList=$(grep "^databaseList" ${iniFile} | cut -d"=" -f2)
AnchorMinBlastIdentity=$(grep "^AnchorMinBlastIdentity" ${iniFile} | cut -d"=" -f2)
wordSizeAnchors=$(grep "^wordSizeAnchors" ${iniFile} | cut -d"=" -f2)
bypass=$(grep "^bypass" $1 | cut -d"=" -f2)


cd ${dir0}

echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Experiment name: ${expName}"
echo -e "Main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "pipeline folder: $dirpipe"
echo -e "Anchor cutoff: $cutoff"
echo -e "databaseList: $databaseList"
echo -e "Identity threshold for BLASTn: $AnchorMinBlastIdentity"
echo -e "BLASTn WORDSIZE parameter value: $wordSizeAnchors"
echo -e "Number of processors used: $procNumber"
if [ -n "${bypass}" ]; then
	echo -e "Reference folder for blastn (bypass option): ${bypass}"
	if [ ! -d "$bypass" ]; then
  		echo -e "\n---\nERROR\nReference folder for bypassing blastn calculations doesn't exists.\nCheck: ${bypass}"
  		exit 1
	fi
else
	echo -e "bypass blastn calculations: FALSE"
fi

echo -e "-----------------------------------\n\n"



########## CHECK IF DATABASES INDEX EXIST WITHIN db folder ###################
set +e
checkDB=$(ls -1d -- ${dir0}/db/*_index/ | sed "s/\/$//" | wc -l)
set -e
if [ ${checkDB} -eq 0 ]
then
	echo -e "You need database index (ex: nt_index) in the db folder. For example:\n${dir0}/db/nt_index"
	echo -e "\nScript will stop now." >&2
	exit 1
fi


########## CHECK IF BLASTn EXISTS ###################
set +e

checkBLASTn=$(command -v blastn | wc -l)
if [ ${checkBLASTn} -eq 0 ]
then
	echo -e "You need BLASTn install in your path. Install it and come back"
	echo -e "\nScript will stop now." >&2
	exit 1
fi


mkdir -p ${dir0}/anchorSelection
cd ${dir0}/anchorSelection
ln -nsf ${runDir}/uniqueContigs/${expName}.trim.contigs.good.unique.fasta
ln -nsf ${runDir}/uniqueContigs/${expName}.trim.contigs.good.count_table
echo -e "Anchor selection: creation of anchor.fasta from chosen anchor cutoff (${cutoff})"
awk -F"\t" -v var="${cutoff}" '$2>=var' ${expName}.trim.contigs.good.count_table | awk 'FNR>1' | cut -f1 > anchors_seqList.txt
${dirpipe}/c/faSomeRecords ${expName}.trim.contigs.good.unique.fasta anchors_seqList.txt anchors.fasta
#create a tsv version
${dirpipe}/shell/fastaToTab.sh anchors.fasta anchors_fasta.tsv


for db in ${databaseList}
do
	HCZipFile="${dir0}/anchorSelection/BLASTn_${db}.zip"
	if [[ -e "${HCZipFile}" && -s "${HCZipFile}" ]];
	then 
		echo -e "${HCZipFile} exists and not empty. Skipping blastn HC vs databases"
	else
		if [ -n "${bypass}" ]; then
			echo -e "A blastn (anchor sequences vs ${db}) with a lower cutoff has already been run, we'll use its output"
			mkdir -p ${dir0}/anchorSelection/BLASTn_${db}/blastResults_${db}
			cd ${dir0}/anchorSelection/BLASTn_${db}/blastResults_${db}
			#get number of output from bypass folder
			find ${bypass}/run/anchorSelection/BLASTn_${db}/blastResults_${db} -type f -name "output_${db}_*.txt" -print > __outputFiles
			i=1
			while read output
			do
				join  -1 1 -2 1 -t $'\t' <(sort ${dir0}/anchorSelection/anchors_seqList.txt) <(sort -t $'\t' -k1,1 ${output}) > output_${db}_${i}.txt
				((i = i + 1))
			done<__outputFiles
			rm -f __outputFiles
		else
			echo -e "Creating blastn files for anchors: ${db}"
			rm -rf ${dir0}/anchorSelection/BLASTn_${db}
			cd ${dir0}/anchorSelection
			python ${dirpipe}/python/blastnPrep.py -i anchors.fasta -p ${dir0}/anchorSelection/BLASTn_${db} -cc 10000 -dp ${dir0}/db -dn ${db}
			#customize blastn
			sed -i "s/num_threads 12 -word_size 30 -perc_identity 80 -max_target_seqs 1/num_threads ${procNumber} -word_size ${wordSizeAnchors} -perc_identity ${AnchorMinBlastIdentity}/" ${dir0}/anchorSelection/BLASTn_${db}/scripts/blastn_vs_${db}.sh
			bash ${dir0}/anchorSelection/BLASTn_${db}/scripts/blastn_vs_${db}.sh
			cd ${dir0}/anchorSelection
			zip -r BLASTn_${db}.zip BLASTn_${db}
			rm -rf BLASTn_${db}/fastachunks
		fi
	fi
done


#CLEAN WORKDIR
cd ${dir0}/anchorSelection/
rm -f ${expName}.trim.contigs.good.unique.fasta ${expName}.trim.contigs.good.count_table 
mv ${dir0}/anchorSelection ${dir0}/run/
mkdir -p ${successDir}
touch ${successDir}/anchorSelection.ok


echo -e "\n---\nanchorSelection.sh exiting normally.\nWell done!"

