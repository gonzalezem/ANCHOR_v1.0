#!/bin/bash
set -e

: <<'END'
What does it need:
	1. NCBI BLASTn.
What is it doing? blastn low-count sequences vs. anchor sequences
END



dir0="YourMainPath"
iniFile="${dir0}/metadata/pipe.ini"
dirpipe="${dir0}/pipelineScripts"
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"
procNumber=$(grep "^procNumber" ${iniFile} | cut -d"=" -f2)
cutoff=$(grep "^cutoff" ${iniFile} | cut -d"=" -f2)
expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
lowCountSeqThreshold=$(grep "^lowCountSeqThreshold" ${iniFile} | cut -d"=" -f2)
wordSizeLowCountSeq=$(grep "^wordSizeLowCountSeq" ${iniFile} | cut -d"=" -f2)


echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Experiment name: ${expName}"
echo -e "main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "pipeline folder: $dirpipe"
echo -e "Anchor cutoff: $cutoff"
echo -e "Identity and coverage thresholds for BLASTn: $lowCountSeqThreshold"
echo -e "BLASTn WORDSIZE parameter value: $wordSizeLowCountSeq"
echo -e "Number of processors used: $procNumber"
echo -e "-----------------------------------\n\n"





blastnZipFile="${dir0}/lowCountSeqSelection/BLASTn_anchors.zip"
if [[ -e "${blastnZipFile}" && -s "${blastnZipFile}" ]];
then 
	echo -e "Skipping blastn low count sequences vs. anchors!"
else 
	echo -e "Create an index for anchor sequences"
	mkdir -p ${dir0}/lowCountSeqSelection/anchors_index
	cd ${dir0}/lowCountSeqSelection/anchors_index
	ln -nsf ${runDir}/anchorSelection/anchors.fasta
	makeblastdb -in anchors.fasta -dbtype 'nucl' -out anchors -input_type fasta
	rm -f ${dir0}/lowCountSeqSelection/anchors.fasta

	echo -e "Extracting low count sequences from main fasta"
	cd ${dir0}/lowCountSeqSelection
	ln -nsf ${runDir}/uniqueContigs/${expName}.trim.contigs.good.unique.fasta
	ln -nsf ${runDir}/uniqueContigs/${expName}.trim.contigs.good.count_table
	awk -F"\t" -v var="${cutoff}" '$2<var' ${expName}.trim.contigs.good.count_table | awk 'FNR>1' | cut -f1 > lowCountSeq_seqList.txt
	${dirpipe}/c/faSomeRecords ${expName}.trim.contigs.good.unique.fasta lowCountSeq_seqList.txt lowCountSeq.fasta

	echo -e "Blasting low count sequences vs. anchor sequences"
	python ${dirpipe}/python/blastnPrep.py -i lowCountSeq.fasta -p ${dir0}/lowCountSeqSelection/BLASTn_anchors -cc 10000 -dp ${dir0}/lowCountSeqSelection -dn anchors
	#customize blastn
	sed -i "s/num_threads 12 -word_size 30 -perc_identity 80/num_threads ${procNumber} -word_size ${wordSizeLowCountSeq} -perc_identity ${lowCountSeqThreshold}/" ${dir0}/lowCountSeqSelection/BLASTn_anchors/scripts/blastn_vs_anchors.sh
	sed -i "s/qlen slen stitle qseq sseq/qlen slen/" ${dir0}/lowCountSeqSelection/BLASTn_anchors/scripts/blastn_vs_anchors.sh
	bash ${dir0}/lowCountSeqSelection/BLASTn_anchors/scripts/blastn_vs_anchors.sh
	cd ${dir0}/lowCountSeqSelection
	zip -r BLASTn_anchors.zip BLASTn_anchors
	rm -f ${expName}.trim.contigs.good.unique.fasta ${expName}.trim.contigs.good.count_table lowCountSeq.fasta
	rm -rf BLASTn_anchors/fastachunks
fi


echo -e "Parsing BLASTn results"
cd ${dir0}/lowCountSeqSelection
rm -f blastnRes_lowCountSeq_vs_anchors.txt
for output in ${dir0}/lowCountSeqSelection/BLASTn_anchors/blastResults_anchors/*.txt
do
	#add coverage (coverage on query(coverageQ) and on subject (coverageS)) columns
	awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14,  $4/$13*100, $4/$14*100 }' ${output} > _temp0

	#I may have a coverage >100 when there are gaps within the alignment. In this case I'll lower the coverages to 100% to keep the gap penalty reflected upon the identity value only.
	awk -v OFS='\t' '{if ($15>100.0){$15=100.0} print}' _temp0 > _temp1
	awk -v OFS='\t' '{if ($16>100.0){$16=100.0} print}' _temp1 > _temp0

	#use low count sequence identity and coverage thresholds to filter the output
	awk -F"\t" -v var="${lowCountSeqThreshold}" '($3>=var && $15>=var  && $16>=var)' _temp0 > _temp1
	#We'll select only the highest score (identity + coverageQ + coverageS)
	awk -v OFS='\t' '{ print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $3+$15+$16 }' _temp1 > _temp2
	sort -k1,1 -k17,17rg _temp2 > _temp3
	#remove duplicates if some are remaining (there is no way of choosing one before another now)
	awk -F"\t" '!seen[$1]++' _temp3 | cut -f-16 >> blastnRes_lowCountSeq_vs_anchors.txt
	rm -f _temp*
done
sed -i "1s/^/queryid\tsubjectid\tidentity\talignmentlength\tmismatches\tgapopens\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tquerylength\tsubjectlength\tAlignment_vs_Q_SeqLength\tAlignment_vs_S_SeqLength\n/" blastnRes_lowCountSeq_vs_anchors.txt


echo -e "Anchor count inflation due to low count sequences"
lastFILE="${dir0}/lowCountSeqSelection/anchor_countInflation_from_lowCountSeq.txt"
if [[ -e "${lastFILE}" && -s "${lastFILE}" ]];
then 
	echo -e "${lastFILE} exists and not empty. Skipping."
else 
	cd ${dir0}/lowCountSeqSelection
	cut -f2 blastnRes_lowCountSeq_vs_anchors.txt | awk 'FNR>1' | sort | uniq -c | sed "s/^  *//" | sed "s/ /\t/" | awk -F"\t" '{ print $2"\t"$1 }' > ${dir0}/lowCountSeqSelection/anchor_countInflation_from_lowCountSeq.txt
fi


#CLEAN WORKDIR
cd ${dir0}/lowCountSeqSelection/
mv ${dir0}/lowCountSeqSelection ${dir0}/run/
mkdir -p ${successDir}
touch ${successDir}/lowCountSeqSelection.ok

echo -e "\n---\nlowCountSeqSelection.sh exiting normally.\nWell done!"
