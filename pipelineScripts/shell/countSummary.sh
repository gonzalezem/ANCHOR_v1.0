#!/bin/bash
set -e


: <<'END'
what it does:
1. Create expected count from starting files (unique counts, blastn) in order to (i.e. later) confirm the counts obtained in the results.
END





dir0="YourMainPath"
iniFile="${dir0}/metadata/pipe.ini"
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"
expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
cutoff=$(grep "^cutoff" ${iniFile} | cut -d"=" -f2)
lowCountSeqThreshold=$(grep "^lowCountSeqThreshold" ${iniFile} | cut -d"=" -f2)


echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Experiment name: ${expName}"
echo -e "Main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "Database list: $databaseList"
echo -e "Anchor count threshold: $cutoff"
echo -e "Identity and coverage thresholds for BLASTn: $lowCountSeqThreshold"
echo -e "-----------------------------------\n\n"





cd ${dir0}

echo -e "-----\nComputing counts, might take some time..."
mkdir -p ${dir0}/countSummary
cd ${dir0}/countSummary
ln -nsf ${runDir}/uniqueContigs/${expName}.trim.contigs.good.count_table __raw_counts
#1
totalAmplicons=$(cut -f2 __raw_counts | awk 'FNR>1' | awk ' { sum+=$1 } END { print sum }')
echo -e "totalAmplicons\t${totalAmplicons}" > statistics.txt
uniqueSequences=$(cut -f2 __raw_counts | awk 'FNR>1' | wc -l)
echo -e "uniqueSequences\t${uniqueSequences}" >> statistics.txt
#2
highCounterSeqs=$(cut -f2 __raw_counts | awk 'FNR>1' | awk -F"\t" -v var="${cutoff}" '$1>=var' | wc -l)
echo -e "highCounterSeqs\t${highCounterSeqs}" >> statistics.txt
anchorsCumCount=$(cut -f2 __raw_counts | awk 'FNR>1' | awk -F"\t" -v var="${cutoff}" '$1>=var' | awk ' { sum+=$1 } END { print sum }')
echo -e "anchorsCumCount\t${anchorsCumCount}" >> statistics.txt
#3
LowCountSeqeqs=$(cut -f2 __raw_counts | awk 'FNR>1' | awk -F"\t" -v var="${cutoff}" '$1<var' | wc -l)
echo -e "LowCountSeqeqs\t${LowCountSeqeqs}" >> statistics.txt
LowCountSeqCumCount=$(cut -f2 __raw_counts | awk 'FNR>1' | awk -F"\t" -v var="${cutoff}" '$1<var' | awk ' { sum+=$1 } END { print sum }')
echo -e "LowCountSeqCumCount\t${LowCountSeqCumCount}" >> statistics.txt
#3 (check)
singleCountLowCountSeq=$(cut -f2 __raw_counts | awk 'FNR>1' | awk '$1<2' | wc -l)
multipleCountLowCountSeq=$(cut -f2 __raw_counts | awk 'FNR>1' | awk -F"\t" -v var="${cutoff}" '($1<var && $1>1)' | wc -l)
singleCountLowCountSeqCumCount=$(cut -f2 __raw_counts | awk 'FNR>1' | awk '$1<2' | awk ' { sum+=$1 } END { print sum }')
multipleCountLowCountSeqCumCoun=$(cut -f2 __raw_counts | awk 'FNR>1' | awk -F"\t" -v var="${cutoff}" '($1<var && $1>1)' | awk ' { sum+=$1 } END { print sum }')
echo -e "multipleCountLowCountSeq\t${multipleCountLowCountSeq}" >> statistics.txt
echo -e "multipleCountLowCountSeqCumCoun\t${multipleCountLowCountSeqCumCount}" >> statistics.txt
cut -f-2 __raw_counts | awk 'FNR>1' | awk -F'\t' '$2<2' > __singleCountLowCountSeq_List_and_counts
cut -f-2 __raw_counts | awk 'FNR>1' | awk -F"\t" -v var="${cutoff}" '($2<var && $2>1)' > __multipleCountLowCountSeq_List_and_counts


cd ${dir0}/countSummary
ln -nsf ${runDir}/uniqueContigs/${expName}_mapping_file.txt __${expName}_mapping_file
cut -f-2 __raw_counts | awk 'FNR>1' | awk -F"\t" -v var="${cutoff}" '$2>=var' | cut -f1 | sort | uniq > __anchorsList
join -t $'\t' -1 1 -2 1 <(sort __anchorsList) <(sort -t $'\t' -k1,1 __${expName}_mapping_file) > __anchors_mapping_file
cut -f-2 __raw_counts | awk 'FNR>1' | awk -F"\t" -v var="${cutoff}" '$2<var' | cut -f1 | sort | uniq > __LowCountSeqList
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k1,1 __LowCountSeqList) <(sort -t $'\t' -k1,1 __${expName}_mapping_file) > __LowCountSeq_mapping_file
	
	
	
#EXPECTED FROM BLAST RESULTS:
cd ${dir0}/countSummary
ln -nsf ${runDir}/lowCountSeqSelection/blastnRes_lowCountSeq_vs_anchors.txt __lowCountSeq_vs_anchors
cut -f-2 __lowCountSeq_vs_anchors | awk 'FNR>1' | sort > __lowCountSeq_vs_anchors_BlastMappingFile
cat __anchors_mapping_file <(awk -F"\t" '{ print $2"\t"$1 }' __lowCountSeq_vs_anchors_BlastMappingFile) >__anchors_BlastMapping_file_added_LowCountSeq #2 cols:  anchors    LowCountSeq
#etarct the multiple LowCountSeq from ${lowCountSeqThreshold} blast
join -t $'\t' -1 1 -2 1 <(cut -f1 __lowCountSeq_vs_anchors_BlastMappingFile | sort) <(cut -f1 __multipleCountLowCountSeq_List_and_counts | sort) > __multipleCountLC_in_blast
join -t $'\t' -1 1 -2 1 <(cut -f1 __lowCountSeq_vs_anchors_BlastMappingFile | sort) <(cut -f1 __singleCountLowCountSeq_List_and_counts | sort) > __singleCountLC_in_blast
multipleCountLC_in_blast_CumCount=$(join -t $'\t' -1 1 -2 1 <(sort __multipleCountLC_in_blast) <(sort -t $'\t' -k1,1 __multipleCountLowCountSeq_List_and_counts) | cut -f2 | paste -sd+ - | bc)
echo -e "multipleCountLC_in_blast_CumCount ${lowCountSeqThreshold}%\t${multipleCountLC_in_blast_CumCount}" >> statistics.txt
singleCountLC_in_blast_CumCount=$(wc -l __singleCountLC_in_blast | cut -d" " -f1)
echo -e "singleCountLC_in_blast_CumCount ${lowCountSeqThreshold}%\t${singleCountLC_in_blast_CumCount}" >> statistics.txt
totalExpectedCountsFromBlastOutput=$((${singleCountLC_in_blast_CumCount} + ${multipleCountLC_in_blast_CumCount} + ${anchorsCumCount} ))
echo -e "totalExpectedCountsFromBlastOutput ${lowCountSeqThreshold}%\t${totalExpectedCountsFromBlastOutput}" >> statistics.txt

#creating a universal mapping file with columns: HC    multipleCountLC    singleCountLC
#a few steps:
#1. create a maping file: multipleCountLC    100% identical sequence name
join -t $'\t' -1 1 -2 1 <(sort -t $'\t' -k1,1 __multipleCountLC_in_blast) <(sort -t $'\t' -k1,1 __LowCountSeq_mapping_file) > __multipleCountLC_in_blast_mapping_file #2 cols: multipleCountLC    100% mapping sequence name
#2. creating a 3 column mapping file: HC   multipleCountLC     100% identical sequence name
join -t $'\t' -1 1 -2 2 -o 2.1 2.2 1.2 <(sort -t $'\t' -k1,1 __multipleCountLC_in_blast_mapping_file) <(sort -t $'\t' -k2,2 __anchors_BlastMapping_file_added_LowCountSeq) > __HC_vs_multipleCountLC_in_blast_mapping_file
#3. creating a 3 column mapping file: HC   -     singleCounterLC
join -t $'\t' -1 2 -2 1 -o 1.1 2.1 <(sort -t $'\t' -k2,2 __anchors_BlastMapping_file_added_LowCountSeq) <(sort __singleCountLC_in_blast) | sed "s/\t/\t-\t/"> __anchors_BlastMapping_file_added_LowCountSeq_NoMultipleCountLC
#4. Creating a 3 column mapping file: HC    -    100% identical sequence name
sed "s/\t/\t-\t/" __anchors_mapping_file > __anchors_modified_mapping_file
#5. Finally geting what we wanted
cat __anchors_modified_mapping_file __HC_vs_multipleCountLC_in_blast_mapping_file __anchors_BlastMapping_file_added_LowCountSeq_NoMultipleCountLC | sort -t $'\t' -k1,1 | sed "1s/^/anchors\tmultipleCountLowCountSeq\tsingleCountLowCountSeq\n/" > afterBlast_mappingFile.txt

rm -f ${dir0}/countSummary/__*


#CLEAN WORKDIR
cd ${dir0}/countSummary
mkdir -p ${runDir}/Summary
cp ${dir0}/countSummary/statistics.txt ${runDir}/Summary/counts_overview.txt
mv ${dir0}/countSummary ${dir0}/run/
mkdir -p ${successDir}
touch ${successDir}/countSummary.ok


echo -e "\n---\ncountSummary.sh is exiting normally.\nWell done!"

