#!/bin/bash
set -e


: <<'END'
What you need:
1. anchor_table.txt from chimeraFlag.sh

what it does:
1. creates the OTU count + taxonomy file


END

dir0="/media/emmanuel/storage2/16S_desrosiers_AZI"
iniFile="${dir0}/metadata/pipe.ini"
dirpipe="${dir0}/pipelineScripts"
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"
expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
cutoff=$(grep "^cutoff" ${iniFile} | cut -d"=" -f2)
AnchorMinBlastIdentity=$(grep "^AnchorMinBlastIdentity" ${iniFile} | cut -d"=" -f2)
lowCountSeqThreshold=$(grep "^lowCountSeqThreshold" ${iniFile} | cut -d"=" -f2)
primerSelectionBypass=$(grep "^primerSelectionBypass" ${iniFile} | cut -d"=" -f2)
databaseList=$(grep "^databaseList" ${iniFile} | cut -d"=" -f2)
design=${dir0}/metadata/design.txt



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
echo -e "Design file: ${design}"
echo -e "-----------------------------------\n\n"



if [[ -e "${design}" && -s "${design}" ]];
then 
	conditions=$(head -n1 ${design} | cut -f2- | sed "s/\t/\n/g" | grep -v -e "Samples" -e "Replicates" | sed "s/\n/ /")
else 
	echo -e "${design} does not exist or is empty. Rerun with a design file (${design})"
	exit 1
fi



mkdir -p ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/Excel
cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}
cp ${runDir}/OTUandCounts/OTU_table.txt ./
cp -r ${runDir}/Summary/ ./
cp ${runDir}/chimeraFlag/anchor_table.txt ./

if [ ! "${primerSelectionBypass}"  == "YES" ]; then
	echo -e  "Primer selection was selected."
	echo -e "Adding subject sequences column (db_sequences) to OTU and Anchor tables"
	#create a list of hits:
	RefAnchorIDColNumber=$(head -n1 OTU_table.txt | sed "s/\t/\n/g" | grep -n "RefAnchorID" | cut -d":" -f1)
	subjectidColNumber=$(head -n1 OTU_table.txt | sed "s/\t/\n/g" | grep -n "subjectid$" | cut -d":" -f1)
	cut -f${RefAnchorIDColNumber},${subjectidColNumber} OTU_table.txt | awk 'FNR>1' | cut -d";" -f1 > _temp1
	rm -f bestBlastnHits_per_Anchor.txt
	for db in ${databaseList}
	do
		blastDir="${runDir}/anchorSelection_primerIntegration/BLASTn_${db}/blastResults_${db}"
		anyhits=$(grep "${db}|" _temp1 | wc -l )
		if [ ${anyhits} -gt 0 ];then
			grep "${db}|" _temp1 | sed  "s/\t${db}|/\t/" > _temp2
			while read refAnchor subjectID
			do
				for output in ${blastDir}/output_*.txt
				do
					anyhits2=$(grep "${subjectID}" ${output} | grep "^${refAnchor}"$'\t' | wc -l)
					if [ ${anyhits2} -gt 0 ];then
						grep "${subjectID}" ${output} | grep "^${refAnchor}"$'\t' | sort | uniq | sort -t $'\t' -k12,12nr | head -n1 | cut -f1,2,12,17 | sed "s/$/\t${db}/" >> bestBlastnHits_per_Anchor.txt
					fi
				done
			done<_temp2
		fi
	done
	#Now I have one hit per representative Anchor, I'll join the 2 tables
	#in case I have multiple anchors in different runs, I'll select for the best bitscore
	sort -t $'\t' -k3,3nr bestBlastnHits_per_Anchor.txt | awk -F"\t" '!seen[$1,$3,$5]++' | cut -f-2,4- > _tmp1
	cp _tmp1 bestBlastnHits_per_Anchor.txt
	for db in ${databaseList}
	do
		anyhits=$(grep "${db}$" _tmp1 | wc -l )
		if [ ${anyhits} -gt 0 ];then
			grep "${db}$" _tmp1 > _tmp2
			while read anchorID subjectID fakeSeq db
			do
				crib="${runDir}/anchorSelection_primerIntegration/${db}/crib_original_seq_primer_transf_seq.txt"
				anyhits2=$(grep "${fakeSeq}$" ${crib} | wc -l)
				if [ ${anyhits2} -gt 0 ];then
					original=$(grep "${fakeSeq}$" ${crib} | head -n1 | cut -f2)
					sed -i "s/${fakeSeq}/${original}/" bestBlastnHits_per_Anchor.txt
				fi
			done<_tmp2
		fi
	done
	rm -f _tmp1 _tmp2
	sed -i "1s/^/RefAnchorID\tsubjectID\tdb_seq\n/" bestBlastnHits_per_Anchor.txt
	sequenceColNumber=$(head -n1 OTU_table.txt | sed "s/\t/\n/g" | grep -n "sequence" | cut -d":" -f1)
	#cut -f-${sequenceColNumber} OTU_table.txt > _temp3
	join -1 ${RefAnchorIDColNumber} -2 1 -t $'\t' -a 1 -a 2 <(sort -t $'\t' -k${RefAnchorIDColNumber},${RefAnchorIDColNumber} OTU_table.txt) <(cut -f1,3 bestBlastnHits_per_Anchor.txt | sort -t $'\t' -k1,1) > _temp4
	rest=$((${RefAnchorIDColNumber} + 1))
	paste <(cut -f2-${RefAnchorIDColNumber} _temp4) <(cut -f1 _temp4) <(cut -f${rest}- _temp4) | sort | uniq > _temp5
	grep "^OTU" _temp5 > __header
	grep -v "^OTU" _temp5 > _body
	cat __header _body > _temp6

	sequenceColNumber=$(head -n1 _temp6 | sed "s/\t/\n/g" | grep -n "sequence" | cut -d":" -f1)
	db_seqColNumber=$(head -n1 _temp6 | sed "s/\t/\n/g" | grep -n "db_seq" | cut -d":" -f1)
	startCountCols=$((${sequenceColNumber} + 1))
	stopCol=$((${db_seqColNumber} - 1))
	paste <(cut -f-${sequenceColNumber} _temp6) <(cut -f${db_seqColNumber} _temp6) <(cut -f${startCountCols}-${stopCol} _temp6) > _temp7
	#remove duplicates (that can happen when there are ambiguous bases within the repository sequence)
	awk -F"\t" '!seen[$1,$4,$5,$6,$7,$9,$10,$22,$112]++' _temp7 > _temp8
	mv _temp8 OTU_table.txt
	rm -f _* bestBlastnHits_per_Anchor.txt
fi

#sort OTU table per sample name
cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}
if [ ! "${primerSelectionBypass}"  == "YES" ]; then
	sequenceCol=$(head -n1 OTU_table.txt | sed "s/\t/\n/g" | grep -n "db_seq" | cut -d":" -f1)
else
	sequenceCol=$(head -n1 OTU_table.txt | sed "s/\t/\n/g" | grep -n "sequence" | cut -d":" -f1)
fi
totalcountsCol=$(head -n1 OTU_table.txt | sed "s/\t/\n/g" | grep -n "totalcounts" | cut -d":" -f1)
startCountCols=$((${sequenceCol} + 1))
endCountCols=$((${totalcountsCol} - 1))
cut -f${startCountCols}-${endCountCols} OTU_table.txt > _counts1
cat <<EOF > transposer.py
#!/bin/python
import pandas as pd
inputfile="_counts1"
data = pd.read_csv(inputfile, sep='\t')
df=data.T.sort_index(axis=0).T
df.to_csv('_counts2',sep='\t', index=False)
EOF
python transposer.py
paste <(cut -f1-${sequenceCol} OTU_table.txt) _counts2 <(cut -f${totalcountsCol} OTU_table.txt) > _counts3
mv _counts3 OTU_table.txt
rm -f  transposer.py _counts[1-2] 


#transform files into Excel-ready files
python ${dirpipe}/python/tsvToxlsx.py -i OTU_table.txt
mv  OTU_table.xlsx ./Excel

python ${dirpipe}/python/tsvToxlsx.py -i anchor_table.txt
mv  anchor_table.xlsx ./Excel

#create a fasta file from anchor_table.txt
sequenceColNumber=$(head -n1 anchor_table.txt | sed "s/\t/\n/g" | grep -n "sequence" | cut -d":" -f1)
cut -f1,${sequenceColNumber} anchor_table.txt | awk 'FNR>1' | sed "s/^/>/" | sed "s/\t/\n/" > anchors.fasta
#create a fasta file from OTU_table.txt
sequenceColNumber=$(head -n1 OTU_table.txt | sed "s/\t/\n/g" | grep -n "sequence" | cut -d":" -f1)
cut -f1,${sequenceColNumber} OTU_table.txt | awk 'FNR>1' | sed "s/^/>/" | sed "s/\t/\n/" > OTU.fasta


echo -e "Create microbiomeanalyst.ca-ready files"
mkdir -p ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/MicrobiomeAnalyst/raw_data
cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/MicrobiomeAnalyst/raw_data
ln -nsf ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/OTU_table.txt anchorOTU_table
#OTU TABLE
#counts start right after sequence column and ends before totalcounts column
sequenceCol=$(head -n1 anchorOTU_table | sed "s/\t/\n/g" | grep -n "db_seq" | cut -d":" -f1)
totalcountsCol=$(head -n1 anchorOTU_table | sed "s/\t/\n/g" | grep -n "totalcounts" | cut -d":" -f1)
startCountCols=$((${sequenceCol} + 1))
endCountCols=$((${totalcountsCol} - 1))
cut -f1,${startCountCols}-${endCountCols} anchorOTU_table | sed "1s/^OTU/#NAME/" > OTU_table.txt
#TAXONOMY TABLE
#Taxonomy starts at Domain column and ends at Species column
domainCol=$(head -n1 anchorOTU_table | sed "s/\t/\n/g" | grep -n "Domain" | cut -d":" -f1)
speciesCol=$(head -n1 anchorOTU_table | sed "s/\t/\n/g" | grep -n "Species" | cut -d":" -f1)
cut -f1,${domainCol}-${speciesCol} anchorOTU_table | sed "1s/^OTU/#TAXONOMY/" > Taxonomy_table.txt
#METADATA FILE
sed "1s/^[^\t]*\t/#NAME\t/" ${dir0}/metadata/design.txt > Metadata_file.txt
rm -f anchorOTU_table

#create a sample_data per condition (MA doesn't support NAs) and remove NAs
for cond in ${conditions}
do
	colNumber=$(head -n1 Metadata_file.txt | sed "s/\t/\n/g" | grep -n "${cond}" | cut -d":" -f1)
	cut -f1,${colNumber} Metadata_file.txt | grep -v -e $'\t'"NA$" -e $'\t'"N/A$" > Metadata_${cond}.txt
done




echo -e "Create phyloseq-ready files"
mkdir -p ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/Phyloseq/raw_data
cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/Phyloseq/raw_data
ln -nsf ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/OTU_table.txt anchorOTU_table
#Sample Data
sed "1s/.[^\t]*\t/\t/" ${dir0}/metadata/design.txt > sample_data.txt
#Taxonomy table
cut -f1,${domainCol}-${speciesCol} anchorOTU_table | sed "1s/^OTU\t/\t/" > taxonomy_table.txt
#Dada2 expect Kingdom instead of Domain:
sed -i "1s/Domain/Kingdom/" taxonomy_table.txt
sed -i "1s/^OTU//" taxonomy_table.txt
#OTU table
cut -f1,${startCountCols}-${endCountCols} anchorOTU_table | sed "1s/^OTU\t/\t/" > otu_table.txt
sed -i "1s/^OTU//" otu_table.txt
rm -f anchorOTU_table


echo -e "Create metagenomeSeq-ready files"
mkdir -p ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/metagenomeSeq/raw_data
cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/metagenomeSeq/raw_data
ln -nsf ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/OTU_table.txt anchorOTU_table
#Sample Data
sed "1s/.[^\t]*\t/\t/" ${dir0}/metadata/design.txt > sample_data.txt
#Taxonomy table
cut -f1,${domainCol}-${speciesCol} anchorOTU_table | sed "1s/^OTU\t/\t/" > taxonomy_table.txt
#OTU table
cut -f1,${startCountCols}-${endCountCols} anchorOTU_table | sed "1s/^OTU\t/\t/" > otu_table.txt
rm -f anchorOTU_table
#metagenomeSeq expect the exact same order between the OTU and the taxonomy tables (OTU name order)
head -n1 otu_table.txt > _head
awk 'FNR>1' otu_table.txt | sort -t $'\t' -k1,1 > _body
cat _head _body > otu_table.txt
head -n1 taxonomy_table.txt > _head
awk 'FNR>1' taxonomy_table.txt | sort -t $'\t' -k1,1 > _body
cat _head _body > taxonomy_table.txt
rm -f _head _body




echo -e "Create STAMP-ready files"
mkdir -p ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/STAMP/raw_data
cd ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/STAMP/raw_data
cp ${dir0}/metadata/design.txt metadata.txt
ln -nsf ${dir0}/Results_${expName}_${cutoff}_${AnchorMinBlastIdentity}_${lowCountSeqThreshold}/OTU_table.txt anchorOTU_table
#Let's do a ANCHOR profile (where all sequences are differentiated and thus there is no blanks for unresolved taxonomy). 
#STAMP wants 7 levels so we will remove the Domain and put the OTU name as the last level thanks to paste commanc
phylumCol=$(head -n1 anchorOTU_table | sed "s/\t/\n/g" | grep -n "Phylum" | cut -d":" -f1)
paste <(cut -f${phylumCol}-${speciesCol} anchorOTU_table) <(cut -f1 anchorOTU_table) <(cut -f${startCountCols}-${endCountCols} anchorOTU_table) > profile_ANCHOR.txt
#define a classic profile table where we go from domain to species with blanks where taxonomy is not resolved
cut -f1,${startCountCols}-${endCountCols} anchorOTU_table > otu_table.txt
cut -f1,${domainCol}-${speciesCol} anchorOTU_table > taxonomy_table.txt
head -n1 otu_table.txt > _otuHead
head -n1 taxonomy_table.txt > _taxHead
awk 'FNR>1' otu_table.txt > _otu
awk 'FNR>1' taxonomy_table.txt > _tax
#I'll transpose the tables and remove duplicates
cat <<EOF > remove_duplication_per_line.py
#!/bin/python
import pandas as pd
inputfile="_tax"
data = pd.read_csv(inputfile, sep='\t', header=None, index_col=[0])
data.replace(to_replace="_M..*", value="", inplace=True, regex=True)
df = data.apply(pd.Series.duplicated, axis=1)
data = data.where(~df,"unclassified")
data.to_csv('_noDupl',sep='\t')
EOF
python remove_duplication_per_line.py
#Remove empty when the taxonomy is redundant (e.g. Bacteria	Firmicutes	Clostridia	Clostridiales	empty	Eubacterium	Eubacterium_eligens)
sed  "s/\(\t[A-Za-z]*\)\tunclassified\(\t[A-Z]\)/\1\1\2/" _noDupl  > _newTax
join -1 1 -2 1 -t $'\t' <(sort -t $'\t' -k1,1 _newTax) <(sort -t $'\t' -k1,1 _otu) > _temp1
join -1 1 -2 1 -t $'\t' <(sort -t $'\t' -k1,1 _taxHead) <(sort -t $'\t' -k1,1 _otuHead) > _temp2
#I'll reove the OTU level
cut -f2- _temp1 > _temp3
cut -f2- _temp2 > _temp4
cat _temp4 _temp3 > profile_classic.txt
rm -f _* otu_table.txt taxonomy_table.txt remove_duplication_per_line.py anchorOTU_table



#CLEAN WORKDIR
mkdir -p ${successDir}
touch ${successDir}/compileResults.ok


echo -e "\n---\ncompileResults.sh is exiting normally.\nWell done!"

