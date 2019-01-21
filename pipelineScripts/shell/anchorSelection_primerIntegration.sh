#!/bin/bash
set -ex

: <<'END'
Goal: Through a 2ndary blast we refine the best hits of the first blast taking ambiguous primers into account (ambiguous primer bases are neutralized in order not to influence the blast output)
END





dir0="YourMainPath"
iniFile="${dir0}/metadata/pipe.ini"
dirpipe="${dir0}/pipelineScripts"
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"
primers="${dir0}/metadata/primers.txt"
expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
databaseList=$(grep "^databaseList" ${iniFile} | cut -d"=" -f2)
AnchorMinBlastIdentity=$(grep "^AnchorMinBlastIdentity" ${iniFile} | cut -d"=" -f2)
wordSizeAnchors=$(grep "^wordSizeAnchors" ${iniFile} | cut -d"=" -f2)
procNumber=$(grep "^procNumber" ${iniFile} | cut -d"=" -f2)
bypass=$(grep "^bypass" $1 | cut -d"=" -f2)
primerSelectionBypass=$(grep "^primerSelectionBypass" $1 | cut -d"=" -f2)


echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Experiment name: ${expName}"
echo -e "main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "Pipeline folder: $dirpipe"
echo -e "Anchor identity and coverage threshold: $AnchorMinBlastIdentity"
echo -e "Database list: $databaseList"
echo -e "-----------------------------------\n\n"


if [ "${primerSelectionBypass}"  == "YES" ]; then
	touch ${successDir}/anchorSelection_primerIntegration.ok
	exit 0
fi



FILE="${primers}"
if [ ! -e "$FILE" ];
then
   echo -e "Primer file $FILE does not exist. You need it to run. Find it and place it here: ${primers}" >&2
   exit 1
fi

mkdir -p ${dir0}/anchorSelection_primerIntegration
cd ${dir0}/anchorSelection_primerIntegration

#PARSE PRIMERS
awk 'FNR>1' ${primers} | sed "s/\t/\n/g" | sort | uniq > __primers
#reverse complement
rm -f __reverseCompl
awk 'FNR>1' ${primers} | sed "s/\t/\n/g" | sort | uniq | while read L; do  echo "$L" | rev | tr "ATGCYRSWKMBVDHN" "TACGRYSWMKVBHDN" ; done>>__reverseCompl
cat __primers __reverseCompl | perl -ne '(!/^\s+$/)&&print' | sed -e "s/H/\[ACT\]/g" | sed -e "s/V/\[ACG\]/g" | sed -e "s/N/\[ACGT\]/g" | sed -e "s/W/\[AT\]/g" | sed -e "s/R/\[AG\]/g" | sed -e "s/Y/\[CT\]/g" | sed -e "s/S/\[GC\]/g" | sed -e "s/K/\[GT\]/g" | sed -e "s/M/\[AC\]/g" | sed -e "s/B/\[CGT\]/g" | sed -e "s/D/\[AGT\]/g" | sed "s/\t/\n/g" | sort | uniq > __allPrimers
#Now do the same for reverse primers
#First, I'll have to reverse the primers as I am going to use the rev command below for convenience
awk 'FNR>1' ${primers} | sed "s/\t/\n/g" | sort | uniq | rev > __rev_primers
cat __reverseCompl | rev > __rev_reverseCompl
cat __rev_primers __rev_reverseCompl | perl -ne '(!/^\s+$/)&&print' | sed -e "s/H/\[ACT\]/g" | sed -e "s/V/\[ACG\]/g" | sed -e "s/N/\[ACGT\]/g" | sed -e "s/W/\[AT\]/g" | sed -e "s/R/\[AG\]/g" | sed -e "s/Y/\[CT\]/g" | sed -e "s/S/\[GC\]/g" | sed -e "s/K/\[GT\]/g" | sed -e "s/M/\[AC\]/g" | sed -e "s/B/\[CGT\]/g" | sed -e "s/D/\[AGT\]/g" | sed "s/\t/\n/g" | sort | uniq  > __rev_allPrimers


for db in ${databaseList}
do
	cribFile="${dir0}/anchorSelection_primerIntegration/${db}/crib_original_seq_primer_transf_seq.txt.zip"
	if [[ -e "${cribFile}" && -s "${cribFile}" ]];then
		echo -e "${cribFile} exists and not empty. Skipping..."
	else

		echo -e "\n---\nParsing ${db} BLAST output"
		mkdir -p ${dir0}/anchorSelection_primerIntegration/${db}
		echo -e "Calculate coverage on query"
		cd ${dir0}/anchorSelection_primerIntegration/${db}
		rm -f blastnRes_anchors_vs_${db}.txt  crib_original_seq_primer_transf_seq.txt ${db}.fasta

		for output in ${dir0}/run/anchorSelection/BLASTn_${db}/blastResults_${db}/output_*.txt
		do
			#I'll remove the subject description (15th column) for now as some hits have so much inside that field that it reaches the limits of awk. I'll include  them later
			cut -f-14,16- ${output} > _temp1
			#extract subjectid and subject description
			cut -f2,15 ${output} | sort | uniq > _temp2
			#limit very very large annotation definition and add database information within the subject id
			cut -d, -f-2 _temp2 | sed "s/^/${db}|/" > _temp3
			#Calculate coverage on query
			awk -v OFS='\t' -v var="${db}" '{ print $1, var"|"$2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, "NoDescription", $4/$13*100, $14, $4/$14*100, $15, $16}' _temp1 > _temp4
			#I may have a coverage >100 when there are gaps within the alignment. In this case I'll lower the coverages to 100% to keep the gap penalty reflected upon the identity value only.
			awk -v OFS='\t' '{if ($15>100.0){$15=100.0} print}' _temp4 > _temp5
			awk -v OFS='\t' '{if ($17>100.0){$17=100.0} print}' _temp5 > _temp4
			#joining with the shortened subject descriptions
			join  -1 2 -2 1 -t $'\t' -a1 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 1.18 1.19 2.2 -e ERROR <(sort -t $'\t' -k2,2 _temp4) <(sort -t $'\t' -k1,1 _temp3) > _temp6 
			#For 16SMicrobial, change 16SMicrobial|gi|1277396189|ref|NR_151900.1| to 16SMicrobial|NR_151900.1
			if [ "${db}" == "16SMicrobial" ]; then
				paste <(cut -f1 _temp6) <(cut -f2 _temp6 | sed "s/|$//" | rev | cut -d"|" -f1 | rev | sed "s/^/16SMicrobial|/") <(cut -f3- _temp6) > _temp7
				mv _temp7 _temp6
			fi
			echo "keep only best hits (potential hits threshold is the rounded floor of the higest id & cov) to ${threshold}% identity and coverage hits"
			awk -F"\t" -v var="${AnchorMinBlastIdentity}" '($3>=var && $15>=var)' _temp6 > _temp7

			#cut: 2-> sseqid  19->sseq 
			cut -f2,19 _temp7 | sort | uniq > _temp8
			#add an index (for accession ids that are represented multiple times with different sequences)
			awk -F'\t' -v OFS='' '{ print "_"NR,"_"$0 }' _temp8 > _temp9
			#Extract all output sequences that have the exact primer
			ln -nsf ${dir0}/anchorSelection_primerIntegration/__allPrimers
			ln -nsf ${dir0}/anchorSelection_primerIntegration/__rev_allPrimers
			rm -f __blastn_with_fwd
			while read primer
			do
				anyhits=$(grep $'\t'"${primer}" _temp9 | wc -l)
				if [ ${anyhits} -gt 0 ];then
					grep $'\t'"${primer}" _temp9 >>__blastn_with_fwd
				fi
			done<__allPrimers
			rm -f __good_blastn_output
			while read primer
			do
				anyhits=$(grep "${primer}$" __blastn_with_fwd | wc -l)
				if [ ${anyhits} -gt 0 ];then
					grep "${primer}$" __blastn_with_fwd >>__good_blastn_output
				fi
			done<__allPrimers
			sort __good_blastn_output | uniq > tmp
			mv tmp __good_blastn_output
			#extract all sequences that did not contain extact primer (i.e. variation within the primer region at non ambiguous sites)
			join -1 1 -2 1 -t $'\t' -v1 <(sort -t $'\t' -k1,1 _temp9) <(sort -t $'\t' -k1,1 __good_blastn_output) > __mismatches_blastn_output

			#Now we have sequences with the exact primer
			#1. we will replace the ambiguous nucleotides with a XXXXX in the db sequences, build a new database
			#2. Same thing with the query sequences 
			#3 blastn again and replace the old blastn stats (before selecting with the primers) with the new blastn stats
			#To do 1. I will grep each primer one by one and replace the positions of the ambiguous by an XXXXX
			#First, let's do the forward sequences
			rm -f _for_new_blast_db_forward
			cp __good_blastn_output __no_hits
			while read primer
			do
				echo -e "${primer}"
				#I copied all the dood sequences into __no_hits to be able to make a check (i.e. if __no_hits is empty, then no need to go further)
				if [[ -e "__no_hits" && -s "__no_hits" ]];then 
					anyhits=$(grep $'\t'"${primer}" __no_hits | wc -l)
					if [ ${anyhits} -gt 0 ];then
						#echo -e "got hits!"
						#We have a hit for a forward primer.
						#we'll replace whatever is in between [ ] by a XXXXX
						newPrimer=$(echo -e "${primer}" | sed "s/\[[ACGT][ACGT]*\]/XXXXX/g")
						#now record all positions of XXXXX within the new primer
						positionAmbiguous=$(echo -e "${newPrimer}" | grep -aob 'XXXXX' | cut -d":" -f1)
						grep $'\t'"${primer}" __no_hits > __hits
						#extract the sequences that did not match to the primer
						anySeqLeft=$(grep -v $'\t'"${primer}" __no_hits | wc -l)
						#echo -e "${anySeqLeft}"
						if [ ${anySeqLeft} -gt 0 ];then
							grep -v $'\t'"${primer}" __no_hits > temp__no_hits
							mv temp__no_hits __no_hits
						else
							echo -e "No more sequence left!"
							rm -f __no_hits
						fi
						cut -f2 __hits > __hits_transformed
						for pos in ${positionAmbiguous}
						do
							((pos = pos + 1))
							sed -i "s/./XXXXX/${pos}" __hits_transformed
						done
						paste <(cut -f1 __hits) __hits_transformed >> _for_new_blast_db_forward
					fi
				fi
			done<__allPrimers
			rm -f _for_new_blast_db
			cp _for_new_blast_db_forward __no_hits
			while read primer
			do
				echo -e "${primer}"
				#I copied all the good sequences into __no_hits to be able to make a check (i.e. if __no_hits is empty, then no need to go further)
				if [[ -e "__no_hits" && -s "__no_hits" ]];then 
					anyhits=$(cat  __no_hits | rev | grep "^${primer}" | wc -l)
					#echo -e "Number of Hits: ${anyhits}"
					if [ ${anyhits} -gt 0 ];then
						#echo -e "got hits!"
						#We have a hit for a forward primer.
						#we'll replace whatever is in between [ ] by a XXXXX
						newPrimer=$(echo -e "${primer}" | sed "s/\[[ACGT][ACGT]*\]/XXXXX/g")
						#now record all positions of XXXXX within the new primer
						positionAmbiguous=$(echo -e "${newPrimer}" | grep -aob 'XXXXX' | cut -d":" -f1)
						cat  __no_hits | rev | grep "^${primer}" > __hits
						#extract the sequences that did not match to the primer
						anySeqLeft=$(cat __no_hits | rev | grep -v "^${primer}" | wc -l)
						#echo -e "${anySeqLeft}"
						if [ ${anySeqLeft} -gt 0 ];then
							cat __no_hits | rev | grep -v "^${primer}" | rev > temp__no_hits
							mv temp__no_hits __no_hits
						else
							echo -e "No more sequence left!"
							rm -f __no_hits
						fi
						cut -f1 __hits > __hits_transformed
						for pos in ${positionAmbiguous}
						do
							((pos = pos + 1))
							sed -i "s/./XXXXX/${pos}" __hits_transformed
						done
						paste __hits_transformed <(cut -f2 __hits) | rev >> _for_new_blast_db
					fi
				fi
			done<__rev_allPrimers
			sort _for_new_blast_db | uniq > tmp
			mv tmp _for_new_blast_db
			#Let's make a crib to still have access to real sequences (and not sequences with A)
			join -1 1 -2 1 -t $'\t' <(sort -t $'\t' -k1,1 __good_blastn_output) <(sort -t $'\t' -k1,1 _for_new_blast_db) | sed "s/XXXXX/A/g" | sed "1s/^/subjectid\toriginal_sequence\tsequence_modified_at_primer_position\n/" >> crib_original_seq_primer_transf_seq.txt
			
			#Make a db sequence fasta
			sed "s/^/>/" _for_new_blast_db | sed "s/\t/\n/" >> ${db}.fasta
		done
		mkdir -p ${dir0}/anchorSelection_primerIntegration/${db}/${db}_index
		mv ${db}.fasta ${dir0}/anchorSelection_primerIntegration/${db}/${db}_index
		cd ${dir0}/anchorSelection_primerIntegration/${db}/${db}_index
		#replace XXXXX with A so blastn Doesn't freak out (edge effect)
		sed -i "s/XXXXX/A/g" ${db}.fasta
		makeblastdb -in ${db}.fasta -dbtype 'nucl' -out ${db} -input_type fasta
		rm -f ${db}.fasta
	fi
done


mkdir -p ${dir0}/anchorSelection_primerIntegration/anchorSequences
anchorCrib="${dir0}/anchorSelection_primerIntegration/anchorSequences/crib_original_AnchorSeq_primer_transf_seq.txt"
if [[ -e "${anchorCrib}" && -s "${anchorCrib}" ]];then
	echo -e "${anchorCrib} exists and not empty. Skipping blastn"
else
	#CHANGE the query sequence with A at primer ambiguous nucleotide positiions
	cd ${dir0}/anchorSelection_primerIntegration/anchorSequences
	ln -nsf ${runDir}/anchorSelection/anchors_fasta.tsv
	ln -nsf ${dir0}/anchorSelection_primerIntegration/__allPrimers
	ln -nsf ${dir0}/anchorSelection_primerIntegration/__rev_allPrimers
	#Let's change every ambiguous nucleotide with an A
	rm -f _for_new_blast_db_forward
	cp anchors_fasta.tsv __no_hits
	
	while read primer
	do
		echo -e "${primer}"
		#I copied all the dood sequences into __no_hits to be able to make a check (i.e. if __no_hits is empty, then no need to go further)
		if [[ -e "__no_hits" && -s "__no_hits" ]];then 
			anyhits=$(grep $'\t'"${primer}" __no_hits | wc -l)
			if [ ${anyhits} -gt 0 ];then
				#echo -e "got hits!"
				#We have a hit for a forward primer.
				#we'll replace whatever is in between [ ] by a XXXXX
				newPrimer=$(echo -e "${primer}" | sed "s/\[[ACGT][ACGT]*\]/XXXXX/g")
				#now record all positions of XXXXX within the new primer
				positionAmbiguous=$(echo -e "${newPrimer}" | grep -aob 'XXXXX' | cut -d":" -f1)
				grep $'\t'"${primer}" __no_hits > __hits
				#extract the sequences that did not match to the primer
				anySeqLeft=$(grep -v $'\t'"${primer}" __no_hits | wc -l)
				#echo -e "${anySeqLeft}"
				if [ ${anySeqLeft} -gt 0 ];then
					grep -v $'\t'"${primer}" __no_hits > temp__no_hits
					mv temp__no_hits __no_hits
				else
					echo -e "No more sequence left!"
					rm -f __no_hits
				fi
				cut -f2 __hits > __hits_transformed
				for pos in ${positionAmbiguous}
				do
					((pos = pos + 1))
					sed -i "s/./XXXXX/${pos}" __hits_transformed
				done
				paste <(cut -f1 __hits) __hits_transformed >> _transformed_anchorSeqs_forward
			fi
		fi
	done<__allPrimers
	rm -f _transformed_anchorSeqs
	cp _transformed_anchorSeqs_forward __no_hits
	while read primer
	do
		echo -e "${primer}"
		#I copied all the good sequences into __no_hits to be able to make a check (i.e. if __no_hits is empty, then no need to go further)
		if [[ -e "__no_hits" && -s "__no_hits" ]];then 
			anyhits=$(cat  __no_hits | rev | grep "^${primer}" | wc -l)
			#echo -e "Number of Hits: ${anyhits}"
			if [ ${anyhits} -gt 0 ];then
				#echo -e "got hits!"
				#We have a hit for a forward primer.
				#we'll replace whatever is in between [ ] by a XXXXX
				newPrimer=$(echo -e "${primer}" | sed "s/\[[ACGT][ACGT]*\]/XXXXX/g")
				#now record all positions of XXXXX within the new primer
				positionAmbiguous=$(echo -e "${newPrimer}" | grep -aob 'XXXXX' | cut -d":" -f1)
				cat  __no_hits | rev | grep "^${primer}" > __hits
				#extract the sequences that did not match to the primer
				anySeqLeft=$(cat __no_hits | rev | grep -v "^${primer}" | wc -l)
				#echo -e "${anySeqLeft}"
				if [ ${anySeqLeft} -gt 0 ];then
					cat __no_hits | rev | grep -v "^${primer}" | rev > temp__no_hits
					mv temp__no_hits __no_hits
				else
					echo -e "No more sequence left!"
					rm -f __no_hits
				fi
				cut -f1 __hits > __hits_transformed
				for pos in ${positionAmbiguous}
				do
					((pos = pos + 1))
					sed -i "s/./XXXXX/${pos}" __hits_transformed
				done
				paste __hits_transformed <(cut -f2 __hits) | rev >> _transformed_anchorSeqs
			fi
		fi
	done<__rev_allPrimers
	#Make a query sequence fasta
	rm -f anchors.fasta
	join -1 1 -2 1 -t $'\t' <(sort -t $'\t' -k1,1 anchors_fasta.tsv) <(sort -t $'\t' -k1,1 _transformed_anchorSeqs) | sed "s/XXXXX/A/g" | sed "1s/^/subjectid\toriginal_sequence\tsequence_modified_at_primer_position\n/" >> crib_original_AnchorSeq_primer_transf_seq.txt
	
	
	sed "s/^/>/" _transformed_anchorSeqs | sed "s/\t/\n/" > anchors.fasta
	#replace XXXXX with A so blastn Doesn't freak out (edge effect)
	sed -i "s/XXXXX/A/g" anchors.fasta
fi


echo -e "Starting new blast"
cd ${dir0}/anchorSelection_primerIntegration
for db in ${databaseList}
do
	echo -e "${db}"
	ZipFile="${dir0}/anchorSelection_primerIntegration/BLASTn_${db}.zip"
	if [[ -e "${ZipFile}" && -s "${ZipFile}" ]];then
		echo -e "${ZipFile} exists and not empty. Skipping blastn"
	else
		mkdir -p ${dir0}/anchorSelection_primerIntegration/BLASTn_${db}/blastResults_${db}
		cd ${dir0}/anchorSelection_primerIntegration/BLASTn_${db}
		ln -nsf ${dir0}/anchorSelection_primerIntegration/anchorSequences/anchors.fasta anchors.fasta
		python ${dirpipe}/python/blastnPrep.py -i anchors.fasta -p ${dir0}/anchorSelection_primerIntegration/BLASTn_${db} -cc 10000 -dp ${dir0}/anchorSelection_primerIntegration/${db} -dn ${db}
		sed -i "s/num_threads 12 -word_size 30 -perc_identity 80 -max_target_seqs 1/num_threads ${procNumber} -word_size ${wordSizeAnchors} -perc_identity ${AnchorMinBlastIdentity}/" ${dir0}/anchorSelection_primerIntegration/BLASTn_${db}/scripts/blastn_vs_${db}.sh
		bash ${dir0}/anchorSelection_primerIntegration/BLASTn_${db}/scripts/blastn_vs_${db}.sh
		rm anchors.fasta
		for output in ${dir0}/anchorSelection_primerIntegration/BLASTn_${db}/blastResults_${db}/output_*.txt
		do
			sed -i "s/\t_[0-9][0-9]*_${db}|/\t/g" ${output}
		done
		cd ${dir0}/anchorSelection_primerIntegration
		zip -r BLASTn_${db}.zip BLASTn_${db}
		rm -rf BLASTn_${db}/fastachunks
		cd ${dir0}/anchorSelection_primerIntegration/${db}
		rm -f _*
		zip crib_original_seq_primer_transf_seq.txt.zip crib_original_seq_primer_transf_seq.txt
	fi
done




##LOST DUDES
#extract the anchors that did have a blast output in the first run but don't have one now (i.e. variation in a(some) non-ambiguous nucleotide base(s) of the primer(s))
cd ${dir0}/anchorSelection_primerIntegration
mkdir -p ${dir0}/anchorSelection_primerIntegration/lost_anchors
cd ${dir0}/anchorSelection_primerIntegration/lost_anchors
for db in ${databaseList}
do
	lostdudesFasta="${dir0}/anchorSelection_primerIntegration/lost_anchors/lost_dudes.fasta"
	if [[ -e "${lostdudesFasta}" && -s "${lostdudesFasta}" ]];then
		echo -e "${lostdudesFasta} exists and not empty. Skipping..."
	else
		cd ${dir0}/anchorSelection_primerIntegration/lost_anchors
		previousBlastn="${dir0}/run/anchorSelection/BLASTn_${db}/blastResults_${db}"
		newBlastn="${dir0}/anchorSelection_primerIntegration/BLASTn_${db}/blastResults_${db}"
	
		rm -f __previous_anchorList __new_anchorList
		for output in ${previousBlastn}/output_*.txt
		do
			cut -f1 ${output} >> __previous_anchorList
		done
		sort __previous_anchorList | uniq > tmp
		mv tmp __previous_anchorList
		for output in ${newBlastn}/output_*.txt
		do
			cut -f1 ${output} >> __new_anchorList
		done
		sort __new_anchorList | uniq > tmp
		mv tmp __new_anchorList
		join -v1 __previous_anchorList __new_anchorList > __lostDudes
		lostdudesFile="${dir0}/anchorSelection_primerIntegration/lost_anchors/__lostDudes"
		if [[ -e "${lostdudesFile}" && -s "${lostdudesFile}" ]];then
			echo -e "${lostdudesFile} exists and not empty."
	
			ln -nsf ${dir0}/anchorSelection_primerIntegration/__allPrimers
			ln -nsf ${dir0}/anchorSelection_primerIntegration/__rev_allPrimers
			cat __allPrimers __rev_allPrimers > _primers
			rm -f __lost_query_information
			while read lostDude
			do
				lostSeq=$(grep "^${lostDude}"$'\t' ${runDir}/anchorSelection/anchors_fasta.tsv | cut -f2)
				while read primer
				do
					anyhit_fwd=$(echo -e "${lostSeq}" | grep "^${primer}" | wc -l)
					if [ ${anyhit_fwd} -gt 0 ];then
						forwardPrim=$(echo -e ${primer} | sed -e "s/\[ACT\]/H/g" | sed -e "s/\[ACG\]/V/g" | sed -e "s/\[ACGT\]/N/g" | sed -e "s/\[AT\]/W/g" | sed -e "s/\[AG\]/R/g" | sed -e "s/\[CT\]/Y/g" | sed -e "s/\[GC\]/S/g" | sed -e "s/\[GT\]/K/g" | sed -e "s/\[AC\]/M/g" | sed -e "s/\[CGT\]/B/g" | sed -e "s/\[AGT\]/D/g")
						linux_primer_fwd="${primer}"
					fi
					anyhit_rev=$(echo -e "${lostSeq}" | rev | grep "^${primer}" | wc -l)
					if [ ${anyhit_rev} -gt 0 ];then
						revPrim=$(echo -e ${primer} | sed -e "s/\[ACT\]/H/g" | sed -e "s/\[ACG\]/V/g" | sed -e "s/\[ACGT\]/N/g" | sed -e "s/\[AT\]/W/g" | sed -e "s/\[AG\]/R/g" | sed -e "s/\[CT\]/Y/g" | sed -e "s/\[GC\]/S/g" | sed -e "s/\[GT\]/K/g" | sed -e "s/\[AC\]/M/g" | sed -e "s/\[CGT\]/B/g" | sed -e "s/\[AGT\]/D/g" | rev)
						linux_primer_rev=$(echo -e "${primer}" | rev | sed "s/\]/__/g" | sed "s/\[/\]/g" | sed "s/__/\[/g")
					fi
				done<_primers
				echo -e  "${lostDude}\t${forwardPrim}\t${revPrim}\t${lostSeq}\t${linux_primer_fwd}\t${linux_primer_rev}" >> __lost_query_information
				#ow extarct lost subject ids
				for output in ${previousBlastn}/output_*.txt
				do
					cut -f1,2 ${output} | grep "^${lostDude}"$'\t' | cut -f2 >> __lost_subject_ids
				done
			done<__lostDudes
			
			
			#Build all possible sequences from degenrate primers and lost anchors
			cut -f2 __lost_query_information | sort | uniq > forwardPrimers
			cut -f3 __lost_query_information | sort | uniq > reversePrimers
			#etarct all possible forward primers
			bash ${dirpipe}/shell/lift_degeneracy_from_primers.sh forwardPrimers
			mv output.txt allPossible_fwd
			#etarct all possible reverse primers
			bash ${dirpipe}/shell/lift_degeneracy_from_primers.sh reversePrimers
			mv output.txt allPossible_rev
			#extract between primers seqeunces from all the lost anchors
			rm -f __lost_query_information2
			while read anchor fwd reverse lostSeq linux_fwd linux_rev
			do
				btwn_primers=$(echo -e "${lostSeq}" | sed "s/^${linux_fwd}//" | sed "s/${linux_rev}$//")
				echo -e  "${anchor}\t${fwd}\t${reverse}\t${lostSeq}\t${btwn_primers}" >> __lost_query_information2
			done<__lost_query_information
			#Now that we have it, build all possible sequences
			
			
			while read anchor fwd reverse lostSeq btwn_primers
			do
				echo -e "${anchor}\t${btwn_primers}" > _temp1
				sed "s/^/${anchor}\t/" allPossible_fwd > _temp2
				sed "s/^/${anchor}\t/" allPossible_rev > _temp3
				join -1 1 -2 1 -t $'\t' -a 1 -a 2 -o 1.1 2.2 1.2 _temp1 _temp2 > _temp4
				join -1 1 -2 1 -t $'\t' -a 1 -a 2 -o 1.1 1.2 1.3 2.2 _temp4 _temp3 > _temp5
				#add an index to the anchor name
				awk -F'\t' -v OFS='\t' '{ print $0,"_TOREMOVE"NR }'  _temp5 > _temp6
				awk -F"\t" '{print $1$5 "\t" $2$3$4 }' _temp6 > _temp7
				sed "s/^/>/" _temp7 | sed "s/\t/\n/" > fasta_${anchor}
			done<__lost_query_information2
			cat fasta_* > lost_dudes.fasta
			rm -f reversePrimers forwardPrimers _primers __rev_allPrimers __allPrimers fasta_* _temp[1-7]
		else
			echo -e "${lostdudesFile} is empty. Skipping ${db} database"
		fi
	fi	

	#RERUN BLAST for the lost dudes
	mkdir -p ${dir0}/anchorSelection_primerIntegration/lost_anchors/BLASTn_${db}/blastResults_${db}
	cd ${dir0}/anchorSelection_primerIntegration/lost_anchors/BLASTn_${db}
	ZipFile="${dir0}/anchorSelection_primerIntegration/lost_anchors/BLASTn_${db}.zip"
	if [[ -e "${ZipFile}" && -s "${ZipFile}" ]];then
		echo -e "${ZipFile} exists and not empty. Skipping blastn on ${db}"
	else
		lostdudesFasta="${dir0}/anchorSelection_primerIntegration/lost_anchors/lost_dudes.fasta"
		if [[ -e "${lostdudesFasta}" && -s "${lostdudesFasta}" ]];then
			newBlastn="${dir0}/anchorSelection_primerIntegration/BLASTn_${db}/blastResults_${db}"
			ln -nsf ${dir0}/anchorSelection_primerIntegration/lost_anchors/lost_dudes.fasta anchors.fasta
			python ${dirpipe}/python/blastnPrep.py -i anchors.fasta -p ${dir0}/anchorSelection_primerIntegration/lost_anchors/BLASTn_${db} -cc 10000 -dp ${dir0}/db -dn ${db}
			sed -i "s/num_threads 12 -word_size 30 -perc_identity 80 -max_target_seqs 1/num_threads ${procNumber} -word_size ${wordSizeAnchors} -perc_identity ${AnchorMinBlastIdentity}/" ${dir0}/anchorSelection_primerIntegration/lost_anchors/BLASTn_${db}/scripts/blastn_vs_${db}.sh
			echo -e "Running a last blast for the lost anchor sequences in ${db} database"
			bash ${dir0}/anchorSelection_primerIntegration/lost_anchors/BLASTn_${db}/scripts/blastn_vs_${db}.sh
			rm anchors.fasta
			#remove the index tag introduced earlier (ex: tag is _TOREMOVE1 in 68F_R357_1_1016_TOREMOVE1)
			for output in ${dir0}/anchorSelection_primerIntegration/lost_anchors/BLASTn_${db}/blastResults_${db}/output_*.txt
			do
				sed -i "s/_TOREMOVE[0-9][0-9]*//" ${output}
			done
			#Add the output to the previous blastn output folder so they will be scanned together in the next script
			startNumber=$(ls -1 ${dir0}/anchorSelection_primerIntegration/BLASTn_${db}/blastResults_${db}/output_*.txt | rev | cut -d"/" -f1 | cut -d"_" -f1 | rev | cut -d"." -f1 | sort -nr | head -n1)
			(( startNumber = startNumber + 1 ))
			prefix=$(ls -1 ${dir0}/anchorSelection_primerIntegration/BLASTn_${db}/blastResults_${db}/output_*.txt | head -n1 | rev | cut -d"/" -f1 | cut -d"." -f2- | sed "s/[0-9][0-9]*//" | rev )
			for output in ${dir0}/anchorSelection_primerIntegration/lost_anchors/BLASTn_${db}/blastResults_${db}/output_*.txt
			do
				cp ${output} ${newBlastn}/${prefix}${startNumber}.txt
				(( startNumber = startNumber + 1 ))
			done
			cd ${dir0}/anchorSelection_primerIntegration/lost_anchors
			zip -r BLASTn_${db}.zip BLASTn_${db}
			rm -rf BLASTn_${db}/fastachunks
		fi
	fi
done

cd ${dir0}/anchorSelection_primerIntegration/anchorSequences
zip crib_original_AnchorSeq_primer_transf_seq.txt.zip crib_original_AnchorSeq_primer_transf_seq.txt
zip anchors.fasta.zip anchors.fasta
rm -f anchors_fasta.tsv _* anchors.fasta crib_original_AnchorSeq_primer_transf_seq.txt
cd ${dir0}/anchorSelection_primerIntegration
rm -f _*

#CLEAN WORKDIR
mv ${dir0}/anchorSelection_primerIntegration ${dir0}/run/
mkdir -p ${successDir}
touch ${successDir}/anchorSelection_primerIntegration.ok

echo -e "\n---\nanchorSelection_primerIntegration exiting normally.\nWell done!"