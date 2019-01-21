#!/bin/bash
set -e

: <<'END'
What you need: a taxonomy folder within the main directory. This folder should contain 2 other folders: nt (with accession2taxid.txt.gz, names.dmp.gz and nodes.dmp.gz) and curated (with all databases parsed taxonomy into 7 columns; Ex: nt_curated.txt)
Goal: Parse the anchor sequences vs. public repositories blast output: extract hits with the percentage identity and the coverage (alignment length / contig length) above a given threshold.
END


dir0="YourMainPath"
iniFile="${dir0}/metadata/pipe.ini"
dirpipe="${dir0}/pipelineScripts"
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"
expName=$(grep "^expName" ${iniFile} | cut -d"=" -f2)
databaseList=$(grep "^databaseList" ${iniFile} | cut -d"=" -f2)
curatedTax=${dir0}/taxonomy/curated
ntNCBITaxonomyFiles=${dir0}/taxonomy/nt


echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Experiment name: ${expName}"
echo -e "Main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "Pipeline folder: $dirpipe"
echo -e "Database list: $databaseList"
echo -e "Parsed taxonomy files folder: $curatedTax"
echo -e "nt NCBI Taxonomy files folder: $ntNCBITaxonomyFiles"
echo -e "-----------------------------------\n\n"

###CHECK FOLDER
checknt=$(echo -e ${databaseList} | grep "nt" | wc -l)
if [ ${checknt} -gt 0 ]
then
	if [ ! -d "$ntNCBITaxonomyFiles" ]; then
  		echo -e "Couldn't find the NCBI nt taxonomy files (in ${ntNCBITaxonomyFiles}).\nDownload:\n\t ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz\n\t ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz\n unzip them into ${ntNCBITaxonomyFiles}"
	fi
fi


mkdir -p ${curatedTax}

for db in ${databaseList}
do
	touch ${curatedTax}/${db}_curated.txt
	mkdir -p ${dir0}/AnchorTaxonomyParser/${db}
	cd ${dir0}/AnchorTaxonomyParser/${db}
	#remove the databse tag before the accession ID to recover original accession ID
	sed "s/${db}|//" ${runDir}/parseAnchorBlastOutput/${db}/${db}_speciesID_list.txt | sort | uniq > ${db}_speciesID_list.txt

	if [ "${db}" == "nt" ]; then
		FILE="${db}_taxonomy.ok"
		if [ -e "$FILE" ];
		then
			echo -e "Skipping nt"
		else
			#this script will create a file nt_taxonomy_parsed.dat with 3 columns: accessionNumber	TaxID	Taxonomy
			python ${dirpipe}/python/NCBI_accession2tax_nt.py -i ${db}_speciesID_list.txt -t ${ntNCBITaxonomyFiles}
		
			#we will replace orginal taxonomy with nr parsed taxonomy (7 columns). If some accession number were not parsed, these will have NOTFOUND instead
			join -1 2 -2 1 -t $'\t' -a 1 -o 1.1 1.2 2.2 -e NOTFOUND <(awk 'FNR>1' ${db}_taxonomy_parsed.dat | sort -t $'\t' -k2,2) <(sort -t $'\t' -k1,1 ${curatedTax}/${db}_curated.txt) > ${db}_curated_parsed.txt
			checkNotFound=$(grep "NOTFOUND" ${db}_curated_parsed.txt | wc -l)
			if [ ${checkNotFound} -gt 0 ]; then
				grep "NOTFOUND" ${db}_curated_parsed.txt > list_TOPARSE.txt
				cut -f2- ${db}_taxonomy_parsed.dat | awk 'FNR>1' | sort | uniq > uniq_taxid.txt
				join -1 1 -2 2 -t $'\t' -a 2 -o 2.2 1.2 <(sort -t $'\t' -k1,1 uniq_taxid.txt) <(sort -t $'\t' -k2,2 list_TOPARSE.txt) | sort | uniq | sort -t $'\t' -k1,1 > ${db}_taxonomy_TOPARSE.txt
				rm -f ${db}_curated_parsed.txt uniq_taxid.txt list_TOPARSE.txt
				echo -e "\n---\nNow parse Taxonomy manually. You need 7 taxons per hit.\n 1.REVIEW ${db}_curated_parsed.txt and change it so it has 7 taxon per line (Kingdom, Phylum, Class, Order, Family, Genus, Species)\n 2. Add it to ${curatedTax}//${db}_curated.txt (Caution: taxonomy needs to be consistent between records wihtin the same database but also between databses) 3. Rerun pipeline"
				exit 1
			else
				cut -f1,3 ${db}_curated_parsed.txt > _temp
				mv _temp ${db}_curated_parsed.txt
				rm -f uniq_taxid.txt list_TOPARSE.txt ${db}_taxonomy_parsed.dat ${db}_speciesID_list.txt
				touch ${db}_taxonomy.ok
			fi
		sed -i "s/^/${db}|/" ${dir0}/AnchorTaxonomyParser/${db}/${db}_curated_parsed.txt
		fi
	else
		FILE="${db}_taxonomy.ok"
		if [ -e "$FILE" ];
		then
			echo -e "Skipping ${db}"
		else
			cd ${dir0}/AnchorTaxonomyParser/${db}
			# join with the curated 16SMicrobial taxonomy database
			join -1 1 -2 1 -t $'\t' -a 1 -o 1.1 2.2 -e NOTFOUND <(sort ${db}_speciesID_list.txt) <(sort -t $'\t' -k1,1 ${curatedTax}/${db}_curated.txt) > ${db}_curated_parsed.txt
			checkNotFound=$(grep "NOTFOUND" ${db}_curated_parsed.txt | wc -l)
			if [ ${checkNotFound} -gt 0 ]; then
				grep "NOTFOUND" ${db}_curated_parsed.txt > list_TOPARSE.txt
				set +e
				join -1 1 -2 1 -t $'\t' -a 2 -o 2.1 1.2 <(sort -t $'\t' -k1,1 /media/emmanuel/storage1/db/16S/${db}TaxonomyCrib.txt) <(sort -t $'\t' -k2,2 list_TOPARSE.txt) | sort | uniq | sort -t $'\t' -k2,2 > ${db}_taxonomy_TOPARSE.txt
				echo -e "\n---\nThe taxonomy file is incomplete for ${db} database. Do this:\n1. Download the database:\n\ta. NCBI16S microbial :ftp://ftp.ncbi.nlm.nih.gov/blast/db/16SMicrobial.tar.gz\n\tb. SILVA: https://www.arb-silva.de/fileadmin/silva_databases\n\tc. RDP: https://rdp.cme.msu.edu/download/current_Bacteria_unaligned.fa.gz\n\td. http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/current_GREENGENES_gg16S_unaligned.fasta.gz\n2. Extract taxonomy into 7 taxons  7 taxon per line (Kingdom, Phylum, Class, Order, Family, Genus, Species)\n3. Add it to ${curatedTax}//${db}_curated.txt (Caution: taxonomy needs to be consistent between records wihtin the same database but also between databses)\n4. Rerun pipeline"
				rm -f list_TOPARSE.txt
				exit 1
			else
				rm -f ${db}_speciesID_list.txt
				touch ${db}_taxonomy.ok
			fi
			sed -i "s/^/${db}|/" ${dir0}/AnchorTaxonomyParser/${db}/${db}_curated_parsed.txt
		fi
	fi
done

#CLEAN WORKDIR
cd ${dir0}/AnchorTaxonomyParser/
mv ${dir0}/AnchorTaxonomyParser ${dir0}/run/
mkdir -p ${successDir}
touch ${successDir}/AnchorTaxonomyParser.ok


echo -e "\n---\nAnchorTaxonomyParser.sh is exiting normally.\nWell done!"