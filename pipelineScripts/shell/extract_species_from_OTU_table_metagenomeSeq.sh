#!/bin/bash
set -e

#extract_species_only_from_OTU_table.sh
#output: otu_table_species.txt and 

cp otu_table.txt _temp
head -n1 _temp > _head
awk 'FNR>1' _temp | grep -e "_[a-z][a-z]*_[0-9]" > _temp1
awk 'FNR>1' _temp | grep -e "_MS_" > _temp2
cat _head _temp1 _temp2 > otu_table_species.txt

cp taxonomy_table.txt _temp
head -n1 _temp > _head
awk 'FNR>1' _temp | grep -e "_[a-z][a-z]*_[0-9]" > _temp1
awk 'FNR>1' _temp | grep -e "_MS_" > _temp2
cat _head _temp1 _temp2 > taxonomy_table_species.txt

rm -f _temp _temp1 _temp2 _head
echo -e "\n---\nNormal exit"

