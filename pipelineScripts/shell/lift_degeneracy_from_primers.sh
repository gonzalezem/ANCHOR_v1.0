#!/bin/bash
set -e

: <<'END'
What is expected: a one column file (my_primer_file) with no header.
Example:
TNANACATGCAAGTCGRRCG
CTGCTGCCTYCCGTA

What it will do: provide a file (output.txt) with all the sequence combination of the degenrate primers in input file

COMMAND:
lift_degeneracy_from_primers.sh my_primer_file

OUTPUT: 
output.txt

END

cp $1 input.txt

remove_ambiguity()
{
  SEQUENCE=$1
  AMBINUCL=$2
  rm -f _all_combinations
  if [ ${AMBINUCL} == "N" ];then
 	mod1=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/A/")
 	echo -e  "${mod1}" >> _all_combinations
 	mod2=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/C/")
 	echo -e  "${mod2}" >> _all_combinations
 	mod3=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/G/")
 	echo -e  "${mod3}" >> _all_combinations
 	mod4=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/T/")
 	echo -e  "${mod4}" >> _all_combinations
 fi
 if [ ${AMBINUCL} == "Y" ];then
 	mod2=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/C/")
 	echo -e  "${mod2}" >> _all_combinations
 	mod4=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/T/")
 	echo -e  "${mod4}" >> _all_combinations
 fi
 if [ ${AMBINUCL} == "R" ];then
 	mod1=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/A/")
 	echo -e  "${mod1}" >> _all_combinations
 	mod3=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/G/")
 	echo -e  "${mod3}" >> _all_combinations
 fi
 if [ ${AMBINUCL} == "S" ];then
 	mod2=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/C/")
 	echo -e  "${mod2}" >> _all_combinations
 	mod3=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/G/")
 	echo -e  "${mod3}" >> _all_combinations
 fi
 if [ ${AMBINUCL} == "W" ];then
 	mod1=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/A/")
 	echo -e  "${mod1}" >> _all_combinations
 	mod4=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/T/")
 	echo -e  "${mod4}" >> _all_combinations
 fi
 if [ ${AMBINUCL} == "K" ];then
 	mod3=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/G/")
 	echo -e  "${mod3}" >> _all_combinations
 	mod4=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/T/")
 	echo -e  "${mod4}" >> _all_combinations
 fi
 if [ ${AMBINUCL} == "M" ];then
 	mod1=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/A/")
 	echo -e  "${mod1}" >> _all_combinations
 	mod2=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/C/")
 	echo -e  "${mod2}" >> _all_combinations
 fi
 if [ ${AMBINUCL} == "B" ];then
 	mod2=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/C/")
 	echo -e  "${mod2}" >> _all_combinations
 	mod3=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/G/")
 	echo -e  "${mod3}" >> _all_combinations
 	mod4=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/T/")
 	echo -e  "${mod4}" >> _all_combinations
 fi
 if [ ${AMBINUCL} == "D" ];then
 	mod1=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/A/")
 	echo -e  "${mod1}" >> _all_combinations
 	mod3=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/G/")
 	echo -e  "${mod3}" >> _all_combinations
 	mod4=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/T/")
 	echo -e  "${mod4}" >> _all_combinations
 fi
 if [ ${AMBINUCL} == "H" ];then
 	mod1=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/A/")
 	echo -e  "${mod1}" >> _all_combinations
 	mod2=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/C/")
 	echo -e  "${mod2}" >> _all_combinations
 	mod4=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/T/")
 	echo -e  "${mod4}" >> _all_combinations
 fi
 if [ ${AMBINUCL} == "V" ];then
 	mod1=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/A/")
 	echo -e  "${mod1}" >> _all_combinations
 	mod2=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/C/")
 	echo -e  "${mod2}" >> _all_combinations
 	mod3=$(echo -e ${SEQUENCE} | sed "s/${AMBINUCL}/G/")
 	echo -e  "${mod3}" >> _all_combinations
 fi
}

while read fwd 
do
	echo -e ${fwd} > __lift_degeneracy_from_primers__temp1
	ambiguousNuclList="Y R S W K M B D H V N"
	COUNT=1
	i=1
	for ambi in ${ambiguousNuclList}
	do
		anyhit_fwd=$(cat __lift_degeneracy_from_primers__temp1| grep "${ambi}" | wc -l)
		while [ ${anyhit_fwd} -gt 0 ]
		do
			while read seqToChange
			do
				remove_ambiguity ${seqToChange} ${ambi}
				mv _all_combinations __lift_degeneracy_from_primers__tmp_${COUNT}
				(( COUNT = COUNT + 1))
			done<__lift_degeneracy_from_primers__temp1
			cat __lift_degeneracy_from_primers__tmp_* > __lift_degeneracy_from_primers__temp1
			rm -f __lift_degeneracy_from_primers__tmp_*
			anyhit_fwd=$(cat __lift_degeneracy_from_primers__temp1| grep "${ambi}" | wc -l)
		done
	done
	cp __lift_degeneracy_from_primers__temp1 _allPossible_fwd_${COUNT}
	(( COUNT = COUNT + 1))
	rm -f __lift_degeneracy_from_primers__temp1
done<input.txt
cat _allPossible_fwd_* > output.txt

rm -f __lift_degeneracy_from_primers__temp1 _allPossible_fwd_* input.txt
echo -e "lift_degeneracy_from_primers.sh: ok" && exit 0
