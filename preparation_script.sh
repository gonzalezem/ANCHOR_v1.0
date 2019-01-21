#!/bin/bash
set -e
dir0=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
dirReads="$1"
dirExp="$2"
mycond="$3"




echo -e "\n---\nParameters:\n\t- ANCHOR program directory: ${dir0}\n\t- Raw reads location: ${dirReads}\n\t- Experiment directory: ${dirExp}\n\t- Design file location: ${mycond}\n---\n"



cd ${dir0}


find ${dir0}/pipelineScripts -type d -exec chmod +rwx {} \;
find ${dir0}/pipelineScripts -type f -exec chmod +rwx {} \;


ANCHORSCRIPTS="${dir0}/pipelineScripts"
if [ -d "$ANCHORSCRIPTS" ]; then
	echo -e "ANCHOR folder is checked."
else
	echo -e "---\nerror\nBe sure you are inside ANCHOR directory. Here is your present location: ${dir0}\nDownload ANCHOR from there: https://github.com/gonzalezem/ANCHOR/tree/master/ANCHOR_v1.0" && exit 1
fi

if [[ -e "${mycond}" && -s "${mycond}" ]];
then 
	echo -e "Design file is checked."
else 
	echo -e "---\nerror\nYou need a condition file containing at least 2 columns and the first one being named Samples. This is what was given by argument 3: ${mycond} See instructions and come back" && exit 1
fi

if [ -n "${dirExp}" ]; then
	mkdir -p ${dirExp}/metadata
	mkdir -p ${dirExp}/taxonomy/curated
	for i in ${dir0}/taxonomy/curated/*.gzip; do newfile=$(basename $i .gzip); gunzip -c $i > ${dirExp}/taxonomy/curated/$newfile.txt; done
	cp ${mycond} ${dirExp}/metadata/design.txt
	cp ${dir0}/iniFileTemplate/pipe.ini ${dirExp}/metadata/pipe.ini
	cd ${dirExp}
	ln -nsf ${dir0}/pipelineScripts
	ln -nsf ${dir0}/db
	ln -nsf ${dir0}/pipelineScripts/shell/bashMe.sh
else
	echo -e "---\nerror\nArgument 2 (folder from where ANCHOR will be run) was empty. " && exit 1
fi

if [ -d "${dirReads}" ]; then
	cd ${dirExp}
	ln -nsf ${dirReads} raw_reads
else
	echo -e "---\nerror\nWe couldn't find illumina reads directory. Did you add it as the first argument? This is what was given: ${dirReads}" && exit 1
fi


#check mothur 
target="${dir0}/pipelineScripts/mothur"
if find "${target}" -mindepth 1 -print -quit 2>/dev/null | grep -q .; then
	checkMothurName=$(find ${target}/. -type f | grep "mothur" | wc -l)
	if [ ${checkMothurName} -gt 0 ]; then
		echo -e "Mothur program was found in ${target}"
	else
		echo -e "---\nerror\nWe couldn't find mothur within the following folder: ${target}\nThe folder is not empty but does not contain a file called mothur. If you downloaded mothur with another name (like Mothur), please rename it and run this script again." && exit 1
	fi
else
	checkMothur=$(which mothur | wc -l)
	if [ ${checkMothur} -gt 0 ]; then
		echo -e "Mothur program was found in the system"
		mymothur=$(which mothur)
		ln -nsf ${mymothur} ${target}/mothur
	else
		echo -e "---\nerror\nWe couldn't find mothur within the following folder: ${target} or within the system. Did you download it?" && exit 1
	fi
fi



#check usearch9 
target="${dir0}/pipelineScripts/usearch9"
if find "${target}" -mindepth 1 -print -quit 2>/dev/null | grep -q .; then
	checkusearch9Name=$(find ${target}/. -type f | grep "usearch9" | wc -l)
	if [ ${checkusearch9Name} -gt 0 ]; then
		echo -e "usearch9 program was found in ${target}"
	else
		echo -e "---\nerror\nWe couldn't find usearch9 within the following folder: ${target}\nThe folder is not empty but does not contain a file called usearch9. If you downloaded usearch9 with another name (like Usearch9), please rename it and run this script again." && exit 1
	fi
else
	checkusearch9=$(which usearch9 | wc -l)
	if [ ${checkusearch9} -gt 0 ]; then
		echo -e "usearch9 program was found in the system"
		myusearch9=$(which usearch9)
		ln -nsf ${myusearch9} ${target}/usearch9
	else
		echo -e "---\nerror\nWe couldn't find usearch9 within the following folder: ${target} or within the system. Did you download it?" && exit 1
	fi
fi


#check blastn 
checkblastn=$(which blastn | wc -l)
if [ ${checkblastn} -gt 0 ]; then
	echo -e "blastn program was found in the system"
else
	echo -e "---\nerror\nWe couldn't find blastn within the system. Did you download it?" && exit 1
fi




echo -e "\n---\nANCHOR is ready to run! What you needto do now:\n\t1. Go to ${dirExp}/metadata and customize the file pipe.ini to your need\n\t2. When done, go to ${dirExp} and run the command \"bash bashMe.sh\"\n"
