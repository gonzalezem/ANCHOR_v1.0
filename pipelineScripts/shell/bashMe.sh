#!/bin/bash
set -e
dir0=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
pipeIniFile="${dir0}/metadata/pipe.ini"
dirpipe="${dir0}/pipelineScripts/shell"

mkdir -p ${dir0}/run/logs
mkdir -p ${dir0}/run/successfulRuns
successDir="${dir0}/run/successfulRuns"






######################################   PRIMER REMOVAL   #################################################################
FILE="${successDir}/primerSelection.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping PRIMER SELECTION step!"
else
    echo -e "\n-----\nIn 5 seconds: Running PRIMER SELECTION step"
    sleep 5
    echo -e "Process\tDuration (minutes)" > ${dir0}/run/logs/runtime.log
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/primerSelection.sh > ${dir0}/run/primerSelection.sh
    bash ${dir0}/run/primerSelection.sh ${pipeIniFile} > ${dir0}/run/logs/primerSelection.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "PRIMER SELECTION\t${runtime}" >> ${dir0}/run/logs/runtime.log

fi



######################################   MERGE READS  #################################################################
FILE="${successDir}/readsMerge.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping MERGING READS step!"
else
    echo -e "\n-----\nIn 5 seconds: Running MERGING READS step"
    sleep 5
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/readsMerge.sh > ${dir0}/run/readsMerge.sh
    bash ${dir0}/run/readsMerge.sh ${pipeIniFile} > ${dir0}/run/logs/readsMerge.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "MERGING READS\t${runtime}" >> ${dir0}/run/logs/runtime.log
fi



######################################   UNIQUE CONTIGS  #################################################################
FILE="${successDir}/uniqueContigs.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping UNIQUE CONTIGS step!"
else
    echo -e "\n-----\nIn 5 seconds: Running UNIQUE CONTIGS step"
    sleep 5
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/uniqueContigs.sh > ${dir0}/run/uniqueContigs.sh
    bash ${dir0}/run/uniqueContigs.sh ${pipeIniFile} > ${dir0}/run/logs/uniqueContigs.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "UNIQUE CONTIGS\t${runtime}" >> ${dir0}/run/logs/runtime.log
fi



######################################   ANCHOR SELECTION #################################################################
FILE="${successDir}/anchorSelection.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping ANCHOR SELECTION step!"
else
    echo -e "\n-----\nIn 5 seconds: Running ANCHOR SELECTION step"
    sleep 5
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/anchorSelection.sh > ${dir0}/run/anchorSelection.sh
    bash ${dir0}/run/anchorSelection.sh ${pipeIniFile} > ${dir0}/run/logs/anchorSelection.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "ANCHOR SELECTION\t${runtime}" >> ${dir0}/run/logs/runtime.log
fi


######################################   ANCHOR SELECTION WITH PRIMER INTEGRATION #################################################################
FILE="${successDir}/anchorSelection_primerIntegration.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping ANCHOR SELECTION WITH PRIMER INTEGRATION step!"
else
    echo -e "\n-----\nIn 5 seconds: Running ANCHOR SELECTION WITH PRIMER INTEGRATION step"
    sleep 5
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/anchorSelection_primerIntegration.sh > ${dir0}/run/anchorSelection_primerIntegration.sh
    bash ${dir0}/run/anchorSelection_primerIntegration.sh ${pipeIniFile} > ${dir0}/run/logs/anchorSelection_primerIntegration.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "ANCHOR SELECTION WITH PRIMER INTEGRATION\t${runtime}" >> ${dir0}/run/logs/runtime.log
fi



######################################   LOW COUNT SEQUENCES SELECTION #################################################################
FILE="${successDir}/lowCountSeqSelection.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping LOW COUNT SEQUENCES SELECTION step!"
else
    echo -e "\n-----\nIn 5 seconds: Running LOW COUNT SEQUENCES SELECTION step"
    sleep 5
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/lowCountSeqSelection.sh > ${dir0}/run/lowCountSeqSelection.sh
    bash ${dir0}/run/lowCountSeqSelection.sh ${pipeIniFile} > ${dir0}/run/logs/lowCountSeqSelection.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "LOW COUNT SEQUENCES SELECTION\t${runtime}" >> ${dir0}/run/logs/runtime.log
fi




######################################   PARSE ANCHOR BLAST OUTPUT #################################################################
FILE="${successDir}/parseAnchorBlastOutput.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping PARSE ANCHOR BLAST OUTPUT step!"
else
    echo -e "\n-----\nIn 5 seconds: Running PARSE ANCHOR BLAST OUTPUT step"
    sleep 5
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/parseAnchorBlastOutput.sh > ${dir0}/run/parseAnchorBlastOutput.sh
    bash ${dir0}/run/parseAnchorBlastOutput.sh ${pipeIniFile} > ${dir0}/run/logs/parseAnchorBlastOutput.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "PARSE ANCHOR BLAST OUTPUT\t${runtime}" >> ${dir0}/run/logs/runtime.log
fi



######################################   PARSE TAXONOMY  #################################################################
FILE="${successDir}/AnchorTaxonomyParser.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping PARSE TAXONOMY step!"
else
    echo -e "\n-----\nIn 5 seconds: Running PARSE TAXONOMY step"
    sleep 5
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/AnchorTaxonomyParser.sh > ${dir0}/run/AnchorTaxonomyParser.sh
    bash ${dir0}/run/AnchorTaxonomyParser.sh ${pipeIniFile} > ${dir0}/run/logs/AnchorTaxonomyParser.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "PARSE TAXONOMY\t${runtime}" >> ${dir0}/run/logs/runtime.log
fi


######################################   MERGE BLASTN AND TAXONOMY   #################################################################
FILE="${successDir}/mergeBlastnTaxonomy.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping MERGE BLASTN AND TAXONOMY step!"
else
    echo -e "\n-----\nIn 5 seconds: Running MERGE BLASTN AND TAXONOMY step"
    sleep 5
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/mergeBlastnTaxonomy.sh > ${dir0}/run/mergeBlastnTaxonomy.sh
    bash ${dir0}/run/mergeBlastnTaxonomy.sh ${pipeIniFile} > ${dir0}/run/logs/mergeBlastnTaxonomy.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "MERGE BLASTN AND TAXONOMY\t${runtime}" >> ${dir0}/run/logs/runtime.log
fi



######################################   ANCHOR PARSER  #################################################################
FILE="${successDir}/anchorParser.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping ANCHOR PARSER step!"
else
    echo -e "\n-----\nIn 5 seconds: Running ANCHOR PARSER step"
    sleep 5
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/anchorParser.sh > ${dir0}/run/anchorParser.sh
    bash ${dir0}/run/anchorParser.sh ${pipeIniFile} > ${dir0}/run/logs/anchorParser.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "ANCHOR PARSER\t${runtime}" >> ${dir0}/run/logs/runtime.log
fi

######################################   COUNT SUMMARY  #################################################################
FILE="${successDir}/countSummary.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping COUNT SUMMARY step!"
else
    echo -e "\n-----\nIn 5 seconds: Running COUNT SUMMARY step"
    sleep 5
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/countSummary.sh > ${dir0}/run/countSummary.sh
    bash ${dir0}/run/countSummary.sh ${pipeIniFile} > ${dir0}/run/logs/countSummary.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "COUNT SUMMARY\t${runtime}" >> ${dir0}/run/logs/runtime.log
fi



######################################   CHIMERA FLAG  #################################################################
FILE="${successDir}/chimeraFlag.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping CHIMERA FLAG step!"
else
    echo -e "\n-----\nIn 5 seconds: Running CHIMERA FLAG step"
    sleep 5
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/chimeraFlag.sh > ${dir0}/run/chimeraFlag.sh
    bash ${dir0}/run/chimeraFlag.sh ${pipeIniFile} > ${dir0}/run/logs/chimeraFlag.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "CHIMERA FLAG\t${runtime}" >> ${dir0}/run/logs/runtime.log
fi



######################################   ESTIMATE OTUs AND COUNTS  #################################################################
FILE="${successDir}/OTUandCounts.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping ESTIMATE OTUs AND COUNTS step!"
else
    echo -e "\n-----\nIn 5 seconds: Running ESTIMATE OTUs AND COUNTS step"
    sleep 5
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/OTUandCounts.sh > ${dir0}/run/OTUandCounts.sh
    bash ${dir0}/run/OTUandCounts.sh ${pipeIniFile} > ${dir0}/run/logs/OTUandCounts.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "ESTIMATE OTUs AND COUNTS\t${runtime}" >> ${dir0}/run/logs/runtime.log
fi



######################################   COUNT CHECK  #################################################################
FILE="${successDir}/countCheck.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping COUNT CHECK step!"
else
    echo -e "\n-----\nIn 5 seconds: Running COUNT CHECK step"
    sleep 5
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/countCheck.sh > ${dir0}/run/countCheck.sh
    bash ${dir0}/run/countCheck.sh ${pipeIniFile} > ${dir0}/run/logs/countCheck.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "COUNT CHECK\t${runtime}" >> ${dir0}/run/logs/runtime.log
fi



######################################   COMPILE RESULTS  #################################################################
FILE="${successDir}/compileResults.ok"
if [ -e "${FILE}" ];then
    echo -e "\n-----\nSkipping COMPILE RESULTS step!"
else
    echo -e "\n-----\nIn 5 seconds: Running COMPILE RESULTS step"
    sleep 5
    start=`date +%s`
    sed "s#^dir0=..*#dir0=\"${dir0}\"#" ${dirpipe}/compileResults.sh > ${dir0}/run/compileResults.sh
    bash ${dir0}/run/compileResults.sh ${pipeIniFile} > ${dir0}/run/logs/compileResults.log 2>&1
    runtime=$((($(date +%s)-$start)/60))
    echo -e "COMPILE RESULTS\t${runtime}" >> ${dir0}/run/logs/runtime.log
fi

cp ${dir0}/run/logs/runtime.log Results_*/Summary

echo -e "ANCHOR is over. Well done!"
