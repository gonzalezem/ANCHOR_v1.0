#!/bin/bash
set -e


: <<'END'
what it does:
Sanity check: Confirms that the counts in the final OTU table and the expected counts from the Unique sequence count table (cf. uniqueContigs.sh), the 2 BLASTn runs (cf. anchorSelection.sh and lowCountSeqSelection.sh) are the same. 
END




dir0="YourMainPath"
iniFile="${dir0}/metadata/pipe.ini"
runDir="${dir0}/run"
successDir="${runDir}/successfulRuns"


echo -e "\n-----------------------------------\nPARAMETERS\n"
echo -e "Main directory: $dir0"
echo -e "Parameter file: ${iniFile}"
echo -e "Database list: $databaseList"
echo -e "Design file: $Conditions"
echo -e "-----------------------------------\n\n"



mkdir -p ${dir0}/countCheck/ERRORS
cd ${dir0}/countCheck/ERRORS
ln -nsf ${runDir}/countSummary/afterBlast_mappingFile.txt
ln -nsf ${runDir}/OTUandCounts/OTU_table.txt
#extract the "totalcounts" column from OTU_table.txt
countCol=$(head -n1 OTU_table.txt | sed "s/\t/\n/g" | grep -n "totalcounts" | cut -d":" -f1)
numberOfSeqpart4Parsed=$(cut -f${countCol} OTU_table.txt |awk 'FNR>1' | paste -sd+ - | bc)
numberOfSeqInMappingFile=$(awk 'FNR>1' afterBlast_mappingFile.txt | wc -l)
if [ ${numberOfSeqpart4Parsed} -eq ${numberOfSeqInMappingFile} ]; then
    echo -e "Number of final sequences has been successfully checked."
    rm -rf ${dir0}/countCheck/ERRORS
    touch ${dir0}/countCheck/OTU_table.ok
else
	echo -e "\n\n------------------------------------------------------------------\nERROR! Check these:\n"
    echo -e "${runDir}/OTUandCounts/OTU_table.txt has:    ${numberOfSeqpart4Parsed} sequences"
    echo -e "${runDir}/countSummary/afterBlast_mappingFile.txt has:    ${numberOfSeqInMappingFile} sequences"
	exit 1
fi



#CLEAN WORKDIR
cd ${dir0}/countCheck
mv ${dir0}/countCheck ${dir0}/run/
mkdir -p ${successDir}
touch ${successDir}/countCheck.ok


echo -e "\n---\ncountCheck.sh is exiting normally.\nWell done!"

