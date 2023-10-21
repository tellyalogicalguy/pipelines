#!/bin/bash
set -e
set -u
set -o pipefail

#check number of arguments
if [[ $# -eq 1 ]]
then
	SLRM_NAME=($(basename "$1"))
	PREFIX="$(sed "s:$SLRM_NAME::" <(echo $1))"
else
	echo "Incorrect number of arguments. Supply SLRM_NAME"
	exit 1
fi


echo "SLRM_NAME is: "$SLRM_NAME
echo "in folder: "$PREFIX

REGEX_PATTERN="${SLRM_NAME}.*_slurm-([0-9]+).out"


if [ $(find ${PREFIX} -maxdepth 1 -regextype posix-extended -regex "${PREFIX}${REGEX_PATTERN}" | wc -w) -gt 0 ]
then
	FILENAME=($(find ${PREFIX} -maxdepth 1 -regextype posix-extended -regex "${PREFIX}${REGEX_PATTERN}" | xargs -n 1 basename))
else
	echo "There is no file with the pattern: "$REGEX_PATTERN
	exit 1
fi

echo "There are ${#FILENAME[@]} files with the pattern: "$REGEX_PATTERN

SEFF_OUT_FILE="./seff_outputs/${SLRM_NAME}_seff.log"

for FILENAME in ${FILENAME[@]}
do
	if [[ $FILENAME =~ $REGEX_PATTERN ]]
	then
        	SLRM_ID="${BASH_REMATCH[1]}"
		echo "Running seff for ${SLRM_ID} (${SLRM_NAME})..." | tee -a ${SEFF_OUT_FILE}
		#run seff with passed out slurm variables
		seff $SLRM_ID >> ${SEFF_OUT_FILE}
		echo "SPG_ChIPseq_PE_pipeline wall times:" >> ${SEFF_OUT_FILE}
		echo "-----------------------------------" >> ${SEFF_OUT_FILE}
		cat <(grep 'SPG_ChIPseq_PE_pipeline' ${PREFIX}${FILENAME} || echo "No pipeline outputs") >> ${SEFF_OUT_FILE}
		echo "=================================================" >> ${SEFF_OUT_FILE} 
	fi
done

cat ${SEFF_OUT_FILE}

