#!/bin/bash

#SBATCH --cpus-per-task=15
#SBATCH --mem=40g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --qos=short

hostname
set -ux

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"

TIME=$(date "+%s")

###################################################################################################
#setup-phase

#create path variables
topOPENdir=$OPENdir
locTMP=${TMPdir}mask_repeats/

#create directories
mkdir $locTMP

#load tools
source ${SCRIPTdir}tools

THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))

###################################################################################################

cd ${locTMP}

cp $assemblyFASTA ${locTMP}${assemblyNAME}

#mask repeats using Repeatmasker
RepeatMasker -qq -e rmblast -pa $THREADS -dir ${locTMP} -species "Drosophila melanogaster" ${locTMP}${assemblyNAME}

#zip output to use with kent-script
NAME=$(basename $assemblyFASTA)
cat ${locTMP}${assemblyNAME}.out | gzip > ${locTMP}${assemblyNAME}.out.gz

#generate the indivdual tracks using a script form kentUtils
export LC_COLLATE=C

cp ${HUBdir}/${assemblyNAME}.chrom.sizes ../../${assemblyNAME}.chrom.sizes
asmHubRepeatMasker ${assemblyNAME} ${locTMP}${assemblyNAME}.out.gz ${locTMP}/

###################################################################################################
#column definition file for bigGenePred file
printf 'table bigGenePred
"bigGenePred gene models"
    (
    string chrom;       	  "Reference sequence chromosome or scaffold"
    uint   chromStart;  	  "Start position in chromosome" 
    uint   chromEnd;    	"End position in chromosome"
    string name;        	  "Name or ID of item, ideally both human-readable and unique"
    uint score;         	  "Score (0-1000)"
    char[1] strand;     	  "+ or - for strand"
    uint SmithWatermanScore; "Smith-Waterman score of the match, usually complexity adjusted"
    uint substitutions;     "percent substitutions in matching region compared to the consensus"
    uint deletions;         "percent of bases opposite a gap in the query sequence (deleted bp)"
    uint insertions;     	  "percent of bases opposite a gap in the repeat consensus (inserted bp)"
    string residualSequence;  "no. of bases in query sequence past the ending position of match"
    string class;           "the class of the repeat"
    string family;       	  "the family of the repeat"
    string repeatSTART; 	    "Position in repeat - START"
    string repeatEND;   	    "Position in repeat - END"
    string repeatRESIDUAL;  "Residual bases of repeat"
   )  
' > ${locTMP}bigGenePred.as

###################################################################################################
#add lines to trackDb.txt

#do not continue until file is unblocked by other process
while [[ -f ${TMPdir}wait.txt ]]; do
  sleep 10s
done

#block trackDb from other processes
touch ${TMPdir}wait.txt

#remove old lines if present
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="repeatMasker_" '
  {
    if( $0 !~ NAME ) print
  } 
' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp
mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt

#add repeatMasker relevant tracks
printf "
track repeatMasker_
compositeTrack on
shortLabel RepeatMasker
longLabel Repeating Elements by RepeatMasker
group Annotations
priority 20
visibility dense
type bed 3 .
allButtonPair on
dragAndDrop subTracks
" | sed '/./,$!d' >>${HUBdir}/trackDb.txt

#add entry for all TE-types
PRIO=20
for rmskFILE in ${locTMP}/bbi/*rmsk*; do
  echo $rmskFILE
  TYPE=$(basename $rmskFILE | awk '{n=split($1,splitNAME,/\./); print splitNAME[n-1]}')
  FILE=$(basename $rmskFILE)
  echo $FILE $TYPE

  bigBedToBed $rmskFILE ${locTMP}${TYPE}.bed
  #create final annotation track
  bedToBigBed -tab -as=${locTMP}bigGenePred.as -type=bed6+10 ${locTMP}${TYPE}.bed ${OPENdir}/${assemblyNAME}/${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/${TYPE}.bb


  PRIO=$(( $PRIO + 1 ))
  printf "
  track repeatMasker_${TYPE}
  parent repeatMasker_
  mouseOver \$name family=\$family SWscore=\$SmithWatermanScore
  shortLabel ${TYPE}
  longLabel ${TYPE} Repeating Elements by RepeatMasker
  priority $PRIO
  spectrum on
  maxWindowToDraw 10000000
  colorByStrand 50,50,150 150,50,50
  type bigBed 6 +
  bigDataUrl annotations/${TYPE}.bb
  " >>${HUBdir}/trackDb.txt
done

#unblock file
rm -rf ${TMPdir}wait.txt


###################################################################################################
#finish script

#clean up
if [[ $DEBUG == N ]]; then
  rm -rf $locTMP
fi

#report processing time
PROCESSED_TIME=$(echo -e $(date "+%s") $TIME | awk '{ print ($1-$2)/60 }')
echo "mask repeats=" ${PROCESSED_TIME} >>${LOG}time-log.txt

