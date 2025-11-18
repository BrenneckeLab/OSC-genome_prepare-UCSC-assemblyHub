#!/bin/bash

#SBATCH --cpus-per-task=15
#SBATCH --mem=40g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --qos=short
#SBATCH --time=5:00:00

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
locTMP=${TMPdir}ONT_RNA/${SLURM_ARRAY_TASK_ID}/

#create directories
mkdir -p $locTMP

#load tools
source ${SCRIPTdir}tools

#extract current run
ONT_RNAcurr=$(echo $ONT_RNA |tr '~' '\n' | sed -n ${SLURM_ARRAY_TASK_ID}p)
NAME=$(echo $ONT_RNAcurr | tr '/' '\n' | grep RNA)
###################################################################################################

#mapping reads with minimap2
THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))
minimap2 -t $THREADS --secondary=no  -ax splice $assemblyFASTA ${ONT_RNAcurr}/*.gz  | samtools view -bS -@ 6 > ${locTMP}mapped.ONT.bam

samtools sort --output-fmt BAM -o ${HUBdir}annotations/ONT/${NAME}.bam -@ $THREADS ${locTMP}mapped.ONT.bam
samtools index -@ $THREADS ${HUBdir}annotations/ONT/${NAME}.bam

###################################################################################################
#add tracks to trackhub

#do not continue until file is unblocked by other process
while [[ -f ${TMPdir}wait.txt ]]; do
  sleep 10s
done

#block trackDb from other processes
touch ${TMPdir}wait.txt

#remove old lines if present
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="mappedONT_RNA_$NAME" '
  {
    if( $0 ~ NAME ) {a=b; } else print
  }
' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp
mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt

#add tracks to trackDb

printf "

  track mappedONT_RNA_$NAME
  type bam
  shortLabel ONTreads_primaryALN
  longLabel ONT RNAseq reads from $NAME
  maxWindowToDraw 10000000
  maxItems 1000000
  bamColorMode strand
  visibility squish
  parent mappedONT_RNA
  bigDataUrl annotations/ONT/${NAME}.bam

" >>${HUBdir}/trackDb.txt

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
echo "map ONT reads=" ${PROCESSED_TIME} >>${LOG}time-log.txt

