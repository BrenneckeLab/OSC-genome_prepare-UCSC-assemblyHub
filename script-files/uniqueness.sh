#!/bin/bash

#SBATCH --cpus-per-task=25
#SBATCH --mem=80g
#SBATCH --partition=c
#SBATCH -e "%x.e.%A-%a.txt"
#SBATCH -o "%x.o.%A-%a.txt"
#SBATCH --qos=short
#SBATCH --time=8:00:00


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
locTMPtop=${TMPdir}uniqueness/
locTMP=${locTMPtop}${SLURM_ARRAY_TASK_ID}

#create directories
mkdir -p $locTMP

#load tools
source ${SCRIPTdir}tools

THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))

###################################################################################################
###################################################################################################

#generate reads from genome
if [[ ! -s ${locTMP}uniq-reads.txt || $FORCE == Y ]]; then 
  parallel -j 2 --linebuffer "
    if [[ {} == \"sense\" ]]; then COMMAND=\"seqkit seq\"; else COMMAND=\"seqkit seq --reverse --complement\"; fi

    \$COMMAND $assemblyFASTA  |
      seqkit sliding --threads 4 -s 1 -W $SLURM_ARRAY_TASK_ID | 
      seqkit fx2tab  " ::: sense antisense | 
    mawk -v OFS="\t" '{
      X[$2]+=1
    }
    END{
      for(SEQ in X){
        if(X[SEQ] == 1 ){
          print SEQ
        }
      }
    }'  > ${locTMP}uniq-reads.txt
fi

#generate bowtie index
if [[ ! -s ${locTMP}index.1.ebwt ]]; then
  bowtieBuild --threads $THREADS --noref $assemblyFASTA ${locTMP}index
fi

###################################################################################################
#do not continue until file is unblocked by other process
while [[ -f ${TMPdir}wait.txt ]]; do
  SLEEPtime=$(shuf -i 1-60 -n 1)
  sleep ${SLEEPtime}s
done

#block trackDb from other processes
touch ${TMPdir}wait.txt

#remove old lines if present
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="uniquness_${SLURM_ARRAY_TASK_ID}mer" '
  {
    if( $0 !~ NAME ) print
  }
' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp

mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt
rm -rf ${TMPdir}wait.txt

###################################################################################################
if [[ $SLURM_ARRAY_TASK_ID -eq 25 ]]; then PRIO=60; elif [[ $SLURM_ARRAY_TASK_ID -eq 50 ]]; then PRIO=64; else PRIO=67;fi
for MM in 0 1 2 3; do
  time bowtie --best --strata -m 1 -rS -v $MM --threads $THREADS -x ${locTMP}index ${locTMP}uniq-reads.txt | 
  samtools view -bS | 
  bedtools bamtobed -i - |
  mawk -v OFS="\t" '{
    if($6 == "+"){
      $3=$2+1
    }else{
      $2=$3-1
    }
    print
  }' |
  LC_COLLATE=C sort -k1,1 -k2,2n --parallel=$THREADS -S60G |
  mawk -v OFS="\t" -v TMP=$locTMP '{
    if($6=="+"){
      print > TMP "sense.bed"
    }else{
      print > TMP "antisense.bed"
    }
  }' > ${locTMP}${MM}.bed

  bedtools merge -i ${locTMP}sense.bed | mawk -v OFS="\t" '{print $0,$1":"$2"-"$3":!:+",0,"+"}'> ${locTMP}${MM}.bed
  bedtools merge -i ${locTMP}antisense.bed | mawk -v OFS="\t" '{print $0,$1":"$2"-"$3":!:-",0,"-"}' >> ${locTMP}${MM}.bed

  LC_COLLATE=C sort -k1,1 -k2,2n --parallel=$THREADS -S60G ${locTMP}${MM}.bed > ${locTMP}${MM}.sort.bed
  bedToBigBed ${locTMP}${MM}.sort.bed ${HUBdir}/${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/uniqueness/${SLURM_ARRAY_TASK_ID}mer.${MM}MM.bb

  #do not continue until file is unblocked by other process
  while [[ -f ${TMPdir}wait.txt ]]; do
    sleep 10s
  done

  #block trackDb from other processes
  touch ${TMPdir}wait.txt

  if [[ $MM -eq 0 ]]; then
    VISIBILITY=pack
  else
    VISIBILITY=hide
  fi
  
  printf "
    track uniquness_${SLURM_ARRAY_TASK_ID}mer_$MM
    shortLabel uniq_${SLURM_ARRAY_TASK_ID}mer_$MM
    longLabel unique regions for ${SLURM_ARRAY_TASK_ID}mers alligned allowing ${MM}MM -best -stratum
    priority $PRIO
    visibility $VISIBILITY
    maxWindowToDraw 10000000
    colorByStrand 115,235,174 140,101,211
    type bigBed 6
    parent uniqueness
    bigDataUrl annotations/uniqueness/${SLURM_ARRAY_TASK_ID}mer.${MM}MM.bb

  " >>${HUBdir}/trackDb.txt

  PRIO=$(( $PRIO + 1 )) 
#unblock file
rm -rf ${TMPdir}wait.txt

done