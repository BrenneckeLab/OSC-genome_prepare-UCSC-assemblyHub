#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --qos=short
#SBATCH --time=0:10:00


hostname
set -u

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"

TIME=$(date "+%s")

###################################################################################################
#setup-phase

#load tools
source ${SCRIPTdir}tools

###################################################################################################
#populate hub.txt

printf "
hub ${assemblyNAME}
shortLabel ${assemblyNAME}
longLabel ${assemblyNAME}
genomesFile genomes.txt
email dominik.handler@imba.oeaw.ac.at
" | sed '/./,$!d' >${OPENdir}hub.txt

###################################################################################################
#populate genomes.txt

defaultPATH=$(sort -k2,2rn ${HUBdir}/${assemblyNAME}.chrom.sizes | head -n 1 | sed 's/>//' | tr ' ' '\t' | cut -f 1)

printf "
genome $assemblyNAME
trackDb ${assemblyNAME}/trackDb.txt
groups groups.txt
description $assemblyNAME
twoBitPath ${assemblyNAME}/${assemblyNAME}.2bit
organism Dmel
defaultPos $defaultPATH:1000-20000
orderKey 4700
#blat brennecke-blat-1.vbc.ac.at 8001
#transBlat brennecke-blat-1.vbc.ac.at 8001
" | sed '/./,$!d' >${OPENdir}genomes.txt

###################################################################################################
#populate groups.txt

printf "
name Annotations
label Annotations
priority 10
defaultIsClosed 0 

name Long-Reads
label Long-Reads
priority 20
defaultIsClosed 1

name Debugging
label Debugging
priority 90
defaultIsClosed 1

" | sed '/./,$!d' >${OPENdir}groups.txt

###################################################################################################
#create cytoband

LC_COLLATE=C sort -k1,1 -k2,2n ${HUBdir}/${assemblyNAME}.chrom.sizes |
  awk '
  {
    n=1
    for(i=0; i<=$2-50000; i=i+50000){
      if( n %2 == 0){ 
        X="gneg"
      }else{
        X="gpos100" 
      }
      print $1,i,i+50000,$1"_"n,X
      n=n+1
    }
    print $1,i,$2,$1"_"n,"gneg"
  }' >${TMPdir}cytoBandIdeo.bed

#cojnvert to bigbed for UCSC
bedToBigBed -type=bed4 ${TMPdir}cytoBandIdeo.bed -as=${UTILITYdir}cyto.as ${HUBdir}/${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/cytoBandIdeo.bigBed

###################################################################################################
#add cytoband to hub

#do not continue until file is unblocked by other process
while [[ -f ${TMPdir}wait.txt ]]; do
  sleep 10s
done

#block trackDb from other processes
touch ${TMPdir}wait.txt

#remove old lines if present
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="cytoBand" '
  {
    if( $0 !~ NAME ) print
  }
' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp

mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt

#add cytoBand to the trackDB file
printf "
track cytoBandIdeo
longLabel Chromosome ideogram with cytogenetic bands @ 50kb intervals
shortLabel cytoBandIdeo50
bigDataUrl annotations/cytoBandIdeo.bigBed
visibility dense
type bigBed
" >>${HUBdir}/trackDb.txt

#unblock file
rm -rf ${TMPdir}wait.txt
