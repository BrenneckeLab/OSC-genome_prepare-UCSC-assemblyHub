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
locTMP=${TMPdir}map_reads/

#create directories
mkdir $locTMP

#load tools
source ${SCRIPTdir}tools

###################################################################################################

#mapping reads with minimap2
THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))
minimap2 -ax map-ont --secondary=no -t $THREADS $assemblyFASTA $ONT_DNA >${locTMP}reads.sam

#$%#get flagstats of the mapped reads
#$%samtools flagstat -@ $SLURM_CPUS_PER_TASK ${locTMP}reads.sam >${LOG}map_reads/flagstat.txt
#$%
#$%##BAM generation
#$%samtools sort --output-fmt BAM --reference $assemblyFASTA -o ${HUBdir}/annotations/mapped_ONT.all.bam -@ #$%$SLURM_CPUS_PER_TASK ${locTMP}reads.sam
#$%samtools index -@ $SLURM_CPUS_PER_TASK ${HUBdir}/annotations/mapped_ONT.all.bam

###################################################################################################
#filter for primary alignment and create BAM
##filtering
samtools view -q 10 -h -F 2048 ${locTMP}reads.sam |
  samtools view -h -F0x900- >${locTMP}reads.prim.sam

##BAM generation
samtools sort --output-fmt BAM --reference $assemblyFASTA -o ${HUBdir}/annotations//ONT/mapped_ONT.prim.bam -@ $SLURM_CPUS_PER_TASK ${locTMP}reads.prim.sam
samtools index -@ $SLURM_CPUS_PER_TASK ${HUBdir}/annotations/ONT/mapped_ONT.prim.bam

##convert sam to bed
samtools view -@ $SLURM_CPUS_PER_TASK -bS ${locTMP}reads.prim.sam | 
bedtools bamtobed -split -i -  | 
  LC_COLLATE=C sort --parallel $SLURM_CPUS_PER_TASK -S25G -k1,1 -k2,2n |
  #fix bed-entries extending over the end of the chromosome
  awk -v OFS="\t" -v CHROMsize=${HUBdir}/${assemblyNAME}.chrom.sizes '
  BEGIN{
    while((getline I < CHROMsize ) > 0) {
      #split into array by tabs
      split(I,splitLINE,/\t/)
      SIZE[splitLINE[1]]=splitLINE[2]
    }
  }
  {
    if($3>=SIZE[$1]){
      $3=SIZE[$1]-1
    }
    print 
  }'  >${locTMP}reads.prim.sort.fix.bed

###################################################################################################
#try to determine mis-assemblies

#convert bed to bedgraph and find regions with coverage <10
bedtools genomecov -bga -i ${locTMP}reads.prim.sort.fix.bed -g ${HUBdir}/${assemblyNAME}.chrom.sizes | 
  awk -v OFS="\t" -v CHROMsize=${HUBdir}/${assemblyNAME}.chrom.sizes '
  BEGIN{
    while((getline I < CHROMsize ) > 0) {
      split(I,splitLINE,/\t/)
      SIZE[splitLINE[1]]=splitLINE[2]  
    }
  }
  { 
    if( $4<10 && $2>5000 && $3 < SIZE[$1]-5000) {
      if($3-$2>20){
        print $1,$2,$3,"GAP_cov="$4,0,"+"
      }
    }
  }'  | bedtools merge -d 100 > ${locTMP}gaps.bed
 

#convert bed to big-bed
bedToBigBed ${locTMP}gaps.bed ${HUBdir}/${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/gaps.bb

#create gap-track with resected read-ends
#resect ends of read alignments
awk -v OFS="\t" '
{
  $2=$2+20
  $3=$3-20
  print
}' ${locTMP}reads.prim.sort.fix.bed> ${locTMP}reads.prim.sort.fix.resected.bed 

#convert bed to bedgraph and find regions with coverage <10
bedtools genomecov -bga -i ${locTMP}reads.prim.sort.fix.resected.bed -g ${HUBdir}/${assemblyNAME}.chrom.sizes | 
  awk -v OFS="\t" -v CHROMsize=${HUBdir}/${assemblyNAME}.chrom.sizes '
  BEGIN{
    while((getline I < CHROMsize ) > 0) {
      split(I,splitLINE,/\t/)
      SIZE[splitLINE[1]]=splitLINE[2]  
    }
  }
  { 
    if( $4<10 && $2>5000 && $3 < SIZE[$1]-5000) {
      if($3-$2>20){
        print $1,$2,$3,"GAP_cov="$4,0,"+"
      }
    }
  }' | bedtools merge -d 100 > ${locTMP}gaps.resected.bed


#only keep gaps found by resection
bedtools intersect -a ${locTMP}gaps.bed -b ${locTMP}gaps.resected.bed -v > ${locTMP}gaps.resected.filtered.bed

if [[ ! -s ${locTMP}gaps.resected.filtered.bed ]]; then
  awk -v OFS="\t" '
  {
    if(NR==1){
      print $1,1,10,"fake-gap",0,"+" 
    }
  }' <(head ${locTMP}reads.prim.sort.fix.bed ) > ${locTMP}gaps.resected.filtered.bed
fi

#convert bed to big-bed
bedToBigBed ${locTMP}gaps.resected.filtered.bed ${HUBdir}/${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/gaps.resected.bb


###################################################################################################
#add tracks to trackhub

#do not continue until file is unblocked by other process
while [[ -f ${TMPdir}wait.txt ]]; do
  sleep 10s
done

#block trackDb from other processes
touch ${TMPdir}wait.txt

#remove old lines if present
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="mappedONT_DNA|GAPs" '
  {
    if( $0 ~ NAME ) {a=b; } else print
  }
' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp
mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt

#add tracks to trackDb

printf "

track mappedONT_DNA
compositeTrack on
shortLabel ONT_OSCreads_DNA
longLabel mapped OSC Nanopore reads
group Long-Reads
priority 3
visibility squish
type bam
allButtonPair on

  track reads_for_assembly_bam_prim
  type bam
  shortLabel ONTreads_primaryALN
  longLabel read set that was used for assembly
  maxWindowToDraw 10000000
  maxItems 1000000
  bamColorMode strand
  visibility squish
  parent mappedONT_DNA
  bigDataUrl annotations/ONT/mapped_ONT.prim.bam

track GAPs
compositeTrack on
shortLabel GAPs
longLabel Gaps determined by mapping of nanopore data
group Long-Reads
priority 3
visibility pack
type bed
allButtonPair on

  track GAPs_10
  type bigBed
  shortLabel GAPs
  longLabel gaps determined with a cutoff of <10
  maxWindowToDraw 100000000
  maxItems 1000000
  bamColorMode strand
  visibility pack
  parent GAPs
  bigDataUrl annotations/gaps.bb

  track GAPs_resected_10
  type bigBed
  shortLabel resectedGAPs
  longLabel gaps determined with a cutoff of <10 using read-alignments trimmed back by 100nt
  maxWindowToDraw 100000000
  maxItems 1000000
  bamColorMode strand
  visibility pack
  parent GAPs
  bigDataUrl annotations/gaps.resected.bb
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

