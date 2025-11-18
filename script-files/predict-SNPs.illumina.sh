#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --partition=c
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
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
locTMPtop=${TMPdir}SNV/

#create directories
mkdir $locTMPtop

#load tools
source ${SCRIPTdir}tools

###################################################################################################
###################################################################################################
THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))

#map transcripts to contigs
cd ${LOG}SNV/



###################################################################################################
#ONT Clair3
locTMP=${locTMPtop}SNV_Illumina/
mkdir $locTMP

#---------------------------------------------------------------------------------------------------------
#map reads
TEST=""
TEST=$(seqkit fx2tab $ILLUMINA_DNAseq | head -n 1 | grep "/1" )
if [[ -z $TEST ]]; then 
  PAIRED=N
else
  PAIRED=Y
fi

TEST=""
TEST=$(samtools view ${locTMP}reads.mapped.bam | head -n 100 | wc -l)
if [[ $TEST -lt 10 ]]; then

  if [[ $PAIRED == Y && ! -s ${locTMP}read1.fq ]]; then  
    #split paired end data into individual files
    #!remove head
    seqkit seq $ILLUMINA_DNAseq | paste - - - - - - - - \
      | tee >(cut -f 1-4 | tr "\t" "\n" > ${locTMP}read1.fq) \
      | cut -f 5-8 | tr "\t" "\n" > ${locTMP}read2.fq
  fi

  if [[ $COMPUTING == C ]]; then
    cp ${assemblyFASTA} ${locTMP}input.fa
    assemblyFASTA=${locTMP}input.fa
    VARI=${VARI},assemblyFASTA=${locTMP}input.fa

    bwa index ${assemblyFASTA}

    if [[ $PAIRED == Y ]]; then
      COMMAND="${locTMP}read1.fq ${locTMP}read2.fq | samtools view --threads 3 -F 256 -b - | samtools fixmate -m --threads 3  - - | samtools sort -m 5g -T $locTMP --threads 5 - | samtools markdup  --threads 5 - ${locTMP}reads.mapped.bam "
    else
      COMMAND="$ILLUMINA_DNAseq | samtools view --threads 3 -F 256 -b -| samtools sort -m 5g -T $locTMP -o ${locTMP}reads.mapped.bam  --threads 5 - "
    fi

    
    sbatch --wait --job-name=SNV_ILL_bwa -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=20 --mem=60g --wrap="
      SINGULARITYdir=${SINGULARITYdir}
      TMPdir=$TMPdir
      source ${SCRIPTdir}tools
      THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
      bwa mem -v 2 -t \${THREADS} ${assemblyFASTA} $COMMAND  
      
" 
  else
    echo not coded!!!!!
  fi
  samtools index -@ $THREADS ${locTMP}reads.mapped.bam
fi



#predict SNPs using CLAIR3
mkdir ${locTMP}clair/
sbatch --wait --job-name=SNV_ILL_deepvariant --partition=g --gres=gpu:1  -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=15 --mem=100g --wrap="
  set -ux
  THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
  SINGULARITYdir=$SINGULARITYdir
  TMPdir=$TMPdir
  source ${SCRIPTdir}tools
  cd $locTMP
  export TMPDIR=$locTMP

  samtools faidx $assemblyFASTA
  samtools index -@ \$THREADS ${locTMP}reads.mapped.bam

  deepvariant /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=${assemblyFASTA} --reads=${locTMP}reads.mapped.bam --output_vcf=${locTMP}/output.vcf.gz --output_gvcf=${locTMP}/output.g.vcf.gz --intermediate_results_dir ${locTMP}/intermediate_results_dir --num_shards=\$THREADS 


  rtg vcfstats ${locTMP}/output.vcf.gz > ${locTMP}Clair3.vcfstats.txt

  gunzip -c ${locTMP}/output.vcf.gz | grep \"#\" > ${HUBdir}annotations/VARIATIONS/ILL.deepvariant.het.vcf
  cp ${HUBdir}annotations/VARIATIONS/ILL.deepvariant.het.vcf ${HUBdir}annotations/VARIATIONS/ILL.deepvariant.hom.vcf

  gunzip -c ${locTMP}/output.vcf.gz | grep \"0/1\" >> ${HUBdir}annotations/VARIATIONS/ILL.deepvariant.het.vcf
  gunzip -c ${locTMP}/output.vcf.gz | grep \"1/1\" >> ${HUBdir}annotations/VARIATIONS/ILL.deepvariant.hom.vcf

  bgzip -f ${HUBdir}annotations/VARIATIONS/ILL.deepvariant.het.vcf
  bgzip -f ${HUBdir}annotations/VARIATIONS/ILL.deepvariant.hom.vcf


  tabix ${HUBdir}annotations/VARIATIONS/ILL.deepvariant.het.vcf.gz
  tabix ${HUBdir}annotations/VARIATIONS/ILL.deepvariant.hom.vcf.gz
  " 

wait


###################################################################################################
###################################################################################################
#add to trackDb.txt

#do not continue until file is unblocked by other process
while [[ -f ${TMPdir}wait.txt ]]; do
  sleep 10s
done

#block trackDb from other processes
touch ${TMPdir}wait.txt

#remove old lines if present
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="SNV_Illumina_" '
  {
    if( $0 !~ NAME ) print
  }
' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp

mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt

printf "
  track SNV_Illumina_deepvariant_hom
  shortLabel homSNV_deepvariant_ILL
  longLabel homozygous SNVs called using deepvariant and Illumina reads
  priority 30
  visibility squish
  maxWindowToDraw 10000000
  type vcfTabix
  parent VARIATIONS
  bigDataUrl annotations/VARIATIONS/ILL.deepvariant.hom.vcf.gz

  track SNV_Illumina_deepvariant_het
  shortLabel _deepvariant_ILL
  longLabel heterozygous SNVs called using deepvariant and Illumina reads
  priority 33
  visibility squish
  maxWindowToDraw 10000000
  type vcfTabix
  parent VARIATIONS
  bigDataUrl annotations/VARIATIONS/ILL.deepvariant.het.vcf.gz

" >>${HUBdir}/trackDb.txt


#unblock file
rm -rf ${TMPdir}wait.txt
exit
###################################################################################################
#finish script

#clean up
if [[ $DEBUG == N ]]; then
  rm -rf $locTMP
fi

#report processing time
PROCESSED_TIME=$(echo -e $(date "+%s") $TIME | awk '{ print ($1-$2)/60 }')
echo "map transcripts using BLAT=" ${PROCESSED_TIME} >>${LOG}time-log.txt





exit



###################################################################################################
