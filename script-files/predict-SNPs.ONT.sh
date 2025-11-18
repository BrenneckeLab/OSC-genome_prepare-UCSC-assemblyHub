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


samtools faidx $assemblyFASTA

###################################################################################################
#ONT Clair3
locTMP=${locTMPtop}CLAIR3/
mkdir $locTMP
mkdir
#---------------------------------------------------------------------------------------------------------
#map reads
if [[ ! -s ${locTMP}reads.mapped.bam ]]; then
  ID=$(sbatch --parsable --wait --job-name=minimap_variations -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=20 --mem=50g --wrap="
      THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
      SINGULARITYdir=$SINGULARITYdir
      TMPdir=$TMPdir
      source ${SCRIPTdir}tools
      minimap2 -Lax map-ont -t \$THREADS --secondary=no $assemblyFASTA $ONT_DNA | samtools sort -@ \$THREADS -m 1G -O BAM - > ${locTMP}reads.mapped.bam 
    
      samtools index -@ \$THREADS ${locTMP}reads.mapped.bam

    " )
    DEPEND="--dependency=afterok:$ID"
else
  DEPEND=""
fi


#predict SNPs using CLAIR3
mkdir ${locTMP}clair/
mkdir ${HUBdir}/annotations/VARIATIONS

sbatch --wait  $DEPEND --job-name=SNV_Clair3 -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=20 --mem=50g --qos=short  --time=2:00:00  --wrap="
  set -ux
  THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
  SINGULARITYdir=$SINGULARITYdir
  TMPdir=$TMPdir  
  source ${SCRIPTdir}tools

  clair3 --bam_fn=${locTMP}reads.mapped.bam --ref_fn=$assemblyFASTA --threads=\${THREADS} --platform=\"ont\" --model_path=\"/opt/models/${MODEL_DNA}\" --output=${locTMP}clair/ --include_all_ctgs --enable_long_indel

  rtg vcfstats ${locTMP}/clair//merge_output.vcf.gz > ${locTMP}Clair3.vcfstats.txt

  gunzip -c ${locTMP}/clair//merge_output.vcf.gz | grep \"#\" > ${HUBdir}/annotations/VARIATIONS/ONT.CLAIR3.het.vcf
  cp ${HUBdir}/annotations/VARIATIONS/ONT.CLAIR3.het.vcf ${HUBdir}/annotations/VARIATIONS/ONT.CLAIR3.hom.vcf

  gunzip -c ${locTMP}/clair//merge_output.vcf.gz | grep \"0/1\" >> ${HUBdir}/annotations/VARIATIONS/ONT.CLAIR3.het.vcf
  gunzip -c ${locTMP}/clair//merge_output.vcf.gz | grep \"1/1\" >> ${HUBdir}/annotations/VARIATIONS/ONT.CLAIR3.hom.vcf

  bgzip -f ${HUBdir}/annotations/VARIATIONS/ONT.CLAIR3.het.vcf
  bgzip -f ${HUBdir}/annotations/VARIATIONS/ONT.CLAIR3.hom.vcf

  tabix ${HUBdir}/annotations/VARIATIONS/ONT.CLAIR3.het.vcf.gz
  tabix ${HUBdir}/annotations/VARIATIONS/ONT.CLAIR3.hom.vcf.gz
  " &


sbatch --wait $DEPEND --job-name=SNV_pepper-deepvariant -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=20 --mem=100g --partition=m --qos=short --time=2:00:00 --wrap="
  set -ux
  THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
  SINGULARITYdir=$SINGULARITYdir
  TMPdir=$TMPdir
  source ${SCRIPTdir}tools
  
  ulimit -u 100000
  pepper_deepvariant call_variant -b ${locTMP}reads.mapped.bam  -f $assemblyFASTA -o ${locTMP}pepper-deepvariant -p pepper-deepvariant -t \${THREADS} --ont_r9_guppy5_sup
  
  rtg vcfstats ${locTMP}pepper-deepvariant/pepper-deepvariant.vcf.gz > ${locTMP}pepper-deepvariant.vcfstats.txt

  gunzip -c ${locTMP}pepper-deepvariant/pepper-deepvariant.vcf.gz | grep \"#\" >${HUBdir}/annotations/VARIATIONS/ONT.pepper-deepvariant.hom.vcf
  cp ${HUBdir}/annotations/VARIATIONS/ONT.pepper-deepvariant.hom.vcf ${HUBdir}/annotations/VARIATIONS/ONT.pepper-deepvariant.het.vcf

  gunzip -c ${locTMP}pepper-deepvariant/pepper-deepvariant.vcf.gz | grep \"0/1\" >> ${HUBdir}/annotations/VARIATIONS/ONT.pepper-deepvariant.het.vcf
  gunzip -c ${locTMP}pepper-deepvariant/pepper-deepvariant.vcf.gz | grep \"1/1\" >> ${HUBdir}/annotations/VARIATIONS/ONT.pepper-deepvariant.hom.vcf

  bgzip -f ${HUBdir}/annotations/VARIATIONS/ONT.pepper-deepvariant.het.vcf
  bgzip -f ${HUBdir}/annotations/VARIATIONS/ONT.pepper-deepvariant.hom.vcf

  mkdir ${HUBdir}/annotations/VARIATIONS
  tabix ${HUBdir}/annotations/VARIATIONS/ONT.pepper-deepvariant.het.vcf.gz
  tabix ${HUBdir}/annotations/VARIATIONS/ONT.pepper-deepvariant.hom.vcf.gz
  " &

sbatch --wait $DEPEND --job-name=SV_sniffles -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=20 --mem=50g --qos=short --time=2:00:00 --wrap="
  set -ux
  THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
  SINGULARITYdir=$SINGULARITYdir
  TMPdir=$TMPdir
  source ${SCRIPTdir}tools
  
  samtools calmd -b -@ \$THREADS ${locTMP}reads.mapped.bam $assemblyFASTA > ${locTMP}reads.mapped.md.bam
  sniffles -m ${locTMP}reads.mapped.md.bam -v ${locTMP}SV.sniffles.vcf

  bcftools sort ${locTMP}SV.sniffles.vcf > ${locTMP}SV.sniffles.sort.vcf
  cat ${locTMP}SV.sniffles.sort.vcf | grep \"#\" >${HUBdir}/annotations/VARIATIONS/ONT.sniffles.hom.vcf
  cp ${HUBdir}/annotations/VARIATIONS/ONT.sniffles.hom.vcf ${HUBdir}/annotations/VARIATIONS/ONT.sniffles.het.vcf

  cat ${locTMP}SV.sniffles.sort.vcf  | grep \"0/1\" >> ${HUBdir}/annotations/VARIATIONS/ONT.sniffles.het.vcf
  cat ${locTMP}SV.sniffles.sort.vcf  | grep \"1/1\" >> ${HUBdir}/annotations/VARIATIONS/ONT.sniffles.hom.vcf

  bgzip -f ${HUBdir}/annotations/VARIATIONS/ONT.sniffles.het.vcf
  bgzip -f ${HUBdir}/annotations/VARIATIONS/ONT.sniffles.hom.vcf

  tabix ${HUBdir}/annotations/VARIATIONS/ONT.sniffles.het.vcf.gz
  tabix ${HUBdir}/annotations/VARIATIONS/ONT.sniffles.hom.vcf.gz
  " &


sbatch --wait $DEPEND --job-name=SV_cuteSV -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=20 --mem=50g --qos=short --time=2:00:00 --wrap="
  set -ux
  THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
  SINGULARITYdir=$SINGULARITYdir
  TMPdir=$TMPdir
  source ${SCRIPTdir}tools
  
  cuteSV ${locTMP}reads.mapped.bam $assemblyFASTA ${locTMP}cuteSV.vcf $locTMP --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --genotype --threads \$THREADS
  
  cat ${locTMP}cuteSV.vcf | grep \"#\" >${HUBdir}/annotations/VARIATIONS/ONT.cuteSV.hom.vcf
  cp ${HUBdir}/annotations/VARIATIONS/ONT.cuteSV.hom.vcf ${HUBdir}/annotations/VARIATIONS/ONT.cuteSV.het.vcf

  cat ${locTMP}cuteSV.vcf  | grep \"0/1\" >> ${HUBdir}/annotations/VARIATIONS/ONT.cuteSV.het.vcf
  cat ${locTMP}cuteSV.vcf  | grep \"1/1\" >> ${HUBdir}/annotations/VARIATIONS/ONT.cuteSV.hom.vcf

  bgzip -f ${HUBdir}/annotations/VARIATIONS/ONT.cuteSV.het.vcf
  bgzip -f ${HUBdir}/annotations/VARIATIONS/ONT.cuteSV.hom.vcf

  
  tabix ${HUBdir}/annotations/VARIATIONS/ONT.cuteSV.het.vcf.gz
  tabix ${HUBdir}/annotations/VARIATIONS/ONT.cuteSV.hom.vcf.gz
  " &

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
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="SNV_ONT|SV_ONT" '
  {
    if( $0 !~ NAME ) print
  }
' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp

mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt

printf "
  track SNV_ONT_CLAIR3_hom
  shortLabel homSNV_CLAIR3_ONT
  longLabel homozygous SNVs called using CLAIR3 and ONT reads
  priority 31
  visibility squish
  maxWindowToDraw 10000000
  type vcfTabix
  parent VARIATIONS
  bigDataUrl annotations/VARIATIONS/ONT.CLAIR3.hom.vcf.gz

  track SNV_ONT_CLAIR3_het
  shortLabel hetSNV_CLAIR3_ONT
  longLabel heterozygous SNVs called using CLAIR3 and ONT reads
  priority 34
  visibility squish
  maxWindowToDraw 10000000
  type vcfTabix
  parent VARIATIONS
  bigDataUrl annotations/VARIATIONS/ONT.CLAIR3.het.vcf.gz

  track SNV_ONT_pepper-deepvariant_hom
  shortLabel homSNV_pepper-deepvariant_ONT
  longLabel homozygous SNVs called using pepper-deepvariant and ONT reads
  priority 32
  visibility squish
  maxWindowToDraw 10000000
  type vcfTabix
  parent VARIATIONS
  bigDataUrl annotations/VARIATIONS/ONT.pepper-deepvariant.hom.vcf.gz

  track SNV_ONT_pepper-deepvariant_het
  shortLabel hetSNV_pepper-deepvariant_ONT
  longLabel heterozygous SNVs called using pepper-deepvariant and ONT reads
  priority 35
  visibility squish
  maxWindowToDraw 10000000
  type vcfTabix
  parent VARIATIONS
  bigDataUrl annotations/VARIATIONS/ONT.pepper-deepvariant.het.vcf.gz

  track SV_ONT_cuteSV_hom
  shortLabel homSV_cuteSV_ONT
  longLabel homozygous SVs called using cuteSV and ONT reads
  priority 36
  visibility squish
  maxWindowToDraw 10000000
  type vcfTabix
  parent VARIATIONS
  bigDataUrl annotations/VARIATIONS/ONT.cuteSV.hom.vcf.gz

  track SV_ONT_cuteSV_het
  shortLabel hetSV_cuteSV_ONT
  longLabel heterozygous SVs called using cuteSV and ONT reads
  priority 38
  visibility squish
  maxWindowToDraw 10000000
  type vcfTabix
  parent VARIATIONS
  bigDataUrl annotations/VARIATIONS/ONT.cuteSV.het.vcf.gz

  track SV_ONT_sniffles_hom
  shortLabel homSV_sniffles_ONT
  longLabel homozygous SVs called using sniffles and ONT reads
  priority 37
  visibility squish
  maxWindowToDraw 10000000
  type vcfTabix
  parent VARIATIONS
  bigDataUrl annotations/VARIATIONS/ONT.sniffles.hom.vcf.gz

  track SV_ONT_sniffles_het
  shortLabel hetSV_sniffles_ONT
  longLabel heterozygous SVs called using sniffles and ONT reads
  priority 39
  visibility squish
  maxWindowToDraw 10000000
  type vcfTabix
  parent VARIATIONS
  bigDataUrl annotations/VARIATIONS/ONT.sniffles.het.vcf.gz
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

