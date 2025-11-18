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
locTMP=${TMPdir}de-novo-annotations/

#create directories
mkdir $locTMP

#load tools
source ${SCRIPTdir}tools

###################################################################################################
###################################################################################################
THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))

#map transcripts to contigs
cd ${LOG}create-gene-tracks/de-novo-annotation/

ID=""

###################################################################################################
#prepare input reads

#map Illumina reads using hisat
if [[ ! -s ${locTMP}mapped.illumina.sort.bam.bai || $FORCE == Y ]]; then
  #prepare log directory
  rm -rf ${LOG}create-gene-tracks/de-novo-annotation/hisat*


  ILLUMINA_RNAseq=$(echo $ILLUMINA_RNAseq | tr '~' ' '  )

  newID=$(sbatch --parsable --job-name=hisat -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=20 --mem=50g --wrap="
    set -ux 
    THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
    SINGULARITYdir=$SINGULARITYdir
    TMPdir=$TMPdir
    source ${SCRIPTdir}tools

    seqkit seq ${ILLUMINA_RNAseq} > ${locTMP}input.illumina.fq

    hisat2_build $assemblyFASTA -p \$THREADS ${locTMP}hisat.index
    hisat2 -x ${locTMP}hisat.index -p \$THREADS -U ${locTMP}input.illumina.fq  | samtools view -bS > ${locTMP}mapped.illumina.bam


    samtools sort ${locTMP}mapped.illumina.bam -m1G -@ \$THREADS -T ${locTMP} > ${locTMP}mapped.illumina.sort.bam
    rm -rf ${locTMP}mapped.illumina.sort.bam.bai
    samtools index ${locTMP}mapped.illumina.sort.bam -@ \$THREADS

    bedtools genomecov -bg -split -strand + -g ${HUBdir}/${assemblyNAME}.chrom.sizes -ibam ${locTMP}mapped.illumina.sort.bam |
      LC_COLLATE=C sort -k1,1 -k2,2n --parallel \$THREADS -S1G> ${locTMP}mapped.illumina.sort.sense.bg
    bedtools genomecov -bg -split -strand - -g ${HUBdir}/${assemblyNAME}.chrom.sizes -ibam ${locTMP}mapped.illumina.sort.bam | awk -v OFS=\"\\t\" '{\$4=\$4*-1; print \$0}' |
      LC_COLLATE=C sort -k1,1 -k2,2n --parallel \$THREADS -S1G> ${locTMP}mapped.illumina.sort.antisense.bg

    bedGraphToBigWig ${locTMP}mapped.illumina.sort.sense.bg ${HUBdir}/${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/de-novo-transcripts.reads.Illumina.sense.bw
    bedGraphToBigWig ${locTMP}mapped.illumina.sort.antisense.bg ${HUBdir}/${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/de-novo-transcripts.reads.Illumina.antisense.bw
  "
  )


    ID="${ID}:${newID}"
fi

if [[ ! -s ${locTMP}mapped.ONT.sort.bam.bai || $FORCE == Y ]]; then
  #prepare log directory
  rm -rf ${LOG}create-gene-tracks/de-novo-annotation/minimap*

  #map transcripts to contigs
  cd ${LOG}create-gene-tracks/de-novo-annotation/

  ONT_deNOVO=${ONT_deNOVO//\~/*.fq.gz }
  ONT_deNOVO=${ONT_deNOVO}\*.fq.gz

  newID=$(sbatch --parsable --job-name=minimap -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=20 --mem=50g --wrap="
    set -ux 
    THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
    SINGULARITYdir=$SINGULARITYdir
    TMPdir=$TMPdir
    source ${SCRIPTdir}tools

    minimap2 -t \$THREADS -ax splice $assemblyFASTA $ONT_deNOVO | samtools view -bS > ${locTMP}mapped.ONT.bam
    
    samtools sort ${locTMP}mapped.ONT.bam -m1G -@ \$THREADS -T ${locTMP} > ${locTMP}mapped.ONT.sort.bam
    rm -rf ${locTMP}mapped.ONT.sort.bam.bai
    samtools index ${locTMP}mapped.ONT.sort.bam -@ \$THREADS
  "

  )

    ID="${ID}:${newID}"
fi


###################################################################################################
#assemble transcripts
#prepare log directory

if [[ ! -s ${locTMP}de-novo-mix.without.gtf ]]; then
  rm -rf ${LOG}create-gene-tracks/de-novo-annotation//stringtie*
  cd ${LOG}create-gene-tracks/de-novo-annotation/

  if [[ -z $ID ]]; then DEPEND=""; else DEPEND="--dependency=afterok$ID"; fi

  sbatch --wait $DEPEND --job-name=stringtie -o "%x.o.%A-%a.txt" -e "%x.e.%A-%a.txt" --cpus-per-task=10 --mem=80g --wrap="
    set -ux 
    THREADS=\$(( \$SLURM_CPUS_PER_TASK * 2 ))
    SINGULARITYdir=$SINGULARITYdir
    TMPdir=$TMPdir
    source ${SCRIPTdir}tools

    #? stringtie --mix -p \$THREADS -v -a 4 -o ${locTMP}de-novo-mix.gtf -G ${locTMP}lifted.gtf ${locTMP}mapped.illumina.sort.bam ${locTMP}mapped.ONT.sort.bam 
    stringtie --mix -p \$THREADS -v -a 4 -j 2 -c 3.9 -s 3.9  -o ${locTMP}de-novo-mix.without.gtf ${locTMP}mapped.illumina.sort.bam ${locTMP}mapped.ONT.sort.bam 

  "
fi

###################################################################################################
#prepare genepred
LC_COLLATE=C sort -k1,1 -k4,4n ${locTMP}de-novo-mix.without.gtf > ${locTMP}de-novo-mix.without.sort.gtf
gtfToGenePred  ${locTMP}de-novo-mix.without.sort.gtf ${locTMP}de-novo-mix.without.gp
genePredToBigGenePred ${locTMP}de-novo-mix.without.gp stdout | LC_COLLATE=C sort -k1,1 -k2,2n > ${locTMP}de-novo-mix.without.bigGenePred

#?@ LC_COLLATE=C sort -k1,1 -k4,4n ${locTMP}de-novo-mix.gtf > ${locTMP}de-novo-mix.sort.gtf
#?@ gtfToGenePred  ${locTMP}de-novo-mix.sort.gtf ${locTMP}de-novo-mix.gp
#?@ genePredToBigGenePred ${locTMP}de-novo-mix.gp stdout | LC_COLLATE=C sort -k1,1 -k2,2n > ${locTMP}de-novo-mix.bigGenePred

printf 'table bigGenePred
"bigGenePred gene models"
    (
    string chrom;       	"Reference sequence chromosome or scaffold"
    uint   chromStart;  	"Start position in chromosome" 
    uint   chromEnd;    	"End position in chromosome"
    string name;        	"Name or ID of item, ideally both human-readable and unique"
    uint score;         	"Score (0-1000)"
    char[1] strand;     	"+ or - for strand"
    uint thickStart;    	"Start of where display should be thick (start codon)"
    uint thickEnd;      	"End of where display should be thick (stop codon)"
    uint reserved;       	"RGB value (use R,G,B string in input file)"
    int blockCount;     	"Number of blocks"
    int[blockCount] blockSizes; "Comma separated list of block sizes"
    int[blockCount] chromStarts;"Start positions relative to chromStart"
    string name2;       	"Alternative/human readable name"
    string cdsStartStat; 	"Status of CDS start annotation (none, unknown, incomplete, or complete)"
    string cdsEndStat;   	"Status of CDS end annotation (none, unknown, incomplete, or complete)"
    int[blockCount] exonFrames; "Exon frame {0,1,2}, or -1 if no frame for exon"
    string type;        	"Transcript type"
    string geneName;    	"Primary identifier for gene"
    string geneName2;   	"Alternative/human-readable gene name"
    string geneType;    	"Gene type"
    )  
' > ${locTMP}bigGenePred.as

bedToBigBed -tab -as=${locTMP}bigGenePred.as -type=bed12+8 ${locTMP}de-novo-mix.without.bigGenePred ${HUBdir}/${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/de-novo-transcripts.bb
#?bedToBigBed -tab -as=${locTMP}bigGenePred.as -type=bed12+8 ${locTMP}de-novo-mix.bigGenePred ${HUBdir}/${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/de-novo-transcripts.guided.bb

###################################################################################################
#add to trackDb.txt
cp ${locTMP}mapped.ONT.sort.bam ${HUBdir}/annotations/de-novo-transcripts.reads.ONT.bam
cp ${locTMP}mapped.ONT.sort.bam.bai ${HUBdir}/annotations/de-novo-transcripts.reads.ONT.bam.bai
cp ${locTMP}mapped.illumina.sort.bam ${HUBdir}/annotations/de-novo-transcripts.reads.Illumina.bam
cp ${locTMP}mapped.illumina.sort.bam.bai ${HUBdir}/annotations/de-novo-transcripts.reads.Illumina.bam.bai

#do not continue until file is unblocked by other process
while [[ -f ${TMPdir}wait.txt ]]; do
  sleep 10s
done

#block trackDb from other processes
touch ${TMPdir}wait.txt

#remove old lines if present
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="de-novo-transcripts" '
  {
    if( $0 !~ NAME ) print
  }
' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp

mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt

printf "
track de-novo-transcripts
shortLabel de-novo-transcripts
longLabel transcripts assembled from Illumna and ONT reads
group Debugging
visibility hide
maxWindowToDraw 10000000
colorByStrand 50,50,150 150,50,50
type bigGenePred 
bigDataUrl annotations/de-novo-transcripts.bb

track de-novo-transcripts-ONT
shortLabel ONT-de-novo
longLabel Nanopore reads used for de-novo transcript prediction
group Debugging
visibility hide
maxWindowToDraw 10000000
maxItems 1000000
bamColorMode strand
type bam 
bigDataUrl annotations/de-novo-transcripts.reads.ONT.bam

track de-novo-transcripts-Illumina-coverage
container multiWig
aggregate transparentOverlay
group Debugging
showSubtrackColorOnUi on
shortLabel ILL-de-novo
boxedCfg on
autoScale on
alwaysZero on
longLabel Illumina reads used for de-novo transcript prediction
type bigWig
visibility hide
maxHeightPixels 100:50:8

  track de-novo-transcripts-Illumina-coverage_+
  type bigWig
  bigDataUrl annotations/de-novo-transcripts.reads.Illumina.sense.bw
  shortLabel ILL-de-novo_+
  longLabel Illumina reads used for de-novo transcript prediction_+
  parent de-novo-transcripts-Illumina-coverage
  windowingFunction mean

  track de-novo-transcripts-Illumina-coverage_-
  type bigWig
  bigDataUrl annotations/de-novo-transcripts.reads.Illumina.antisense.bw
  shortLabel ILL-de-novo_-
  longLabel Illumina reads used for de-novo transcript prediction_-
  parent de-novo-transcripts-Illumina-coverage
  windowingFunction mean

" >>${HUBdir}/trackDb.txt

#? track de-novo-transcripts-guided
#? shortLabel de-novo-transcripts-guided
#? longLabel transcripts assembled from Illumna and ONT reads guided by lifted annotations
#? group Debugging
#? visibility hide
#? maxWindowToDraw 10000000
#? colorByStrand 50,50,150 150,50,50
#? type bigGenePred 
#? bigDataUrl annotations/de-novo-transcripts.guided.bb

#unblock file
rm -rf ${TMPdir}wait.txt

###################################################################################################


bedtools bamtobed -i ${locTMP}mapped.ONT.sort.bam | mawk -v OFS="\t" '{
  if($6=="-"){
    $3=$2+1
  }else{
    $2=$3-1
  }
  print
}' |
  LC_COLLATE=C sort --parallel=$THREADS -S10G -k1,1 -k2,2n > ${locTMP}3ends.sort.bed
bedtools genomecov -strand + -bg -i ${locTMP}3ends.sort.bed -g ${HUBdir}/${assemblyNAME}.chrom.sizes > ${locTMP}3end.sense.bg
bedtools genomecov -strand - -bg -i ${locTMP}3ends.sort.bed -g ${HUBdir}/${assemblyNAME}.chrom.sizes | awk -v OFS="\t" '{$4=-$4; print }'> ${locTMP}3end.antisense.bg

bedGraphToBigWig ${locTMP}3end.sense.bg  ${HUBdir}/${assemblyNAME}.chrom.sizes  ${HUBdir}annotations/3ends.sense.bw
bedGraphToBigWig ${locTMP}3end.antisense.bg ${HUBdir}/${assemblyNAME}.chrom.sizes  ${HUBdir}/annotations/3ends.antisense.bw

bedtools bamtobed -i ${locTMP}mapped.ONT.sort.bam | mawk -v OFS="\t" -v X=5 '{
  if($6=="-"){
    $3=$2+X

    if($2-X>=0){
      $2=$3-2
    }else{
      $2=0
    }
  }else{
    if($3-X>=0){
      $2=$3-X 
    }else{
      $2=0
    }
    $3=$2+X
  } 
  print  
}' |
  LC_COLLATE=C sort --parallel=$THREADS -S10G -k1,1 -k2,2n > ${locTMP}3ends.sort.10.bed
bedtools genomecov -strand + -bg -i ${locTMP}3ends.sort.10.bed -g ${HUBdir}/${assemblyNAME}.chrom.sizes > ${locTMP}3end.sense.10.bg
bedtools genomecov -strand - -bg -i ${locTMP}3ends.sort.10.bed -g ${HUBdir}/${assemblyNAME}.chrom.sizes | awk -v OFS="\t" '{$4=-$4; print }'> ${locTMP}3end.antisense.10.bg

bedGraphToBigWig ${locTMP}3end.sense.10.bg  ${HUBdir}/${assemblyNAME}.chrom.sizes  ${HUBdir}annotations/3ends.sense.10.bw
bedGraphToBigWig ${locTMP}3end.antisense.10.bg ${HUBdir}/${assemblyNAME}.chrom.sizes  ${HUBdir}/annotations/3ends.antisense.10.bw

#do not continue until file is unblocked by other process
while [[ -f ${TMPdir}wait.txt ]]; do
  sleep 10s
done

#block trackDb from other processes
touch ${TMPdir}wait.txt

#remove old lines if present
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="curated3end" '
  {
    if( $0 !~ NAME ) print
  }
' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp

mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt

printf "
track curated3end-single
container multiWig
aggregate transparentOverlay
group Debugging
showSubtrackColorOnUi on
shortLabel 3ends_single-nt
boxedCfg on
autoScale on
alwaysZero on
longLabel 3 prime ends of nanopore direct RNA reads - single nucleotide resolution
type bigWig
visibility hide
maxHeightPixels 100:100:8

  track curated3end-single_+
  type bigWig
  color 50,50,150
  bigDataUrl annotations/3ends.sense.bw
  shortLabel 3ends_single-nt_+
  longLabel 3 prime ends of nanopore direct RNA reads - single nucleotide resolution sense
  parent curated3end-single
  windowingFunction mean

  track curated3end-single_-
  type bigWig
  color 150,50,50
  bigDataUrl annotations/3ends.antisense.bw
  shortLabel 3ends_single-nt_-
  longLabel 3 prime ends of nanopore direct RNA reads - single nucleotide resolution antisense
  parent curated3end-single
  windowingFunction mean

track curated3end-10nt
container multiWig
aggregate transparentOverlay
group Debugging
showSubtrackColorOnUi on
shortLabel 3ends_10nt
boxedCfg on
autoScale on
alwaysZero on
longLabel 3 prime ends of nanopore direct RNA reads - +/- 5 nucleotides
type bigWig
visibility hide
maxHeightPixels 100:100:8

  track curated3end-10nt_+
  type bigWig
  color 50,50,150
  bigDataUrl annotations/3ends.sense.10.bw
  shortLabel 3ends_10nt_+
  longLabel 3 prime ends of nanopore direct RNA reads - +/- 5 nucleotides sense
  parent curated3end-10nt
  windowingFunction mean

  track curated3end-10nt_-
  type bigWig
  color 150,50,50
  bigDataUrl annotations/3ends.antisense.10.bw
  shortLabel 3ends_10nt_-
  longLabel 3 prime ends of nanopore direct RNA reads - +/- 5 nucleotides antisense
  parent curated3end-10nt
  windowingFunction mean
" >>${HUBdir}/trackDb.txt

#unblock file
rm -rf ${TMPdir}wait.txt


###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

###################################################################################################
###################################################################################################
###################################################################################################

###################################################################################################
#finish script

#clean up
if [[ $DEBUG == N ]]; then
  rm -rf $locTMP
fi

#report processing time
PROCESSED_TIME=$(echo -e $(date "+%s") $TIME | awk '{ print ($1-$2)/60 }')
echo "map transcripts using BLAT=" ${PROCESSED_TIME} >>${LOG}time-log.txt

