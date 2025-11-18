#!/bin/bash

#SBATCH --cpus-per-task=15
#SBATCH --mem=20g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --qos=short
#SBATCH --time=5:00:00


hostname
set -u

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"

TIME=$(date "+%s")

###################################################################################################
#setup-phase

#create path variables
topOPENdir=$OPENdir
locTMP=${TMPdir}map_transcripts/

#create directories
mkdir $locTMP

#load tools
source ${SCRIPTdir}tools

###################################################################################################
###################################################################################################
#define variables used for Y exclusion
if [[ $Yinclude == Y ]]; then
  PATTERN="±"
  COMMAND='cat'
else
  PATTERN="loc=Y"
  COMMAND='grep -v loc=Y'
fi

THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))
###################################################################################################
#process full-transcripts

# #download transcript-fasta file and pre-process
wget --no-verbose -O ${locTMP}transcripts.fa.gz http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.44_FB2022_01/fasta/dmel-all-transcript-r6.44.fasta.gz

#change the sequence name to the trancript name and filter Y located transcripts if required
gunzip -c ${locTMP}transcripts.fa.gz |
  eval $COMMAND |
  awk -v OFS="\t" -v RS=">" -v locTMP=$locTMP '
  {
    sub("\n", "\t")
    gsub("\n", "")
    if(NR>1){
      for(i=1; i<=NF; i++){
        if($i~"name="){
          split($i,splitNAME,/=|;/)  
        }
      }

      print $1, splitNAME[2] > locTMP "fbtr_to_name.txt"
      print ">"splitNAME[2]"\n"$NF
    }
  }' >${locTMP}transcripts_filtered.fa

#* these parts are only required if minimap is used for mapping of transcripts.
#* #map using minimap2
#* ${SINGULARITYdir}minimap2.simg minimap2 -ax splice -uf -C5 $assemblyFASTA ${locTMP}transcripts_filtered.fa > ${locTMP}mapped.sam

#split transcripts into multiple files for parallel processing
seqkit split2 --line-width 0 --force  -p ${THREADS} ${locTMP}transcripts_filtered.fa

#create input IDs for parallel
IDs=$(seq -w 1 100 | head -n $THREADS)

#run blat-mapping in parallel
parallel --linebuffer -j $THREADS "
  blat $assemblyFASTA ${locTMP}transcripts_filtered.fa.split/transcripts_filtered.part_{}.fa -stepSize=5 -fine -noHead -q=dna ${locTMP}transcripts_{}.psl
" ::: $IDs

#filter blat mappings based on fraction of transcript aligned
awk -v OFS="\t" '{
    if($1 > ($11 * 0.90)) {
      print
    }
  }' <(cat ${locTMP}transcripts_0*.psl) >${locTMP}mapped_blat.psl

#convert BLAT output to bed-file
pslToBed ${locTMP}mapped_blat.psl ${locTMP}mapped.bed

#* #extract unmapped transcripts
#* #only works with minimap mappings but not BLAT
#* samtools view -f 4 ${locTMP}mapped.sam | samtools fasta -0 ${locTMP}unmapped_transcripts.fa -
#*  mv ${locTMP}unmapped_transcripts.fa ${LOG}map_transcripts/unmapped_transcripts.fa

#* #convert mappings to bed12
#* samtools view -bS ${locTMP}mapped.sam |
#*   ${SINGULARITYdir}bedtools.simg bamtobed -bed12 -i - > ${locTMP}mapped.bed

#---------------------------------------------------------------------------------------------------------
#process UTRs

#create star index to map short UTRs later
mkdir ${locTMP}STAR
if [[ ! -s ${locTMP}STAR/SAindex ]]; then
  star --runThreadN $THREADS --runMode genomeGenerate --genomeSAindexNbases 12 --genomeDir ${locTMP}STAR --genomeFastaFiles ${assemblyFASTA}
fi

#five_prime_UTR three_prime_UTR 
for EXT in CDS; do

  mkdir -p ${locTMP}${EXT}_split/
  
  #download sequences
  wget no-verbose -O ${locTMP}${EXT}.fa.gz http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.44_FB2022_01/fasta/dmel-all-${EXT}-r6.44.fasta.gz

  #convert name from FBTR to transcript name
  if [[ $EXT == "CDS" ]]; then
    seqkit seq --only-id --line-width 0  ${locTMP}${EXT}.fa.gz >${locTMP}${EXT}.fa
  else
    gunzip -c ${locTMP}${EXT}.fa.gz |
      eval $COMMAND |
      awk -v OFS="\t" -v fbtr_to_name=${locTMP}fbtr_to_name.txt '
        BEGIN{
          while((getline LINE < fbtr_to_name ) > 0) {
            #split into array by tabs
            split(LINE,splitLINE,/\t| /);
            FBTRtoNAME[splitLINE[1]]=splitLINE[2]
          }
          RS=">"
        }
        {
          sub("\n", "\t")
          gsub("\n", "")
          if(NR>1){
            sub(">","",$1)
            print ">"FBTRtoNAME[$1]"\n"$NF
          }
        }' >${locTMP}${EXT}.fa
  fi
  #separate sequences <50 to map them with star
  ##star has a maximum for input sequence length; also 250 nt long sequences should get mapped by blat
  seqkit fx2tab ${locTMP}${EXT}.fa |
    awk -v locTMP=${locTMP} -v EXT=$EXT '
      {
        if(length($2)< 50){
          print ">"$1"\n"$NF > locTMP EXT"_short.fa"
        }else{
          print ">"$1"\n"$NF 
        }
      }' |
    seqkit split2 --force --out-dir ${locTMP}${EXT}_split/ -p ${THREADS}

  #map and convert to bed
  #create input IDs for parallel
  IDs=$(seq -w 1 100 | head -n $THREADS)

  #run blat-mapping in parallel
  parallel --linebuffer -j $THREADS "
  blat $assemblyFASTA ${locTMP}${EXT}_split/stdin.part_{}.fasta -stepSize=5 -fine -noHead -q=dna ${locTMP}${EXT}_{}.psl" ::: $IDs

  #merge mapped UTRs and convert to bed
  awk -v OFS="\t" '{
    if($1 > ($11 * 0.90)) {
      print
    }
  }' <(cat ${locTMP}${EXT}_*.psl) >${locTMP}${EXT}.psl

  pslToBed ${locTMP}${EXT}.psl ${locTMP}mapped_${EXT}.bed

  #map short sequences using star
  mkdir ${locTMP}star_${EXT}/
  star --runThreadN $THREADS --genomeDir ${locTMP}STAR --readFilesIn ${locTMP}${EXT}_short.fa --outFileNamePrefix ${locTMP}star_${EXT}/ --outSAMmode NoQS --readFilesCommand cat --alignEndsType Local --twopassMode Basic --outReadsUnmapped Fastx --outMultimapperOrder Random --outSAMtype SAM --outFilterMultimapNmax 1000 --winAnchorMultimapNmax 2000 --alignSoftClipAtReferenceEnds No --outFilterMatchNmin 15

  # #move all unmapped reads to the log-folder
  mv ${locTMP}star_${EXT}/Unmapped.out.mate1 ${LOG}map_transcripts/unmapped_STAR_${EXT}.fa
  # mv ${locTMP}unmapped_long_${EXT}.fa ${LOG}map_transcripts/unmapped_long_${EXT}.fa

  #convert star mapped reads to bed and add them to the bed file from the minimap mapping
  samtools view -bS ${locTMP}star_${EXT}/Aligned.out.sam |
    bedtools bamtobed -bed12 -i - >>${locTMP}mapped_${EXT}.bed

  LC_COLLATE=C sort --parallel=$THREADS -k1,1 -k2,2n ${locTMP}mapped_${EXT}.bed >${locTMP}mapped_${EXT}_sort.bed

  #?#don't think this is required
  #?generate bigbed file
  #?bedToBigBed -type=bed12 ${locTMP}mapped_${EXT}_sort.bed ${OPENdir}/${assemblyNAME}.chrom.sizes ${OPENdir}/${EXT}.bb

done
wait


#---------------------------------------------------------------------------------------------------------
#remove short introns as they are most likely wrog due to sequence divergence between transcripts and the assembly

bedtools bedtobam -bed12 -i ${locTMP}mapped.bed -g ${HUBdir}/${assemblyNAME}.chrom.sizes | samtools view -h  | 
  awk -v OFS="\t" '
  {
    if($0~"^@"){
      print
    }else{
      TAG=$6
      gsub(/[A-Z]/," & ",TAG)
      n=split(TAG,splitTAG,/ |\t/)

      newTAG=""
      for(i=2; i<=n; i=i+2){
        newTAG=newTAG splitTAG[i-1] 
        if(splitTAG[i]=="N" && splitTAG[i-1]<20){
          
          newTAG=newTAG "M"
        }else{
          newTAG=newTAG splitTAG[i]
        }
      }
      $6=newTAG
      print
    }
  }' | samtools view -bS | bedtools bamtobed -bed12 -i - >${locTMP}mapped.fused.bed

#---------------------------------------------------------------------------------------------------------
#move thick bounds according to the UTR coordinates
awk -v OFS="\t" -v locTMP=${locTMP} '
BEGIN{
  #define psl column names
    matches=1
    misMatches=2
    repMatches=3
    nCount=4
    qNumInsert=5
    qBaseInsert=6
    tNumInsert=7
    tBaseInsert=8
    strand=9
    qName=10
    qSize=11
    qStart=12
    qEnd=13
    tName=14
    tSize=15
    tStart=16
    tEnd=17
    blockCount=18
    blockSizes=19
    qStarts=20
    tStarts=21

  FILE=locTMP "five_prime_UTR.psl"
  while((getline LINE < FILE ) > 0) {
    split(LINE,splitLINE,/\t| /)
    fiveUTR[splitLINE[qName]]=LINE
  }
  FILE=locTMP "three_prime_UTR.psl"
  while((getline LINE < FILE ) > 0) {
    split(LINE,splitLINE,/\t| /)
    threeUTR[splitLINE[qName]]=LINE
  }
  FILE=locTMP "CDS.psl"
  while((getline LINE < FILE ) > 0) {
    split(LINE,splitLINE,/\t| /)
    n=split(splitLINE[10],X,/-P/)
    ID=X[1]
    for(i=2;i<n;i++){
      ID=ID"-P"X[i]
    }
    ID=ID"-R"X[n]
    CDS[ID][splitLINE[11]]=LINE 
  }
}
{
  if($4 in CDS){
    for (POS in CDS[$4]){
      split(CDS[$4][POS],splitCDS,/ |\t/)
      if( splitCDS[tStart]>=$2 && splitCDS[tEnd]<=$3) {
        split(fiveUTR[$4],split5UTR,/ |\t/)
        split(threeUTR[$4],split3UTR,/ |\t/)
        
        if(splitCDS[qStart]==0  && splitCDS[qEnd]==splitCDS[qSize]){
          CDSstat="fullCDS"
        }else{
          if($6=="+"){
            if(splitCDS[qStart]>0  && splitCDS[qEnd]<splitCDS[qSize]){
              CDSstat="truncated-both"
            }else{
              if(splitCDS[qStart]>0  && splitCDS[qEnd]==splitCDS[qSize]){
                CDSstat="truncated-5end"
              }else{
                CDSstat="truncated-3end"
              }
            }
          }else{
            if(splitCDS[qStart]>0  && splitCDS[qEnd]<splitCDS[qSize]){
              CDSstat="truncated-both"
            }else{
              if(splitCDS[qStart]>0  && splitCDS[qEnd]==splitCDS[qSize]){
                CDSstat="truncated-3end"
              }else{
                CDSstat="truncated-5end"
              }
            }
          }
        }
        $7=splitCDS[tStart]
        $8=splitCDS[tEnd]
      }else{
        CDSstat="noCDS"
      }
    }
  }else{
    CDSstat="noCDS"
  }
  if(CDSstat==""){CDSstat="noCDS"}
  
  print  $0,CDSstat
  print $4, CDSstat > locTMP "CDSstat.txt"
  CDSstat=""
}' ${locTMP}mapped.fused.bed >${locTMP}mapped+UTR.bed


#?@ awk -v OFS="\t" -v locTMP=${locTMP} '
#?@ BEGIN{
#?@   FILE=locTMP "mapped_five_prime_UTR.bed"
#?@   while((getline LINE < FILE ) > 0) {
#?@     #split into array by tabs
#?@     split(LINE,splitLINE,/\t| /);
#?@     if(splitLINE[6] == "+") {
#?@       FIVEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[2]]["CHR"]=splitLINE[1]
#?@       FIVEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[2]]["START"]=splitLINE[2]
#?@       FIVEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[2]]["STOP"]=splitLINE[3]
#?@       FIVEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[2]]["STRAND"]=splitLINE[6]
#?@     }else{
#?@       FIVEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[3]]["CHR"]=splitLINE[1]
#?@       FIVEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[3]]["START"]=splitLINE[2]
#?@       FIVEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[3]]["STOP"]=splitLINE[3]
#?@       FIVEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[3]]["STRAND"]=splitLINE[6]

#?@     }
#?@   }
#?@   FILE=locTMP "mapped_three_prime_UTR.bed"
#?@   while((getline LINE < FILE ) > 0) {
#?@     #split into array by tabs
#?@     split(LINE,splitLINE,/\t| /);
#?@     if (splitLINE[6] == "+") {
#?@       THREEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[3]]["CHR"]=splitLINE[1]
#?@       THREEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[3]]["START"]=splitLINE[2]
#?@       THREEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[3]]["STOP"]=splitLINE[3]
#?@       THREEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[3]]["STRAND"]=splitLINE[6]
#?@     }else{
#?@       THREEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[2]]["CHR"]=splitLINE[1]
#?@       THREEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[2]]["START"]=splitLINE[2]
#?@       THREEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[2]]["STOP"]=splitLINE[3]
#?@       THREEprimeUTR[splitLINE[4]"_"splitLINE[1]"_"splitLINE[2]]["STRAND"]=splitLINE[6]

#?@     }
#?@  }
#?@ }
#?@ {
#?@   # print
#?@   # print FIVEprimeUTR[$4"_"$1"_"$2]["START"],FIVEprimeUTR[$4"_"$1"_"$2]["STOP"]
#?@   # print THREEprimeUTR[$4"_"$1"_"$3]["START"],THREEprimeUTR[$4"_"$1"_"$3]["STOP"]
#?@   # print

#?@   if( $6 == "+" ){
#?@     if($1 == FIVEprimeUTR[$4"_"$1"_"$2]["CHR"] && $2 == FIVEprimeUTR[$4"_"$1"_"$2]["START"] && FIVEprimeUTR[$4"_"$1"_"$2]["STOP"]<THREEprimeUTR[$4"_"$1"_"$3]["START"] && FIVEprimeUTR[$4"_"$1"_"$2]["STOP"]<$3){
#?@       $7=FIVEprimeUTR[$4"_"$1"_"$2]["STOP"]
#?@     }
#?@     if($1 == THREEprimeUTR[$4"_"$1"_"$3]["CHR"] && $3 == THREEprimeUTR[$4"_"$1"_"$3]["STOP"] && FIVEprimeUTR[$4"_"$1"_"$2]["STOP"]<THREEprimeUTR[$4"_"$1"_"$3]["START"] && THREEprimeUTR[$4"_"$1"_"$3]["START"]>$2){
#?@       $8=THREEprimeUTR[$4"_"$1"_"$3]["START"]
#?@     }
#?@   }else{
#?@     if($1 == THREEprimeUTR[$4"_"$1"_"$2]["CHR"] && $2 == THREEprimeUTR[$4"_"$1"_"$2]["START"] && THREEprimeUTR[$4"_"$1"_"$2]["STOP"]<FIVEprimeUTR[$4"_"$1"_"$3]["START"] && THREEprimeUTR[$4"_"$1"_"$2]["STOP"]<$3){
#?@       $7=THREEprimeUTR[$4"_"$1"_"$2]["STOP"]
#?@     }
#?@     if($1 == FIVEprimeUTR[$4"_"$1"_"$3]["CHR"] && $3 == FIVEprimeUTR[$4"_"$1"_"$3]["STOP"] && THREEprimeUTR[$4"_"$1"_"$2]["STOP"]<FIVEprimeUTR[$4"_"$1"_"$3]["START"] && FIVEprimeUTR[$4"_"$1"_"$3]["START"]>$2){
#?@       $8=FIVEprimeUTR[$4"_"$1"_"$3]["START"]
#?@     }
#?@   }
#?@   print
#?@ }' ${locTMP}mapped.bed >${locTMP}mapped+UTR.bed

#---------------------------------------------------------------------------------------------------------
#process bed file to final bigbed

#convert to genepred and check genepred to make sure predictions are within chromosom bounds
#?#some warinings will indicate that the CDS boundary is outside of an exon 
#?#this is a problem with only a few nt UTR in the next exon that are not mappped propperly
#?#error should be minimal
bedToGenePred ${locTMP}mapped+UTR.bed ${locTMP}mapped_transcripts.gp
genePredCheck -chromSizes=${HUBdir}/${assemblyNAME}.chrom.sizes ${locTMP}mapped_transcripts.gp

cat ${locTMP}mapped_transcripts.gp > ${locTMP}all.gp

#? #convert back to bed and sort
#? genePredToBed ${locTMP}mapped_transcripts.gp ${locTMP}mapped_transcripts.bed
#? LC_COLLATE=C sort --parallel=$THREADS -k1,1 -k2,2n ${locTMP}mapped_transcripts.bed >${locTMP}mapped_transcripts_sort.bed

#? #generate bigbed file
#? bedToBigBed -as=${UTILITYdir}genes.as -type=bed12+8 -extraIndex=name ${locTMP}mapped_transcripts_sort.bed ${OPENdir}/${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/genes.bb

###################################################################################################
#add tRNA and miRNA annotations

for TYPE in tRNA miRNA ncRNA pseudogene miscRNA; do

  #download sequence file from flybase
  wget -O ${locTMP}${TYPE}.fa.gz  http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.44_FB2022_01/fasta/dmel-all-${TYPE}-r6.44.fasta.gz

  #rename sequences
  gunzip -c ${locTMP}${TYPE}.fa.gz |
    eval $COMMAND |
    awk -v OFS="\t" -v RS=">" -v locTMP=$locTMP -v TYPE=$TYPE '
    {
      sub("\n", "\t")
      gsub("\n", "")
      if(NR>1){
        for(i=1; i<=NF; i++){
          if($i~"name="){
            split($i,splitNAME,/=|;/)
          }
        }

        print $1, splitNAME[2] > locTMP "fbtr_to_name_" TYPE ".txt"
        print ">"splitNAME[2]"\n"$NF
      }
    }' >${locTMP}filtered_${TYPE}.fa

  seqkit split --force --out-dir ${locTMP}${TYPE}_split/ -p ${THREADS} ${locTMP}filtered_${TYPE}.fa

  #map and convert to bed
  #create input IDs for parallel
  IDs=$(seq -w 1 100 | head -n $THREADS)

  #run blat-mapping in parallel
  parallel --linebuffer -j $THREADS "
  blat $assemblyFASTA ${locTMP}${TYPE}_split/filtered_${TYPE}.part_{}.fa -stepSize=5 -fine -noHead -q=dna ${locTMP}${TYPE}_{}.psl" ::: $IDs

  #merge mapped UTRs and convert to bed
  awk -v OFS="\t" '{
    if($1 > ($11 * 0.90)) {
      print
    }
  }' <(cat ${locTMP}${TYPE}_*.psl) >${locTMP}${TYPE}.psl

  pslToBed ${locTMP}${TYPE}.psl ${locTMP}mapped_${TYPE}.bed
 
  #convert to genepred and check genepred to make sure predictions are within chromosom bounds
  bedToGenePred ${locTMP}mapped_${TYPE}.bed ${locTMP}mapped_${TYPE}.gp
  genePredCheck -chromSizes=${HUBdir}/${assemblyNAME}.chrom.sizes ${locTMP}mapped_${TYPE}.gp

  cat ${locTMP}mapped_${TYPE}.gp >> ${locTMP}all.gp

  #? #convert back to bed and sort
  #? genePredToBed ${locTMP}mapped_${TYPE}.gp ${locTMP}mapped_${TYPE}_checked.bed
  #? LC_COLLATE=C sort --parallel=$THREADS -k1,1 -k2,2n ${locTMP}mapped_${TYPE}_checked.bed >${locTMP}mapped_${TYPE}_sorted.bed

  #? #generate bigbed file
  #? bedToBigBed -as=${UTILITYdir}genes.as -type=bed12+8 -extraIndex=name ${locTMP}mapped_${TYPE}_sorted.bed $?{OPENdir}/${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/${TYPE}.bb

done

###################################################################################################
#process all.gp file

LC_COLLAT=C sort -k2,2 -k4,4n ${locTMP}all.gp > ${locTMP}all.sort.gp

grep FlyBase ${TMPdir}FB.gff > ${locTMP}filtered.gff

#prepare final files for UCSC
genePredToBigGenePred ${locTMP}all.sort.gp stdout | LC_COLLATE=C sort -k1,1 -k2,2n | 
awk -v OFS="\t" -v FBgff=${locTMP}filtered.gff -v TMP=${locTMP} -v CDSstatusFile=${locTMP}CDSstat.txt '
  BEGIN{
    while((getline I < FBgff ) > 0) {
      #split into array by tabs
      n=split(I,splitLINE,/\t/);
      #extract FlyBase transcripts and put into array
      if(splitLINE[2]=="FlyBase" ) {
        split(splitLINE[9],splitTAG,/;/)

        if(splitTAG[1]~"ID=FBtr") {
          split(splitTAG[1],splitID,/=/)
          split(splitTAG[2],splitNAME,/=/)
          split(splitTAG[3],splitPARENT,/=/)
          X[splitNAME[2]]["class"]=splitLINE[3]
          X[splitNAME[2]]["name"]=splitID[2]
          X[splitNAME[2]]["parent"]=splitPARENT[2]
          #?@ print splitID[2],splitNAME[2], splitLINE[3], splitTAG[2], splitTAG[3]
        }else{
          if(splitTAG[1]~"ID=FBgn") {
            split(splitTAG[1],splitID,/=/)
            split(splitTAG[2],splitNAME,/=/)
            split(splitTAG[3],splitPARENT,/=/)
            X[splitID[2]]["class"]=splitLINE[3]
            X[splitID[2]]["name"]=splitNAME[2]
            X[splitID[2]]["parent"]=splitPARENT[2]
            #?@ print splitID[2],splitNAME[2], splitLINE[3], splitTAG[2], splitTAG[3]
          }
        }
      }
    }
    while((getline I < CDSstatusFile ) > 0) {
      split(I, splitLINE,/ |\t/)
      CDSstatus[splitLINE[1]]=splitLINE[2]
    }
  }
  {
    NAME=$4

    #remap existing columns
    $16=$15
    $15=$14
    $14=$13

    #FBtr to 2nd ID
    $13=X[NAME]["name"]
    $14="unk"
    $15="unk"

    #add transcript type
    #add transcript type
    if($4~"hpRNA:"){
      $17="hpRNA"
    }else{
      $17=X[NAME]["class"]
    }

    if($4 in CDSstatus){
      CDSstatusTag=CDSstatus[$4]
    }else{
      CDSstatusTag="noCDS"
    }

    if(CDSstatusTag=="fullCDS"){
      $14="cmpl"
      $15="cmpl"
    }
    if(CDSstatusTag=="noCDS"){
      $14="none"
      $15="none"
    }
    if(CDSstatusTag=="truncated-3end"){
      $14="cmpl"
      $15="incmpl"
    }
    if(CDSstatusTag=="truncated-5end"){
      $14="incmpl"
      $15="cmpl"
    }
    if(CDSstatusTag=="truncated-5end"){
      $14="incmpl"
      $15="cmpl"
    }
    if(CDSstatusTag=="truncated-both"){
      $14="incmpl"
      $15="incmpl"
    }

    #add gene ID and name
    if($13~"FBtr"){
      $18=X[X[NAME]["parent"]]["name"]
      $19=X[NAME]["parent"]
      $20=X[X[NAME]["parent"]]["class"]
      $21="FBtranscript_reMapped"
      $22=CDSstatusTag
    }else{
      $18=$4
      $19=$13
      $20=X[NAME]["class"]
      $21="FBtranscript_reMapped"
      $22=CDSstatusTag
    }
    print $0 
  }
'  > ${locTMP}all.named.bigGenePred

#generate bigbed file
#prepare final files for UCSC


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
    string annotationSource;    	"method used for the creation of this annotation"
    string CDSstatus;    	"status of CDS prediction"
    )  
' > ${locTMP}bigGenePred.as

bedToBigBed -tab -as=${locTMP}bigGenePred.as -extraIndex=name -type=bed12+10 ${locTMP}all.named.bigGenePred ${HUBdir}/${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/mapped-genes.bb

###################################################################################################
#? #---------------------------------------------------------------------------------------------------------
#? #create transcript list
#? awk -v OFS="\t" '
#? {
#?   print $4, $1":"$2"-"$3
#? }' ${locTMP}mapped_transcripts.bed >${OPENdir}gene-locations.txt

#---------------------------------------------------------------------------------------------------------
#create trix files

#download reference file
wget -O ${locTMP}fbgn_annotation.tsv.gz http://ftp.flybase.net/releases/FB2022_01/precomputed_files/genes/fbgn_annotation_ID_fb_2022_01.tsv.gz
gunzip -f ${locTMP}fbgn_annotation.tsv.gz

#generate transcriptID to gene info txt
gunzip -c ${locTMP}transcripts.fa.gz |
  awk -v OFS="\t" -v convFILE=${locTMP}fbgn_annotation.tsv '
  BEGIN{
    while((getline CONV < convFILE) > 0) {
        split(CONV,splitCONV,/\t| /)
      {
        CONVarray[splitCONV[3]]=splitCONV[1]
      }
    }
  }
  {
    if($1~">") {
      gsub(";", "")
      n=split($0,splitLINE,/\t| |=/)
      for(i=1; i<=n; i++) {  
        if(splitLINE[i] == "name") {NAMEcol=i+1}
        if(splitLINE[i]=="ID") {IDcol=i+1}
        if(splitLINE[i]=="parent"){PARENTcol=i+1}
      }
      print splitLINE[NAMEcol], splitLINE[IDcol], splitLINE[PARENTcol], CONVarray[splitLINE[PARENTcol]]
    }
  }' | sed 's/;//g;s/ /\t/g' >${locTMP}genes_names.txt

#index file
ixIxx  ${locTMP}genes_names.txt ${HUBdir}/genes.ix ${HUBdir}/genes.ixx

#index for all annotations 
awk '{
  print $4,$13,$18,$19,$17
}' ${locTMP}all.named.bigGenePred >${locTMP}genes_names.all.txt

#index file
ixIxx  ${locTMP}genes_names.all.txt ${HUBdir}/annotations/mapped-genes.ix ${HUBdir}/annotations/mapped-genes.ixx

###################################################################################################
#add to trackDb.txt

#do not continue until file is unblocked by other process
while [[ -f ${TMPdir}wait.txt ]]; do
  sleep 10s
done

#block trackDb from other processes
touch ${TMPdir}wait.txt

#remove old lines if present
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="mapped-transcripts" '
  {
    if( $0 !~ NAME ) print
  }
' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp

mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt

#add repeatMasker relevant tracks
printf "
  track mapped-transcripts
  shortLabel genes-mapped
  longLabel gene annotation created from FlyBase sequences by alignment
  priority 1
  group Debugging
  spectrum on
  vidibility hide
  maxWindowToDraw 10000000
  colorByStrand 50,50,150 150,50,50
  type bigGenePred
  searchIndex name
  searchTrix annotations/mapped-genes.ix
  bigDataUrl annotations/mapped-genes.bb
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
echo "map transcripts using BLAT=" ${PROCESSED_TIME} >>${LOG}time-log.txt

