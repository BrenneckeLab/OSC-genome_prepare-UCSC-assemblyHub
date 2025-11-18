#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --mem=30g
#SBATCH --partition=c
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
locTMP=${TMPdir}liftAnnotations/

#create directories
mkdir $locTMP

#load tools
source ${SCRIPTdir}tools

###################################################################################################
###################################################################################################
THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))


###################################################################################################
#preparations

#download transcriptome sequence file
wget --no-verbose -O ${locTMP}FB.fa.gz http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.44_FB2022_01/fasta/dmel-all-chromosome-r6.44.fasta.gz
gunzip ${locTMP}FB.fa.gz

#curate gff file
gffread ${TMPdir}FB.gff  -O > ${locTMP}FB.mod.gff 

#generate file defining RNA classes for annotation
printf "tRNA
miRNA
snoRNA
snRNA
rRNA
ncRNA
pre_miRNA
pseudogene
snoRNA
snRNA" >  ${locTMP}types.txt

#define chromosome assignment
printf "2L,2L_RagTag
2R,2R_RagTag
3L,3L_RagTag
3R,3R_RagTag
4,4_RagTag
X,X_RagTag" > ${locTMP}chromosome-asignment.txt 

###################################################################################################
#run liftover
liftoff -g ${locTMP}FB.mod.gff -o ${locTMP}liftoff.out.gff -f ${locTMP}types.txt -u ${locTMP}unmapped.txt -dir ${locTMP}intermediate -p $THREADS -chroms ${locTMP}chromosome-asignment.txt -copies  $assemblyFASTA ${locTMP}FB.fa

#filter weird gff annotations
grep -v insertion_site ${locTMP}liftoff.out.gff > ${locTMP}liftoff.filtered.gff

#fix gff features in the lifted gff
gt gff3 -sort -tidy -retainids  ${locTMP}liftoff.filtered.gff > ${locTMP}liftoff.fixed.gff

#artificially convert pseudogene to mRNA to prevent loss during file conversion
sed -i 's/pseudogene/mRNA/' ${locTMP}liftoff.fixed.gff

#convert gff to genepred format
#allowMinimalGenes to keep non-mRNA features
gff3ToGenePred -geneNameAttr=Name -allowMinimalGenes -bad=${locTMP}bad.gp -unprocessedRootsOut=${locTMP}unproccessed.gff ${locTMP}liftoff.fixed.gff ${locTMP}liftoff.fixed.gp

###################################################################################################
#fix naming and format of the gp file

#reduce gff size for faster processing
grep FlyBase ${TMPdir}FB.gff > ${locTMP}filtered.gff

#prepare final files for UCSC
#this script fixes the naming and brings it in sync with the expected UCSC format
genePredToBigGenePred ${locTMP}liftoff.fixed.gp stdout | LC_COLLATE=C sort -k1,1 -k2,2n | 
awk -v OFS="\t" -v FBgff=${locTMP}filtered.gff -v TMP=${locTMP} '
  BEGIN{
    while((getline I < FBgff ) > 0) {
      #split into array by tabs
      n=split(I,splitLINE,/\t/);
      #extract FlyBase transcripts and put into array
      if(splitLINE[2]=="FlyBase" ) {
        split(splitLINE[9],splitTAG,/;/)

        if(splitTAG[1]~"ID=FBtr|ID=FBgn") {
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
  {
    n=split($4,splitFBTR,/_/)
    #FBtr to 2nd ID
    $13=$4

    #add Name of transcript 
    if(n==1){
      $4=X[splitFBTR[1]]["name"]
    }else{
      $4=X[splitFBTR[1]]["name"]"_"splitFBTR[2]
    }
    #add transcript type
    if(X[splitFBTR[1]]["name"]~"hpRNA:"){
      $17="hpRNA"
    }else{
      $17=X[splitFBTR[1]]["class"]
    }

    #add gene ID and name
    if($13~"FBtr"){
      $18=X[X[splitFBTR[1]]["parent"]]["name"]
      $19=X[splitFBTR[1]]["parent"]
      $20=X[X[splitFBTR[1]]["parent"]]["class"]
      $21="FBtranscript_liftOver"
      if($14 == "cmpl" && $15=="cmpl"){
        $22="fullCDS"
      }else{
        $22="noCDS"
      }
    }else{
      $18=$4
      $19=$13
      $20=X[splitFBTR[1]]["class"]
      $21="FBtranscript_liftOver"
      if($14 == "cmpl" && $15=="cmpl"){
        $22="fullCDS"
      }else{
        $22="noCDS"
      }

    }
    print $0 
    #> TMP X[splitFBTR[1]]["class"] ".gp"
  }
'  > ${locTMP}liftoff.named.bigGenePred

#column definition file for bigGenePred file
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

#create final annotation track
bedToBigBed -tab -as=${locTMP}bigGenePred.as -type=bed12+10 ${locTMP}liftoff.named.bigGenePred ${OPENdir}/${assemblyNAME}/${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/liftoff.bb


###################################################################################################
#add to trackDb.txt

#do not continue until file is unblocked by other process
while [[ -f ${TMPdir}wait.txt ]]; do
  sleep 10s
done

#block trackDb from other processes
touch ${TMPdir}wait.txt

#remove old lines if present
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="liftoff" '
  {
    if( $0 !~ NAME ) print
  }
' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp

mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt

printf "
track Liftoff
shortLabel genes-lifted
longLabel gene annotation created from FlyBase sequences by lifting using liftoff
group Debugging
priority 1
visibility hide
maxWindowToDraw 10000000
colorByStrand 50,50,150 150,50,50
type bigGenePred 
bigDataUrl annotations/liftoff.bb
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

