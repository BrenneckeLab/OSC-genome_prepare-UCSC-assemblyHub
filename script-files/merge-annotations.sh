#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=20g
#SBATCH --partition=c
#SBATCH -e "%x.e.%A-%a.txt"
#SBATCH -o "%x.o.%A-%a.txt"   
#SBATCH --qos=short
#SBATCH --time=3:00:00


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
locTMP=${TMPdir}merge-annotations/

#create directories
mkdir -p $locTMP

#load tools
source ${SCRIPTdir}tools

THREADS=$(( $SLURM_CPUS_PER_TASK * 2 ))

###################################################################################################

echo a
grep beat-VII ${TMPdir}liftAnnotations/liftoff.named.bigGenePred
echo b
grep beat-VII  ${TMPdir}map_transcripts/all.named.bigGenePred


awk -v OFS="\t" -v LIFTED=${TMPdir}liftAnnotations/liftoff.named.bigGenePred -v MAPPED=${TMPdir}map_transcripts/all.named.bigGenePred '
  BEGIN{
    while((getline LINE < LIFTED ) > 0) {
      #split into array by tabs
      n=split(LINE,splitLINE,/\t/);
      ARRAYlifted[splitLINE[13]]=LINE
      ARRAYall[splitLINE[13]]+=1
    }
    while((getline LINE < MAPPED ) > 0) {
      #split into array by tabs
      n=split(LINE,splitLINE,/\t/);
      ARRAYmapped[splitLINE[13]]=LINE
      ARRAYall[splitLINE[13]]+=1
    }
  }
  {
    x=a
  }
  END{
    for(TR in ARRAYall){
      if( (TR in ARRAYlifted) && (TR in ARRAYmapped) ){
        split(ARRAYlifted[TR],splitLIFTED,/\t| /)
        split(ARRAYmapped[TR],splitMAPPED,/\t| /)

        n=split(splitLIFTED[11],EXONSlifted,/,/)
        for(i=1;i<=n;i++){LENGTHlifted+=EXONSlifted[i]}

        m=split(splitMAPPED[11],EXONSmapped,/,/)
        for(i=1;i<=m;i++){LENGTHmapped+=EXONSmapped[i]}

        if(LENGTHlifted <  LENGTHmapped){
          print ARRAYmapped[TR],"longer_in_mapped"
        }else{
          if(LENGTHlifted >  LENGTHmapped){
            print ARRAYlifted[TR],"longer_in_lifted"
          }else{
            print ARRAYlifted[TR],"same-length_selected_lifted"
          }
        }
        LENGTHlifted=0
        LENGTHmapped=0
      }else{
        if( TR in ARRAYlifted ) {
            print ARRAYlifted[TR],"only_present_lifted"
        }else{
            print ARRAYmapped[TR],"only_present_mapped"
        }
      }
    }
  }' <(echo a) | LC_COLLATE=C sort -k1,1 -k2,2n > ${locTMP}merged.bigGenePred

cut -f 22 ${locTMP}merged.bigGenePred | sort | uniq -c | sort -k1,1rn > ${HUBdir}/annotations/annotationMerge.stats.txt

###################################################################################################
#create track

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
    string mergeReason;    	"reason for chosing this variant during annotation-merging"
   )  
' > ${locTMP}bigGenePred.as

bedToBigBed -extraIndex=name -tab -as=${locTMP}bigGenePred.as -type=bed12+12 ${locTMP}merged.bigGenePred ${HUBdir}${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/HQmergedAnnotations.bb

###################################################################################################
#create trix files

#index for all annotations 
awk '{
  print $4,$13,$18,$19,$17
}' ${locTMP}merged.bigGenePred | tr ' ' '\t'  >${locTMP}genes_names.all.txt

#index file
ixIxx  ${locTMP}genes_names.all.txt ${HUBdir}/annotations/merged-genes.ix ${HUBdir}/annotations/merged-genes.ixx


###################################################################################################
#add to trackDb.txt

#do not continue until file is unblocked by other process
while [[ -f ${TMPdir}wait.txt ]]; do
  sleep 10s
done

#block trackDb from other processes
touch ${TMPdir}wait.txt

#remove old lines if present
awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="mergedAnnotations" '
  {
    if( $0 !~ NAME ) print
  }
' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp

mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt

printf "
track HQmergedAnnotations
shortLabel HQ-Annotations
longLabel high quality transcript annotations merged from different approaches
group Annotations
priority 10
visibility pack
maxWindowToDraw 10000000
mouseOver source \$annotationSource reason \$mergeReason
colorByStrand 50,50,150 150,50,50
type bigGenePred 
searchIndex name
searchTrix annotations/merged-genes.ix
bigDataUrl annotations/HQmergedAnnotations.bb
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

