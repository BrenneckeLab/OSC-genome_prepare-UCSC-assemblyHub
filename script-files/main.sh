#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH -e "%x.e.%j.txt"
#SBATCH -o "%x.o.%j.txt"
#SBATCH --qos=short
#SBATCH --time=1:00:00

hostname
set -ux

###################################################################################################
#extract variables
VARI=$(echo "$1" | sed 's/,/\t/g;s/"//g')
eval "$VARI"
VARI=$1

#report all variables to log
echo $1 | tr ',' '\n'

#initiate time-variable
TIME=$(date "+%s")

echo "started annotation of assembly" >${LOG}time-log.txt
###################################################################################################
#setup-phase
HUBdir=${OPENdir}${assemblyNAME}/
VARI=${VARI},HUBdir=${HUBdir}
mkdir -p $HUBdir

#load tools
source ${SCRIPTdir}tools

###################################################################################################
#initiate hub  - performed in any case as it is required to pre-process the input assembly for all downstram tools

#create chromosome length file
if [[ ! -s ${HUBdir}/${assemblyNAME}.chrom.sizes || $FORCE == Y ]]; then
  #add genome sequence to the hub and remove ":" as this can cause problems in certain tools
  seqkit replace --line-width 0 -p "\:|\+" -r "_" $assemblyFASTA |
    #also remove name-exstentions from scaffolding or similar as this can cause problems as well
    seqkit seq -i --line-width 0 >${HUBdir}/${assemblyNAME}.fa

  assemblyFASTA=${HUBdir}/${assemblyNAME}.fa
  VARI=${VARI},assemblyFASTA=${assemblyFASTA}

  seqkit fx2tab -n -l ${assemblyFASTA} | sort -k1,1 | tr -s '\t' '\t' >${HUBdir}/${assemblyNAME}.chrom.sizes
fi

#reset assemblyFASTA to converted version 
assemblyFASTA=${HUBdir}/${assemblyNAME}.fa
VARI=${VARI},assemblyFASTA=${assemblyFASTA}

#only if hub does not already exist or should be wiped
if [[ (! -s ${HUBdir}/${assemblyNAME}.2bit || $FORCE == Y) ]]; then
  #clean directory and create path and file structure
  #@ rm -rf ${OPENdir}/
  mkdir -p ${HUBdir}/annotations
  touch ${HUBdir}/trackDb.txt
  
  #convert fasta to twobit format
  faToTwoBit $assemblyFASTA ${HUBdir}/${assemblyNAME}.2bit

  #populate text files for hub
  ${SCRIPTdir}populate_hub.sh $VARI

  #add scaffold-gaps
  twoBitInfo ${HUBdir}/${assemblyNAME}.2bit -nBed stdout | LC_COLLATE=C sort -k1,1 -k2,2n > ${TMPdir}/scaffold_gaps.bed
  
  if [[ -s ${TMPdir}/scaffold_gaps.bed ]]; then 
    bedToBigBed ${TMPdir}/scaffold_gaps.bed ${HUBdir}/${assemblyNAME}.chrom.sizes ${HUBdir}/annotations/scaffold_gaps.bb
    #do not continue until file is unblocked by other process
    while [[ -f ${TMPdir}wait.txt ]]; do
      sleep 10s
    done

    #block trackDb from other processes
    touch ${TMPdir}wait.txt

    #remove old lines if present
    awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="scaffold_gaps" '
      {
        if( $0 ~ NAME || $0 ~ "GAPs" ) {a=b; } else print
      }
    ' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp
    mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt

printf "
track scaffold_gaps
type bigBed 3
shortLabel scaffold_gaps
longLabel scaffold_gaps in the assembly - filled with Ns
maxWindowToDraw 10000000
maxItems 1000000
visibility pack
group Annotations
bigDataUrl annotations/scaffold_gaps.bb
" >>${HUBdir}/trackDb.txt

    rm -rf ${TMPdir}wait.txt
  fi
fi

###################################################################################################
#generate gene annotations
ID=""

# #download transcript-fasta file and pre-process
if [[ ! -s ${TMPdir}FB.gff ]]; then 
  wget --no-verbose -O ${TMPdir}FB.gff.gz http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.44_FB2022_01/gff/dmel-all-r6.44.gff.gz
  gunzip -c ${TMPdir}FB.gff.gz | 
    awk -v OFS="\t" '{if($1~"^##" || $2=="FlyBase") print}' > ${TMPdir}FB.gff 
fi

#---------------------------------------------------------------------------------------------------------
#lift over flybase annotations using liftover
mkdir ${LOG}create-gene-tracks/

if [[ ( -z $STAGE || $STAGE == *gene-liftOVER* ) ]]; then
  #prepare log directory
  rm -rf ${LOG}create-gene-tracks/liftAnnotations*
  cd ${LOG}create-gene-tracks/

  #map transcripts to contigs
  COMMAND="${SCRIPTdir}liftAnnotations.sh"
  cd ${LOG}create-gene-tracks/

  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch --parsable $COMMAND ${VARI})
    ID="${ID}:${newID}"
  else
    $COMMAND ${VARI}
  fi
fi

#---------------------------------------------------------------------------------------------------------
#map flybase transcripts

if [[ ( -z $STAGE || $STAGE == *gene-mapping* ) ]]; then
  #prepare log directory
  rm -rf ${LOG}create-gene-tracks/map_transcripts*
  cd ${LOG}create-gene-tracks/

  #map transcripts to contigs
  COMMAND="${SCRIPTdir}map_transcripts.sh"
  cd ${LOG}create-gene-tracks/

  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch --parsable $COMMAND ${VARI})
    ID="${ID}:${newID}"
  else
    $COMMAND ${VARI}
  fi
fi

#---------------------------------------------------------------------------------------------------------
#de-novo annotate using stringtie
if [[  ( -z $STAGE || $STAGE == *gene-deNovo* ) ]]; then
  mkdir -p ${LOG}create-gene-tracks/de-novo-annotation
  rm -rf ${LOG}create-gene-tracks/de-novo-annotation/*
  cd ${LOG}create-gene-tracks/de-novo-annotation/

  locTMP=${TMPdir}de-novo-annotation/
  mkdir $locTMP

  #map transcripts to contigs
  COMMAND="${SCRIPTdir}de-novo-annotations.sh"
  cd ${LOG}create-gene-tracks/de-novo-annotation

  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch --parsable $COMMAND ${VARI})
    if [[ -z ${ID+x} ]]; then
      ID=$newID
    else
      ID="${ID}:${newID}"
    fi
  else
    $COMMAND ${VARI} 
  fi
fi

#---------------------------------------------------------------------------------------------------------
#unify annotations
#prepare log directory
rm -rf ${LOG}create-gene-tracks/merge-annotations*
cd ${LOG}create-gene-tracks/


if [[ -z $ID ]]; then DEPEND=""; else DEPEND="--dependency=afterok$ID"; fi

if [[ ( -z $STAGE || $STAGE == *gene-merging* ) ]]; then
  #merge annotations
  #? at the moment only lifted and mapped flybase annotations used for merging to final HQ annotations

  #map transcripts to contigs
  COMMAND="${SCRIPTdir}merge-annotations.sh"
  cd ${LOG}

  if [[ $COMPUTING == C ]]; then
    sbatch $DEPEND --parsable $COMMAND ${VARI}
  else
    $COMMAND ${VARI}
  fi
fi

###################################################################################################
#SNP prediction
if [[ ( ! -s ${HUBdir}annotations/SNV/ONT.cuteSV.hom.vcf.gz || ! -s ${HUBdir}annotations/SNV/ILL.deepvariant.hom.vcf.gz || $FORCE == Y ) && ( -z $STAGE || $STAGE == *VARIANTS_* ) ]]; then
  #prepare log directory
  mkdir $LOG/SNV/
  #@ rm -rf $LOG/SNV/*
  cd $LOG/SNV/

  #preset trackDb if required
  TEST=""
  TEST=$(grep -w "track VARIATIONS" ${HUBdir}trackDb.txt)
  if [[ -z $TEST ]]; then

    mkdir   ${HUBdir}annotations/SNV/

    while [[ -f ${TMPdir}wait.txt ]]; do
    sleep 10s
    done

    #block trackDb from other processes
    touch ${TMPdir}wait.txt

    printf "
    track VARIATIONS
    compositeTrack on
    priority 30
    shortLabel Variation tracks
    longLabel Variation tracks - SNV + SV from ONT and Illumina reads
    type vcfTabix
    group Annotations
    visibility squish
    allButtonPair on
    " | sed 's/  //g' >> ${HUBdir}trackDb.txt

    rm -rf ${TMPdir}wait.txt
  fi
  
  #ONT
  if [[ ( -z $STAGE || $STAGE == *VARIANTS_ONT* ) ]]; then
    COMMAND="${SCRIPTdir}predict-SNPs.ONT.sh"

    if [[ $COMPUTING == C ]]; then
      newID=$(sbatch --parsable $COMMAND ${VARI})
      if [[ -z ${ID+x} ]]; then
        ID=$newID
      else
        ID="${ID}:${newID}"
      fi
    else
      $COMMAND ${VARI}
    fi
  fi

  #Illumina
  if [[ ( -z $STAGE || $STAGE == *VARIANTS_ILLUMINA* ) ]]; then
    COMMAND="${SCRIPTdir}predict-SNPs.illumina.sh"

    if [[ $COMPUTING == C ]]; then
      newID=$(sbatch --parsable $COMMAND ${VARI})
      if [[ -z ${ID+x} ]]; then
        ID=$newID
      else
        ID="${ID}:${newID}"
      fi
    else
      $COMMAND ${VARI}
    fi
  fi
fi

###################################################################################################
#uniqueness tracks

if [[ ( $FORCE == Y || ! -s ${HUBdir}/annotations/uniqueness/25mer.3MM.bb) && ( -z $STAGE || $STAGE == *uniqueness* ) ]]; then
  #preset LOG
  mkdir ${LOG}uniqueness
  rm -rf ${LOG}uniqueness/*
  cd ${LOG}uniqueness/

  #test if composit track is already present in trackDb and if not add it
  TEST=""
  TEST=$(grep -w "track uniqueness" ${HUBdir}trackDb.txt)
  if [[ -z $TEST ]]; then

    mkdir   ${HUBdir}annotations/uniqueness

    while [[ -f ${TMPdir}wait.txt ]]; do
    sleep 10s
    done

    #block trackDb from other processes
    touch ${TMPdir}wait.txt

    printf "
    track uniqueness
    compositeTrack on
    shortLabel unique regions
    longLabel unique-regions in the genome predicted with different read-size and MM settings
    type bigBed
    group Annotations
    colorByStrand 110,118,73 132,109,116
    visibility squish
    allButtonPair on
    " | sed 's/  //g' >> ${HUBdir}trackDb.txt

    rm -rf ${TMPdir}wait.txt
  fi


  #run generation of genome uniqueness calculation
  COMMAND="${SCRIPTdir}uniqueness.sh"

  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch --array=25,50,100 --parsable $COMMAND ${VARI})
    if [[ -z ${ID+x} ]]; then
      ID=$newID
    else
      ID="${ID}:${newID}"
    fi
  else
    $COMMAND ${VARI},SLURM_ARRAY_TASK_ID=25
  fi
fi
###################################################################################################
#map raw-nanopore reads

if [[ (( -n $ONT_DNA && ! -s ${HUBdir}annotations/mapped_ONT.prim.bb ) || $FORCE == Y ) && ( -z $STAGE || $STAGE == *map-ONT_DNA* )  ]]; then
  #prepare log directory
  rm -rf ${LOG}map_reads*

  COMMAND="${SCRIPTdir}map_ONT.DNA.sh"
  cd ${LOG}

  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch --parsable $COMMAND ${VARI})
    if [[ -z ${ID+x} ]]; then
      ID=$newID
    else
      ID="${ID}:${newID}"
    fi
  else
    $COMMAND ${VARI}
  fi
fi

###################################################################################################
#map raw-nanopore reads

if [[ (( -n $ONT_RNA ) || $FORCE == Y ) && ( -z $STAGE || $STAGE == *map-ONT_RNA* )  ]]; then

  TEST=""
  TEST=$(grep -w "track mappedONT_RNA" ${HUBdir}trackDb.txt)
  if [[ -z $TEST ]]; then
    #do not continue until file is unblocked by other process
    while [[ -f ${TMPdir}wait.txt ]]; do
      sleep 10s
    done

    #block trackDb from other processes
    touch ${TMPdir}wait.txt

    #remove old lines if present
    awk -v FS="\n" -v RS="\n\n" -v OFS="\t" -v ORS="\n\n" -v NAME="mappedONT_RNA" '
      {
        if( $0 ~ NAME ) {a=b; } else print
      }
    ' ${HUBdir}/trackDb.txt >${TMPdir}trackDB.tmp
    mv ${TMPdir}trackDB.tmp ${HUBdir}/trackDb.txt

    #add tracks to trackDb

    printf "

    track mappedONT_RNA
    superTrack on hide
    shortLabel ONT_RNAseq
    longLabel mapped OSC Nanopore RNAseq reads
    group Long-Reads
    " >> ${HUBdir}/trackDb.txt

    rm -rf ${TMPdir}wait.txt
  fi

  #prepare log directory
  rm -rf ${LOG}map_ONT.RNA*
  mkdir ${HUBdir}annotations/ONT/

  ONT_RNA=$(echo $ONT_RNA | tr '~' '\t')
  nONT_RNA=$(echo $ONT_RNA | wc -w)

  COMMAND="${SCRIPTdir}map_ONT.RNA.sh"
  cd ${LOG}

  if [[ $COMPUTING == C ]]; then
    sbatch --array=1-${nONT_RNA} $COMMAND ${VARI}
  else
    $COMMAND ${VARI}
  fi
fi

###################################################################################################
#run RepeatMasker over the genome

if [[ ( ! -f ${OPENdir}/${assemblyNAME}/annotations/${assemblyNAME}.rmsk.DNA.bb || $FORCE == Y ) && ( -z $STAGE || $STAGE == *RepeatMasker* ) ]]; then
  #prepare log directory
  rm -rf ${LOG}mask_repeats*

  COMMAND="${SCRIPTdir}mask_repeats.sh"
  cd ${LOG}

  echo $COMPUTING $SYSTEM
  if [[ $COMPUTING == C ]]; then
    newID=$(sbatch --parsable $COMMAND ${VARI})
    if [[ -z ${ID+x} ]]; then
      ID=$newID
    else
      ID="${ID}:${newID}"
    fi
  else
    $COMMAND ${VARI}
  fi
fi



###################################################################################################
