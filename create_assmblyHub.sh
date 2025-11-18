#!/bin/bash

if [[ $LMOD_SYSHOST == CLIP ]]; then
  SYSTEM="CLIP"
  GRIDsystem="SLURM"
elif [[ ${LMOD_SYSHOST} == IMPIMBA-2 ]]; then
  SYSTEM="ii2"
  GRIDsystem="SLURM"
fi

#!#############################################################################################
#!#############################################################################################
#hard-coded things 

STAGEavail="gene-liftOVER gene-mapping gene-deNovo gene-merging RepeatMasker map-ONT_DNA map-ONT_RNA VARIANTS_ONT VARIANTS_ILLUMINA uniqueness"


#FB version in 
  #main.sh
  #liftAnnotations.sh
  #map_transcripts.sh
#chromosome-assignment
  #liftAnnotations.sh
#!#############################################################################################
#!#############################################################################################
###############################################################################################
set -u

# Argument = -i input -c chunksize -D blast-database -v
usage() {
  cat <<EOF
  usage: $0 options
  
  ###############################################################################
  This tool calculates read counts on a template of choice. It reportes the 
  counts for both, the sense and the antisense strand.
  
  usage: [PATH]/selectUTRs [options] -F 
  
  OPTIONS:
      -h  Show this message
      -A  assembly fasta file to be annotated
            results will be depostied into the very same directory of the fasta file
      -N  Name for assembly
      -V  Version of assembly
      -O  Path to directory to place sequence in
      -r  ONT raw DNA read file for aligning to the genome
      -R  ONT raw RNA read file for de-novo gene prediction using Stringtie
            supply path to the directory containing input fastq files
            only use poly-A trimmed reads
            can be comma separated list of input directories containing fastq files
      -o  ONT raw RNA read file for aligning to the genome to get tracks
            supply path to the directory containing input fastq files
            only use poly-A trimmed reads
            can be comma separated list of input directories containing fastq files
      -i  Illumina DNA reads for SNP calling
      -I  Illumina RNA reads for transcriptome assembly using Stringtie
            comma separate list of uncollapsed fastq files
      -Y  set flag if Y-chromosomal transcripts should be included
      -S  run only particular stages
            $STAGEavail
      -C  set flag for local processing (use only if multiple cores available)
      -D  sed debug mode - does not trigger git commit and tmp-files not deleted
      -F  force re-generation of bowtie-indexes and other fixed files
      -W  wipe all data and start fresh
EOF
}

assemblyFASTA=
assemblyNAME=
VERSION=
outPATH=
ONT_DNA=
ONT_deNOVO=
ONT_RNA=
ILLUMINA_DNAseq=
ILLUMINA_RNAseq=
Yinclude=N
STAGE=
COMPUTING=C
DEBUG=N
FORCE=N
WIPE=N

while getopts ÒhA:N:V:O:r:R:o:i:I:YS:CDFW,Ó OPTION; do
  case $OPTION in
  h)
    usage
    exit 1
    ;;
  A)
    assemblyFASTA=$OPTARG
    ;;
  N)
    assemblyNAME=$OPTARG
    ;;
  V)
    VERSION=$OPTARG
    ;;
  O)
    outPATH=$OPTARG
    ;;
  r)
    ONT_DNA=$OPTARG
    ;;
  R)
    ONT_deNOVO=$OPTARG
    ;;
  o)
    ONT_RNA=$OPTARG
    ;;
  i)
    ILLUMINA_DNAseq=$OPTARG
    ;;
  I)
    ILLUMINA_RNAseq=$OPTARG
    ;;
  Y)
    Yinclude=Y
    ;;
  S)
    STAGE=$OPTARG
    ;;
  C)
    COMPUTING=L
    ;;
  D)
    DEBUG=Y
    ;;
  F)
    FORCE=Y
    ;;
  W)
    WIPE=Y
    ;;
  ?)
    usage
    exit
    ;;
  esac
done

###################################################################################################
#fixed variables

#location of reference genome fasta
refFASTA="dmel/dm6/genome_no-Mito_excl-Y.fa"

#basecalling model DNA
MODEL_DNA=r941_prom_sup_g5014

#path to tmp-storage
if [[ $GRIDsystem == SLURM ]]; then
  TMP=
else
  TMP=
fi

#@!@#
###################################################################################################
#variable-setup

#dates
DATE_OF_DAY=$(date +%F)
FULL_DATE=$(date)
#DATE_OF_DAY=2017-02-14

#extract base-path from input fasta
if [[ -z $assemblyNAME ]]; then
  #usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "       Please provide assemblyNAME in option N!\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
  exit
fi

if [[ -z $outPATH ]]; then
  #usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "       Please provide the path to the output directory in option o!\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
  exit
fi

if [[ $assemblyNAME == *"+"* ]]; then
  #usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "       no special characters like + allowed in assembly name!\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
  exit
fi

if [[ ! -z $STAGE ]]; then
    STAGEavail=$(echo $STAGEavail | tr ' ' ',')

  STAGE=$(echo $STAGE | tr ',' '\t')
  for i in $STAGE; do
    if [[ $STAGEavail != *$i* ]];  then
      #usage
      printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
      printf "       Stage $i does not exist!\n"
      printf "$STAGEavail \n"
      printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
      exit
    fi

  done
fi

#extract base-path from input fasta
if [[ -z $assemblyFASTA ]]; then
  #usage
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
  printf "       Please provide assemblyFASTA in option A!\n"
  printf "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n"
  exit
else
  ASSEMBLYdir=$(dirname $assemblyFASTA)
  assemblySTEP=$(echo $ASSEMBLYdir | awk '{n=split($1,X,/\//); print X[n]} ')
  BASEdir=$(dirname $ASSEMBLYdir)
  TMPdir="${TMP}${assemblyNAME}_${VERSION}/"
  OPENdir="${outPATH}/${assemblyNAME}/"

fi


if [[ $WIPE == Y ]]; then
  rm -rf $OPENdir
  rm -rf $TMPdir
fi

mkdir -p $OPENdir
mkdir -p $TMPdir

###################################################################################################
#preset scripts

#determine script-location
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$(cd -P "$(dirname "$SOURCE")" && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done

#generate final SCRIPT_DIR variables
SCRIPTdir_raw="$(cd -P "$(dirname "$SOURCE")" && pwd)"
SCRIPTdir_raw="${SCRIPTdir_raw}/"
SCRIPTdir="${SCRIPTdir_raw}script-files/"
UTILITYdir="${SCRIPTdir_raw}utility-files/"

#move scripts to TMP-directory
cp -r ${SCRIPTdir} ${TMPdir}
SCRIPTdir=${TMPdir}script-files/
mv ${SCRIPTdir}main.sh ${SCRIPTdir}FINann_${assemblyNAME}.sh

###################################################################################################
#push git version and commit

cd ${SCRIPTdir_raw}
echo ${SCRIPTdir_raw}

if [[ $DEBUG != Y ]]; then

  #ask for commit-message
  while true; do
    read -r -p "Plese specify a commit-message: " msg
    case $msg in
    [Nn]) break ;;
    *)
      commitMESSAGE=$msg
      break
      ;;
    esac
  done

  #commit all changes
  git add .
  if [[ -z $commitMESSAGE ]]; then
    git commit -m "automatic commit on submission"
  else
    git commit -m "$commitMESSAGE"
  fi
  #git push --all
fi

commitID=$(git log -1 --pretty=format:"%h")

###################################################################################################
#ml singularity/2.5.2

SINGULARITYdir=${TMP}singularity/

mkdir -p ${SINGULARITYdir}
cd ${SINGULARITYdir}

wget -O repeatmasker.app https://singularity.vbc.ac.at/brenneckelab/repeatmasker/master/repeatmasker.app
wget -O minimap2.app https://singularity.vbc.ac.at/brenneckelab/minimap2/master/minimap2.app
wget -O basicTools.app https://singularity.vbc.ac.at/brenneckelab/basictools/master/basicTools.app
wget -O mosdepth.app https://dominik-handler.github.io/AP_singu/mosdepth/mosdepth.app
wget -O bedops.app https://dominik-handler.github.io/AP_singu/bedops/bedops.app
wget -O busco.app https://singularity.vbc.ac.at/brenneckelab/busco/master/busco.app
wget -O Liftoff.app https://singularity.vbc.ac.at/brenneckelab/liftoff/master/Liftoff.app
wget -O annotate_my_genome.app https://singularity.vbc.ac.at/brenneckelab/annotate_my_genome/master/annotate_my_genome.app
wget -O clair3.app https://singularity.vbc.ac.at/brenneckelab/clair3/master/clair3.app
wget -O rtg-tools.app https://singularity.vbc.ac.at/brenneckelab/rtg-tools/master/rtg-tools.app
wget -O pepper-deepvariant.app https://singularity.vbc.ac.at/brenneckelab/peper_deepvariant/master/pepper-deepvariant.app
wget -O deepvariant.app https://singularity.vbc.ac.at/brenneckelab/deep-variant/master/deepvariant.app
wget -O sniffles.app https://singularity.vbc.ac.at/brenneckelab/sniffles/master/sniffles.app
wget -O cuteSV.app https://singularity.vbc.ac.at/brenneckelab/cutesv/master/cuteSV.app

###################################################################################################
#submit main-run script
LOG=${OPENdir}/LOGs/
mkdir -p ${LOG}
rm -rf ${LOG}FINann_*.txt

#add settings to LOGs
mkdir -p ${LOG}/settings/
cat ${SCRIPTdir_raw}*.sh |
  awk -v RS="#@!@#" '{if (NR==1) print }' >${LOG}/settings/SETTINGS_${commitID}.log

#clear left-over stop-commands
rm -rf ${TMPdir}wait.txt

#cd to log directory for correct log deposition
cd $LOG

#convert variables for submission
ILLUMINA_DNAseq=$(echo $ILLUMINA_DNAseq | tr ',' '~')
ILLUMINA_RNAseq=$(echo $ILLUMINA_RNAseq | tr ',' '~')
ONT_DNA=$(echo $ONT_DNA | tr ',' '~')
ONT_deNOVO=$(echo $ONT_deNOVO | tr ',' '~')
ONT_RNA=$(echo $ONT_RNA | tr ',' '~')
STAGE=$(echo $STAGE | tr '\t' '~')

COMMAND=${SCRIPTdir}FINann_${assemblyNAME}.sh
VARI="OPENdir=${OPENdir},TMPdir=${TMPdir},LOG=${LOG},COMPUTING=${COMPUTING},DEBUG=${DEBUG},SYSTEM=${SYSTEM},FORCE=${FORCE},SCRIPTdir=${SCRIPTdir},UTILITYdir=${UTILITYdir},SINGULARITYdir=${SINGULARITYdir},refFASTA=${refFASTA},assemblyFASTA=${assemblyFASTA},assemblyNAME=${assemblyNAME},ONT_DNA=${ONT_DNA},ONT_deNOVO=${ONT_deNOVO},ONT_RNA=${ONT_RNA},Yinclude=${Yinclude},ILLUMINA_DNAseq=${ILLUMINA_DNAseq},ILLUMINA_RNAseq=${ILLUMINA_RNAseq},MODEL_DNA=${MODEL_DNA},STAGE=${STAGE}"

if [[ $COMPUTING == C ]]; then
  sbatch $COMMAND ${VARI}
else
  if [[ -z ${SLURM_CPUS_PER_TASK+x} ]]; then
    srun --cpus-per-task=10 --mem=20g --time=1:00:00 --qos=short $COMMAND ${VARI}
  else
    $COMMAND ${VARI}
  fi
fi

exit
