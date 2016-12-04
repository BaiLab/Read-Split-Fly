#!/bin/bash

#------Only Configure the thing ONCE, please---
if [ -z "$CONFIGURED" ]; then
CONFIGURED=true
DEBUG=false
BASEDIR=$( cd ${0%/*} >& /dev/null ; pwd -P )     #moved to this method to avoid platform-specific comamnds
BASENAME_SCRIPT="${BASEDIR}/src/utils/basedir.py" # Python clone of GNU basename program
DIRNAME_SCRIPT="${BASEDIR}/src/utils/dirname.py"  # Python clone of GNU dirname program
yell() { echo "$(python $BASENAME_SCRIPT $0): $*" >&2; }
die() { yell "$*"; log "$*"; exit 111; kill $$; }
try() { "$@" || die "cannot $*"; }
dprint() { if $DEBUG; then yell "$*"; fi }
log() { echo "$(python $BASENAME_SCRIPT $0): $*" >> "${LOG_FILE}"; dprint "$@"; }
timing_start() { if [ -z "$START_TIME" ]; then START_TIME=$(date +%s); fi }
timing_end() { if [ ! -z "$START_TIME" ]; then log "Duration: $(( $(date +%s) - $START_TIME )) seconds"; fi }


#-------USER CONFIGURATION------------
#-------Misc. Settings----------------
RM_TEMP_FILES=1                             # set =1 to delete all intermediate files, =0 to keep them
#-------Directories-------------------
BOWTIE_INDEXES=${BASEDIR}/bt/indexes        # Location where you store your bowtie indexes.  Comment out to use
                                            #   BOWTIE_INDEXES environment variable instead
BASE_TEMP_DIR="${BASEDIR}/tmp"
BOWTIE_TEMP_DIR="${BASE_TEMP_DIR}/bowtie"
SPLIT_TEMP_DIR="${BASE_TEMP_DIR}/split"
RSR_TEMP_DIR="${BASE_TEMP_DIR}/splitpairs"
LOG_DIR="${BASEDIR}/logs"
BOWTIE_INDEX_ROOT=""
REFDIR="${BOWTIE_INDEXES}"                  #where the refFlats are
#------Programs----------------------
BOWTIE_PROGRAM="${BASEDIR}/bt/bowtie"                 # by default, use the bowtie in bt/. If you want to use your own version of bowtie,


#------NOT-USER CONFIGURATION--------
#CHANGE THESE AT YOUR OWN RISK

#------Programs----------------------
SPLIT_PROGRAM="${BASEDIR}/srr"                        # Program for splitting reads, compiled from split_read_rsr.c
FORMAT_PROGRAM="${BASEDIR}/sfc"                       # Program for formatting reads, compiled from split_first_column.c
RSR_PROGRAM="${BASEDIR}/sp4"                          # RSR Program ("split pairs"), compiled from splitPairs.cpp
COMPARE_PROGRAM="${BASEDIR}/compare"                  # Program for comparing RSR outputs, compiled from 
BOWTIE_INSPECT_RSR="${BASEDIR}/bt/bowtie-inspect-RSR" # Program adds spliced sequences to results file

#Meta defines
if [ -z "$BOWTIE_PROGRAM" ]; then
export -p BOWTIE_PROGRAM="${BASEDIR}/bt/bowtie"
fi

export -p TERM=vt100                                        # this seems to be necessary for bowtie...
if [ -z "$LOG_FILE" ] || [ -z "$RUN_ID" ]; then
export -p RUN_ID="$(date +%F.%s)"
export -p LOG_FILE="${LOG_DIR}/RSF_${RUN_ID}.log"
fi

BASES_TO_TRIM=0
#this sets the location to find the sub-parts of the pipeline

#rsr_batch_job.sh:  The big guy!
PIPELINE="${BASEDIR}/pipeline.sh"
COMPARE_SCRIPT="${BASEDIR}/compare.sh"

#pipeline.sh: Pipeline Constants
#Define program/script files
ALIGN_SCRIPT="${BASEDIR}/bowtie.sh"
MEASURE_SCRIPT="${BASEDIR}/readlength.sh"
SPLIT_SCRIPT="${BASEDIR}/split.sh"
RSR_SCRIPT="${BASEDIR}/splitPairs.sh"   #_SCRIPT'd for consistency

# BLAST Query Stuff
BLAST_DIR="${BASE_TEMP_DIR}/blast"
BLAST_SCRIPT="${BASEDIR}/blast.sh"

#readlength.sh:  measuring tool
OFFSET=$(( $BASES_TO_TRIM + 1))      #1 + the number of bases to trim
#the 1 is to account for the newline when reading from the file.

#bowtie.sh:  Alignment constants
ENCODING_GUESSER="${BASEDIR}/guess-encoding.py"
#BOWTIE_INDEXES will be dynamically set in bowtie.sh
QUALITY_TESTS=1000

#split.sh:   reads splitting
#everything in this section has since been defined elsewhere. -AC

#sp4:   split pair finder
RSR="${BASEDIR}/sp4"
OUTPUTFILE=""    #dummy value; filled in within the scripts

#rsr_compare.sh :  compare the output of two jobs
#TODO: look to see if this is strictly a duplicate of COMPARE_PROGRAM;
# they are both used in different places
COMPARE_PROG="${BASEDIR}/compare"
LINES_TO_SKIP=22

#cleanup.sh:  gets rid of old files
#none. everything is defined above

#Manage directories
if [ ! -d "$BASE_TEMP_DIR" ]; then
    mkdir "$BASE_TEMP_DIR" || die "Could not make $BASE_TEMP_DIR aborting."
fi
if [ ! -d "$BOWTIE_TEMP_DIR" ]; then
    mkdir "$BOWTIE_TEMP_DIR" || die "Could not make $BOWTIE_TEMP_DIR aborting."
fi
if [ ! -d "$SPLIT_TEMP_DIR" ]; then
    mkdir "$SPLIT_TEMP_DIR" || die "Could not make $SPLIT_TEMP_DIR aborting."
fi
if [ ! -d "$RSR_TEMP_DIR" ]; then
    mkdir "$RSR_TEMP_DIR" || die "Could not make $RSR_TEMP_DIR aborting."
fi
if [ ! -d "$LOG_DIR" ]; then
    mkdir "$LOG_DIR" || die "Could not make $LOG_DIR aborting."
fi

#SANITY CHECKs

#first, BOWTIE_INDEXES for no particular reason
if [ -z "$BOWTIE_INDEXES" ]; then
    die "No BOWTIE_INDEXES. cannot continue."
fi

function sanity_check() {
    if [[ ${BASH_SOURCE[0]} != $0 ]]; then

        #check that the input rna-seq files exist
        deficient=()
        if [[ $3 =~ .*\|.* ]]; then # Paired mode
            OIFS=$IFS
            IFS='|' read pair1 pair2 <<< "$3"
            IFS=$OIFS
            if [[ $pair1 =~ .*,.* ]]; then
                OIFS=$IFS
                IFS=','; for file in $pair1; do if [ ! -f "$file" ]; then deficient+="$file "; fi; done
                IFS=$OIFS
            else
                if [ ! -f "$pair1" ]; then deficient+="$pair1 "; fi  #bug fix reported by Feng in 1/18/16 email
            fi
        else
            if [[ $3 =~ .*,.* ]]; then
                OIFS=$IFS
                IFS=','; for file in $3; do if [ ! -f "$file" ]; then deficient+="$file "; fi; done
                IFS=$OIFS
            else
                if [ ! -f "$3" ]; then deficient+="$3 "; fi
            fi
        fi

    fi

    if (( ${#deficient[@]} > 0 )); then die "Could not find the following input files(${#deficient[@]}): ${deficient[@]}"; fi

    #check for the pipeline components
    deficient=()
    yell "sanity_check::\$BOWTIE_PROGRAM=$BOWTIE_PROGRAM"
    if [ ! -f $BOWTIE_PROGRAM ]; then
        bowtie --version >& /dev/null
        if (( $? > 0 )); then 
            deficient+="$BOWTIE_PROGRAM "
        else
            export -p BOWTIE_PROGRAM=$(which bowtie)
        fi
    fi
    yell "sanity_check::\$BOWTIE_PROGRAM=$BOWTIE_PROGRAM"
    if [ ! -f $SPLIT_PROGRAM ]; then deficient+="$SPLIT_PROGRAM "; fi
    if [ ! -f $FORMAT_PROGRAM ]; then deficient+="$FORMAT_PROGRAM " ; fi
    if [ ! -f $RSR_PROGRAM ]; then deficient+="$RSR_PROGRAM " ; fi
    if [ ! -f $COMPARE_PROGRAM ]; then deficient+="$COMPARE_PROGRAM " ; fi
    if [ ! -f $MEASURE_SCRIPT ]; then deficient+="$MEASURE_SCRIPT " ; fi
    if [ ! -f $ALIGN_SCRIPT ]; then deficient+="$ALIGN_SCRIPT " ; fi
    if [ ! -f $SPLIT_SCRIPT ]; then deficient+="$SPLIT_SCRIPT " ; fi
    if [ ! -f $RSR_SCRIPT ]; then deficient+="$RSR_SCRIPT " ; fi
    if [ ! -f $ENCODING_GUESSER ]; then deficient+="$ENCODING_GUESSER "; fi

    if (( ${#deficient[@]} > 0 )); then die "The following components of the pipeline are missing: ${deficient[@]}"; fi
}

export -f sanity_check
#end if..configured
fi
