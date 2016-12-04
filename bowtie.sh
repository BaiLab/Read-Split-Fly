#!/bin/bash

BASEDIR=$( cd ${0%/*} >& /dev/null ; pwd -P )  #moved to this method to avoid platform-specific comamnds
source "${BASEDIR}/config.sh"

function set_indexes() {
    if [ -n "$BOWTIE_INDEX_ROOT" ]; then
        BOWTIE_INDEXES="${BOWTIE_INDEX_ROOT}/$1"
    fi
    export -p BOWTIE_INDEXES
}

function run_bowtie() {
    phase=$1
    bowtie_index_file=$2
    params="$3"
    inputs="$4"
    outfiles=()
    log "Using BOWTIE_INDEXES in $(printenv BOWTIE_INDEXES)";
    if [ -z "$(printenv BOWTIE_INDEXES)" ]; then
        die "No BOWTIE_INDEXES"
    fi
    outname="$BOWTIE_TEMP_DIR/$5"
    if [ "$phase" == "phase1" ]; then
        log "bowtie $bowtie_index_file $QUALS  $params -q $inputs --un ${outname}.unmapped.txt ${outname}.bowtie.txt"
        try $BOWTIE_PROGRAM $bowtie_index_file $QUALS $params -q $inputs --un "$outname".unmapped.txt "$outname".bowtie.txt >> "$LOG_FILE" 2>&1
        #try $BOWTIE_PROGRAM $bowtie_index_file $QUALS $params -q $inputs --un "$outname".unmapped.txt -S "$outname".bowtie.SAM.txt >> "$LOG_FILE" 2>&1
    elif [ "$phase" == "phase2" ]; then
        log "bowtie $bowtie_index_file $QUALS  $params -q $inputs $outname.bowtie.txt"
        try $BOWTIE_PROGRAM $bowtie_index_file $QUALS $params -q $inputs "$outname".bowtie.txt >> "$LOG_FILE" 2>&1
    else
        die "Don't know what to do on phase $phase"
    fi

    #yes, it literally passes $4 back out, with some adornment
    echo $outname
}

#------MAIN-------

if (( $# < 4 )); then
yell "usage -- $0 genome phase fileSet maxGoodAlignments"
yell "      fileSet should be a comma-separated list of fastQ files"
yell "              Optionally, a second set of comma-separated files"
yell "              May be added, separated by a pipe | for paired-end"
exit 1
fi

try set_indexes $1
genome=$1
if (( $? )) && [ -z "$genome" ]; then
    die "Cannot find genome for $1"
fi

#parameter discovery
if [[ $3 =~ .*\|.* ]]; then # Paired mode
    if $DEBUG; then echo "Paired mode" 1>&2; fi
    OIFS=$IFS
    IFS='|' read pair1 pair2 <<< "$3"
    IFS=$OIFS
    input_params="-1 $pair1 -2 $pair2"
    outname=$(python $BASENAME_SCRIPT "$(echo $pair1 | cut -d, -f1)")
    QUALS=$(awk 'NR % 4 == 0' $(echo "$pair1" | cut -d, -f1) | head -$(( $QUALITY_TESTS * 4)) | python $ENCODING_GUESSER -b -n $QUALITY_TESTS)

else
    input_params=$3
    outname=$(python $BASENAME_SCRIPT "$(echo $3 | cut -d, -f1)")
    QUALS=$(awk 'NR % 4 == 0' $(echo "$3" | cut -d, -f1) | head -$(( $QUALITY_TESTS * 4)) | python $ENCODING_GUESSER -b -n $QUALITY_TESTS)
fi

#number of threads is now 1/4 of total cores on the system--this doesn't seem
# to be slower than using 3/4 of total cores
#bowtie_params="-t --chunkmbs 2048 -p $(( $(grep -c ^processor /proc/cpuinfo) * 1 / 4)) "
bowtie_params="-t --chunkmbs 2048 -p $NUM_THREADS "
if [ "$2" == "phase1" ]; then
	bowtie_params+="-n 3 -e 112"
elif [ "$2" == "phase2" ]; then
	bowtie_params+="--best -k $4 -m $4 -v 0"
fi
log "BOWTIE_INDEX_ROOT: $BOWTIE_INDEX_DIR, BOWTIE_INDEXES: $BOWTIE_INDEXES"
log "bowtying $2 $genome $bowtie_params $input_params"

echo $(run_bowtie $2 "$genome" "$bowtie_params" "$input_params" "$outname")

