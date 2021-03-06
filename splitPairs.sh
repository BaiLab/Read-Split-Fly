#!/bin/bash

BASEDIR=$( cd ${0%/*} >& /dev/null ; pwd -P )
source "${BASEDIR}/config.sh"

OPTSFILE="" #to be filled-in
OUTPUTFILE="" #also to be filled-in

# options file needs:
# RNA-seq file: $2
# max split distance: $6
# rna sample length: $4
# known gene reference: $3
# reference boundries: $3.intronBoundries.exonsgaps
# supporting read tolerance: $7
# output basename: $(basename $2)
function make_options_file() {
    if [ -f "$OPTSFILE" ]; then
        rm "$OPTSFILE"
    fi
    if [ -z $REFDIR ];then
        if [ -n "$BOWTIE_INDEX_ROOT" ]; then
            REFDIR="$BOWTIE_INDEX_ROOT/$1"
        else
            die "Panic! dunno where the ref files are!" 
        fi
    fi
    if [ ! -f "${REFDIR}/${1}.refFlat.txt" ]; then
        die "Error: cannot find refFlat file for $1" 
    fi
    if [ ! -f "${REFDIR}/${1}.refFlat.txt.intronBoundary.exonsgaps" ]; then
        die "Error: cannot find intron/exon boundry file for $1" 
    fi
	touch "$OPTSFILE"
    echo "$2"  > $OPTSFILE  #reads file
    echo "$6" >> $OPTSFILE  #maxSplitDistance
    echo "$3" >> $OPTSFILE  #sampleLength
    echo "${REFDIR}/${1}.refFlat.txt" >> $OPTSFILE  #refFlat
    echo "${REFDIR}/${1}.refFlat.txt.intronBoundary.exonsgaps" >> $OPTSFILE  #intron/exon boundry
    echo "$5" >> $OPTSFILE  #minSplitDistance
    echo "$7" >> $OPTSFILE  #Support tolerance
    if [ -z "$OUTPUTFILE" ]; then
        OUTPUTFILE="${9}/$(echo "$(python $BASENAME_SCRIPT $2)" | cut -d. -f1)"
    fi
    echo "$OUTPUTFILE" >> $OPTSFILE  #results base name
    echo "$8" >> $OPTSFILE   #required supports
}

function dry_run() {
    if [ -f "$OPTSFILE" ]; then
        cat "$OPTSFILE"
    fi
}


function run_rsw() {
    if [ -f "$OPTSFILE" ]; then
        log "OUTPUTFILE = $OUTPUTFILE"
        if [ -z "$LOG_FILE" ]; then
            logfile="/dev/null"
        else
            logfile="${LOG_FILE}"
        fi
        try $RSR_PROGRAM "$OPTSFILE" >> $logfile
        if [ ! -f "${OUTPUTFILE}.results" ]; then
            log "Panic! rsw failed to generate output file. Check stderr." 
            exit 1
        fi
    fi
}

function cleanup() {
#      OUTPUTFILE=$(tail -2 "$1" | head -1 )
 #     mv ${OUTPUTFILE}.* "$2"
      mv "$OPTSFILE" "$1"
}

if (( $# < 9 )); then
    yell "usage: $0 genome readsFile readLength minSplitSize minSplitdistance maxSplitdistance regionBuffer requiredSupports pathToSaveFesults" 
    die "you had $#" 
    exit 1
fi

date=$(date +%s)

genome=$1
log "fyi: \$2 = $2"
OPTSFILE="$RSR_TEMP_DIR/$(python $BASENAME_SCRIPT $2).${date}.options.txt"
OUTPUTFILE="$9/$(echo "$(python $BASENAME_SCRIPT $2)" | cut -d. -f1)"

make_options_file $*

if [ ! -d "${9}" ];then
    log "${9} didn't exist. something wrong?"
    mkdir "${9}"
fi


run_rsw
cleanup "${9}"
#dry_run "$OPTSFILE"
echo "${OUTPUTFILE}.results"
