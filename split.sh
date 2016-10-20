#!/bin/bash

BASEDIR=$( cd ${0%/*} >& /dev/null ; pwd -P )
#bring in the configuration, which should be in the script's directory
source "${BASEDIR}/config.sh"


#inputs: reads_file left_cut reads_length
function n_split() {
    if [ -z "$LOG_FILE" ]; then
        output="/dev/null"
    else
        output=$LOG_FILE
    fi
    $SPLIT_PROGRAM "$1" "$2" "$SPLIT_TEMP_DIR" >> $output
    log "n_split:: ${SPLIT_TEMP_DIR}/$(basename $1).split"
    echo "${SPLIT_TEMP_DIR}/$(basename $1).split"
}

#input: resultsFileBase minimum-cut-size readsLength 
function split_pairs() {
    base=$1
    start=$2
    #no more len, delete in next version
    #if (( $start >= $len / 2 )); then
    #    start=$(( $len - $start ))
    #fi
    log "split_pairs::file=$1,start=$2,TERM=$TERM"
    if [ -f "${base}.unmapped.txt" ]; then
        result=$(try n_split "${base}.unmapped.txt" $start)
    elif [ -f "${base}.unmapped_1.txt" ] && [ -f "${base}.unmapped_2.txt" ]; then
        #tmp=$(try n_split "${base}.unmapped_1.txt" $start $(( $len - $start ))  )  # why was it len-start, does that matter on replicates?
        tmp=$(try n_split "${base}.unmapped_1.txt" $start  )
        mv "$tmp" "${tmp}-1.txt"
        result="${tmp}-1.txt"
        #tmp=$(try n_split "${base}.unmapped_2.txt" $start $(( $len - $start ))  )  # why was it len-start, does that matter on replicates?  why not just len?
        tmp=$(try n_split "${base}.unmapped_2.txt" $start  )
        mv "$tmp" "${tmp}-2.txt"
        result+="|${tmp}-2.txt"
    else
        die "split_pairs::Don't know what do do with $base"
    fi
    log "split_pairs:: $result"
    echo "$result"
}


#-----Main---

if (( $# < 2 )); then
    yell "Usage: $(basename $0) <reads filename base> <min split>"
    yell "The reads filename should be the part before the \"unmapped.txt\" or \"unmapped_1.txt\""
    yell "The script will detect whether there is a pair and split both."
    yell "examples: "
    yell "if your file to split is \"Dtt_KO.txt.unmapped.txt\" and you wanted a minimum split of 11bp, then..."
    yell "         $(basename $0) Dtt_KO.txt 11"
else
    echo $(split_pairs $@)
fi
