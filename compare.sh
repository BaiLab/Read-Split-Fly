#!/bin/bash

BASEDIR=$( cd ${0%/*}  >& /dev/null ; pwd -P )
source "${BASEDIR}/config.sh"

SUFFIX=".results"
OUTNAME="comparison.txt"

function do_compare() {
    tolerance=$1
    shift
    whereto=$(python $DIRNAME_SCRIPT $2)
    $COMPARE_PROG "$LINES_TO_SKIP" $tolerance $OUTNAME $@ >& /dev/null
    mv "${OUTNAME}.comparisonSummary.txt" "$whereto"
    mv "${OUTNAME}.comparedResults.txt" "$whereto"

}

#---------Main-------

if (( $# < 3 )); then
    yell "usage -- $(python $BASENAME_SCRIPT $0)  <support tolerance> <file1> <file2>"
    die  "you had $# : $@"
fi


do_compare $@
