#! /bin/bash

# $1 : target directory
# $2 : e-value for BLAST
#
# Everything in this script currently must be done from the
# relative path of that directory.

BASEDIR=$( cd ${0%/*} >& /dev/null ; pwd -P )
source "${BASEDIR}/config.sh"

make_fasta="${BASEDIR}/blast/makeFASTA"
miRNA_query="${BASEDIR}/blast/miRNABase_homosapeins.txt"
u12db_dir="${BASEDIR}/blast/u12db"
target=""

if [ $# -ne 2 ]; then
	echo "Usage: run_blast.sh <results_directory> <e-value>"
	exit 1
else
	target=$1
	if [ ! -d $target ]; then
		echo "run_blast.sh: Target directory doesn't exist."
		exit 1
	fi
fi
e_value=$2

cd $target
$make_fasta $target B
#$make_fasta $target C <-- notNovel only
#$make_fasta $target D <-- Novel only

fasta_dir="${target}/$(basename $target)"
results_fasta=$(basename $(ls ${fasta_dir} | grep -E 'FASTA.new$'))
results_fasta="${fasta_dir}/${results_fasta}"

blastn -query $miRNA_query -db ${results_fasta} -out ${target}/miRNA.hits -outfmt 6 -task blastn -reward 1 -dust yes -penalty -1 -gapopen 2 -gapextend 2 -evalue $e_value

for fl in "u12db_3pFlank_u12" "u12db_3pFlank_u2" "u12db_3pFull_u12" "u12db_3pFull_u2"
	do
		blastn -query ${u12db_dir}/${fl} -db ${results_fasta} -out ${target}/${fl}.hits -outfmt 6 -task blastn -reward 1 -dust yes -penalty -1 -gapopen 2 -gapextend 2 -evalue $e_value
done

for fl in "u12db_5pFlank_u12" "u12db_5pFlank_u2" "u12db_5pFull_u12" "u12db_5pFull_u2"
	do
		blastn -query ${u12db_dir}/${fl} -db ${results_fasta} -out ${target}/${fl}.hits -outfmt 6 -task blastn -reward 1 -dust yes -penalty -1 -gapopen 2 -gapextend 2 -evalue $e_value
done

for fl in "u12db_branch_extend_u12" "u12db_branch_u12"
	do
		blastn -query ${u12db_dir}/${fl} -db ${results_fasta} -out ${target}/${fl}.hits -outfmt 6 -task blastn -reward 1 -dust yes -penalty -1 -gapopen 2 -gapextend 2 -evalue $e_value
done

rm -fr $(basename $target)
