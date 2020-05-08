#!/bin/bash -l

## VARIABLES
QUEUE="weka"
TOP_DIR=$(pwd)
PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## Threads
threads=8
threadstring="-j \$SLURM_JOB_CPUS_PER_NODE"

#PYTHON_CMD="/gpfs0/apps/x86/anaconda3/bin/python"

# Usage and commands
usageHelp="Usage: ${0##*/} [-d TOP_DIR] [-t THREADS] -h"
dirHelp="* [TOP_DIR] is the top level directory (default \"$TOP_DIR\")\n\
  [TOP_DIR]/fastq must contain the fastq files"
threadsHelp="* [THREADS] is the number of threads for BWA alignment"

while getopts "d:t:h" opt; do
    case $opt in
	h) printHelpAndExit 0;;
        d) TOP_DIR=$OPTARG ;;
	t) threads=$OPTARG 
	   threadstring="-t $threads"
	   ;;
	[?]) printHelpAndExit 1;;
    esac
done

# We assume the files exist in a fastq directory
FASTQ_DIR=${TOP_DIR}"/fastq/*_R*.fastq*"
READ1_STR="_R1"
READ2_STR="_R2"

if [ ! -d "$TOP_DIR/fastq" ]; then
    echo "Directory \"$TOP_DIR/fastq\" does not exist."
    echo "Create \"$TOP_DIR/fastq\" and put fastq files to be aligned there."
    printHelpAndExit 1
else
    if stat -t ${FASTQ_DIR} >/dev/null 2>&1
    then
        echo "(-: Looking for fastq files...fastq files exist"
        testname=$(ls -l ${FASTQ_DIR} | awk 'NR==1{print $9}')
        if [ "${testname: -3}" == ".gz" ]
        then
            read1=${TOP_DIR}"/fastq/*${READ1_STR}*.fastq.gz"
        else
            read1=${TOP_DIR}"/fastq/*${READ1_STR}*.fastq"
        fi
    else
        echo "***! Failed to find any files matching ${FASTQ_DIR}"
	printHelpAndExit 1
    fi
fi

read1files=()
read2files=()
for i in ${read1}
do
    ext=${i#*$READ1_STR}
    name=${i%$READ1_STR*}
    # these names have to be right or it'll break                                                                            
    name1=${name}${READ1_STR}
    name2=${name}${READ2_STR}

    file1=$name1$ext
    file2=$name2$ext
    jid=`sbatch <<- EOF | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -n $threads
	#SBATCH --time=00:30:00
	#SBATCH -p $QUEUE
	#SBATCH -o ${name}_preprocess.out
	#SBATCH -e ${name}_preprocess.err

	cutadapt -a CTGTCTCTTATACACATCT -o ${file1%%.fastq.gz}.trim.fastq.gz $threadstring $file1
	cutadapt -a CTGTCTCTTATACACATCT -o ${file2%%.fastq.gz}.trim.fastq.gz $threadstring $file2
	cutadapt -g file:${PIPELINE_DIR}/primers/artic_primers.fasta -G file:${PIPELINE_DIR}/primers/artic_primers_rc.fasta -o ${file1%%.fastq.gz}.primers.trim.fastq.gz -p ${file2%%.fastq.gz}.primers.trim.fastq.gz $threadstring ${file1%%.fastq.gz}.trim.fastq.gz ${file2%%.fastq.gz}.trim.fastq.gz 
EOF`
    echo "Submitted job $jid to preprocess $file1 and $file2"
done
