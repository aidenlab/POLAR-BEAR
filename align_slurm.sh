#!/bin/bash
##### 
### Polar pipeline script for viral diagnostic
### Given paired end sequencing of putative viral data
###   -> Aligns to a selection of potential viruses, sorts and merges
###   -> In parallel, assembles the data into contigs
###   -> After alignment, runs samtools depth to create stats.csv
###   -> After contigging, creates pairwise alignment
###   -> Final report created as PDF with dotplot from paf and stats
#####

### VARIABLES: SET THESE FOR YOUR SYSTEM

## Software commands; add load modules as necessary
# BWA for alignment to different genomes
LOAD_BWA="spack load bwa@0.7.17 arch=\`spack arch\`"
BWA_CMD="bwa"
# Samtools for file manipulation
LOAD_SAMTOOLS=""
SAMTOOLS_CMD="samtools"
# Megahit for contig creation
LOAD_MEGAHIT=""
MEGAHIT_CMD="/gpfs0/work/brian/scripts/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit"
# Minimap2 for creating dotplot
LOAD_MINIMAP2=""
MINIMAP2_CMD="minimap2"
# Python for creating dotplot
LOAD_PYTHON=""
PYTHON_CMD="/gpfs0/apps/x86/anaconda3/bin/python"

## Queues
QUEUE="weka"
QUEUE_X86="weka"

## Threads
threads=8
threadstring="-t \$SLURM_JOB_CPUS_PER_NODE"

### VARIABLES: Automatically set
TOP_DIR=$(pwd)
PIPELINE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BETACORONA_REF_DIR="${PIPELINE_DIR}/betacoronaviruses/*/*/*.fasta \
                    ${PIPELINE_DIR}/control/*.fasta"
BETACORONA_SMALL="${PIPELINE_DIR}/betacoronaviruses/close/*/*.fasta \
		  ${PIPELINE_DIR}/betacoronaviruses/match/*/*.fasta \
                  ${PIPELINE_DIR}/control/*.fasta" 
# Viral match - assumes only one directory under "match" 
MATCH_REF="${PIPELINE_DIR}/betacoronaviruses/match/*/*.fasta"
MATCH_NAME=$(echo $MATCH_REF | sed 's:.*/::' | rev | cut -c7- | rev )

# Usage and commands
usageHelp="Usage: ${0##*/} [-d TOP_DIR] [-t THREADS] -jkrh"
dirHelp="* [TOP_DIR] is the top level directory (default \"$TOP_DIR\")\n\
  [TOP_DIR]/fastq must contain the fastq files"
threadsHelp="* [THREADS] is the number of threads for BWA alignment"
indexHelp="* -j produce index file for aligned files"
reducedHelp="* -r reduced set for alignment"
stageHelp="* -k start pipeline after alignment"
helpHelp="* -h: print this help and exit"

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$dirHelp"
    echo -e "$threadsHelp"
    echo -e "$indexHelp"
    echo -e "$reducedHelp"
    echo -e "$stageHelp"
    echo "$helpHelp"
    exit "$1"
}

while getopts "d:t:hrjk" opt; do
    case $opt in
	h) printHelpAndExit 0;;
        d) TOP_DIR=$OPTARG ;;
	j) produceIndex=1 ;;
	r) reducedSet=1 ;;
	t) threads=$OPTARG 
	   threadstring="-t $threads"
	   ;;
	k) afteralignment=1 ;;
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
    read1files+=($name1$ext)
    read2files+=($name2$ext)
done

# replace spaces with commas for megahit
read1filescomma=$(echo "${read1files[*]}" | sed 's/ /,/g;s/,$//')
read2filescomma=$(echo "${read2files[*]}" | sed 's/ /,/g;s/,$//')

if [[ "$reducedSet" -eq 1 ]]
then
    REFERENCES=$BETACORONA_SMALL
else
    REFERENCES=$BETACORONA_REF_DIR
fi

WORK_DIR=${TOP_DIR}/work
LOG_DIR=${TOP_DIR}/log
FINAL_DIR=${TOP_DIR}/final


if [[ "$afteralignment" -ne 1 ]]
then
    if ! mkdir "${WORK_DIR}"; then echo "***! Unable to create ${WORK_DIR}! Exiting"; exit 1; fi
    if ! mkdir "${LOG_DIR}"; then echo "***! Unable to create ${LOG_DIR}! Exiting"; exit 1; fi
    if ! mkdir "${FINAL_DIR}"; then echo "***! Unable to create ${FINAL_DIR}! Exiting"; exit 1; fi
    
    for REFERENCE in $REFERENCES
    do
        ######################################################################
        ########## Align 
        ######################################################################
	REFERENCE_NAME=$(echo $REFERENCE | sed 's:.*/::' | rev | cut -c7- | rev )

	echo -e "(-: Aligning files matching $FASTQ_DIR\n to genome $REFERENCE_NAME"

	if ! mkdir "${WORK_DIR}/${REFERENCE_NAME}"; then echo "***! Unable to create ${WORK_DIR}/${REFERENCE_NAME}! Exiting"; exit 1; fi
	if ! mkdir "${WORK_DIR}/${REFERENCE_NAME}/aligned"; then echo "***! Unable to create ${WORK_DIR}/${REFERENCE_NAME}/aligned! Exiting"; exit 1; fi
	if ! mkdir "${WORK_DIR}/${REFERENCE_NAME}/debug"; then echo "***! Unable to create ${WORK_DIR}/${REFERENCE_NAME}/debug! Exiting"; exit 1; fi

	dependsort="afterok"

	for ((i = 0; i < ${#read1files[@]}; ++i)); do
            file1=${read1files[$i]}
            file2=${read2files[$i]}

	    FILE=$(basename ${file1%$READ1_STR*})
	    ALIGNED_FILE=${WORK_DIR}/${REFERENCE_NAME}/aligned/${FILE}"_mapped"
	    
            # Align reads
	    jid=`sbatch <<- ALGNR | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $QUEUE
		#SBATCH -o ${WORK_DIR}/${REFERENCE_NAME}/debug/align-%j.out
		#SBATCH -e ${WORK_DIR}/${REFERENCE_NAME}/debug/align-%j.err
		#SBATCH -t 2880
		#SBATCH -n 1
		#SBATCH -c $threads
		#SBATCH --mem-per-cpu=4G
		#SBATCH -J "align_${FILE}"
		#SBATCH --threads-per-core=1

		$LOAD_BWA
		$BWA_CMD 2>&1 | awk '\\\$1=="Version:"{printf(" BWA %s; ", \\\$2)}'
		echo "Running command $BWA_CMD mem $threadstring $REFERENCE $file1 $file2 > ${ALIGNED_FILE}.sam"
		srun --ntasks=1 $BWA_CMD mem $threadstring $REFERENCE $file1 $file2 > ${ALIGNED_FILE}.sam
		if [ \$? -ne 0 ]                      
		then  
			exit 1                                         
		else  
			echo "(-: Mem align of ${ALIGNED_FILE}.sam done successfully"             
		fi
ALGNR`
	    dependalign="afterok:$jid"

            ######################################################################
            ##########SAM: fixmate, sort
            ######################################################################
	    jid=`sbatch <<- SAMTOBAM | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $QUEUE
		#SBATCH -o ${WORK_DIR}/${REFERENCE_NAME}/debug/sam2bam-%j.out
		#SBATCH -e ${WORK_DIR}/${REFERENCE_NAME}/debug/sam2bam-%j.err
		#SBATCH -t 2880 
		#SBATCH -n 1
		#SBATCH -c $threads
		#SBATCH --mem-per-cpu=4G
		#SBATCH --threads-per-core=1
		#SBATCH -d $dependalign

		$LOAD_SAMTOOLS
		$SAMTOOLS_CMD fixmate -m $ALIGNED_FILE".sam" $ALIGNED_FILE".bam"
		$SAMTOOLS_CMD sort -@ $threads -o $ALIGNED_FILE"_matefixd_sorted.bam" $ALIGNED_FILE".bam"

SAMTOBAM`
	    dependsort="${dependsort}:$jid"
	done
    
        ######################################################################
        ##########SAM: merge, dedup, depth, stats
        ######################################################################
	jid=`sbatch <<- MRKDUPS | egrep -o -e "\b[0-9]+$"
		#!/bin/bash -l
		#SBATCH -p $QUEUE
		#SBATCH -o ${WORK_DIR}/${REFERENCE_NAME}/debug/mrkdups-%j.out
		#SBATCH -e ${WORK_DIR}/${REFERENCE_NAME}/debug/mrkdups-%j.err
		#SBATCH -t 2880 
		#SBATCH -n 1
		#SBATCH -c $threads
		#SBATCH --mem-per-cpu=4G
		#SBATCH --threads-per-core=1
		#SBATCH -d $dependsort

		$LOAD_SAMTOOLS 
		if $SAMTOOLS_CMD merge ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam ${WORK_DIR}/${REFERENCE_NAME}/aligned/*_matefixd_sorted.bam
		then
			rm ${WORK_DIR}/${REFERENCE_NAME}/aligned/*_matefixd_sorted.bam ${WORK_DIR}/${REFERENCE_NAME}/aligned/*_mapped*
		fi

		if $SAMTOOLS_CMD markdup ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged_dups_marked.bam
		then
			rm ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam
		fi
		
		$SAMTOOLS_CMD depth -a ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged_dups_marked.bam > ${WORK_DIR}/${REFERENCE_NAME}/aligned/depth_per_base.txt
		$SAMTOOLS_CMD stats ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged_dups_marked.bam | grep  ^SN | cut -f 2- > ${WORK_DIR}/${REFERENCE_NAME}/aligned/stats.txt
MRKDUPS`
	dependmerge="afterok:$jid"

	if [[ "$REFERENCE" == *match* ]]
	then
	    MATCH_REF=${REFERENCE}
	    MATCH_NAME=${REFERENCE_NAME}
	    dependmatchdone="afterok:$jid"
	fi

	echo $produceIndex
	if [ -n "$produceIndex" ]
	then
	    jid=`sbatch <<- INDEXSAM | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $QUEUE
	#SBATCH -o ${WORK_DIR}/${REFERENCE_NAME}/debug/indexsam-%j.out
	#SBATCH -e ${WORK_DIR}/${REFERENCE_NAME}/debug/indexsam-%j.err
	#SBATCH -t 2880 
	#SBATCH -n 1 
	#SBATCH -c 1
	#SBATCH --mem-per-cpu=10G
	#SBATCH --threads-per-core=1 
	#SBATCH -d $dependmerge

	$LOAD_SAMTOOLS
	$SAMTOOLS_CMD index ${WORK_DIR}/${REFERENCE_NAME}/aligned/sorted_merged.bam

INDEXSAM`
	fi
    done
    slurm_depend_merge="#SBATCH -d $dependmerge"
else
    dependmatchdone="afterok"
    slurm_depend_merge=""
fi

echo "#!/bin/bash -l" > "$WORK_DIR"/collect_stats.sh
echo "#SBATCH -p $QUEUE" >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH -o ${LOG_DIR}/collectstats-%j.out"  >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH -e ${LOG_DIR}/collectstats-%j.err" >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH -t 30" >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH -n 1 " >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH -c 1" >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH --mem=200" >> "$WORK_DIR"/collect_stats.sh
echo "#SBATCH --threads-per-core=1 " >> "$WORK_DIR"/collect_stats.sh
echo "$slurm_depend_merge"  >> "$WORK_DIR"/collect_stats.sh 
echo "echo \"label,percentage\" > $WORK_DIR/stats.csv " >> "$WORK_DIR"/collect_stats.sh
echo "for f in $WORK_DIR/*/aligned/depth_per_base.txt; do"  >> "$WORK_DIR"/collect_stats.sh
echo "if [[ ! \"\$f\" == *\"control\"* ]]; then" >> "$WORK_DIR"/collect_stats.sh 
echo  "awk -v fname=\$(basename \${f%%/aligned*}) 'BEGIN{count=0; onisland=0}\$3>5{if (!onisland){onisland=1; island_start=\$2}}\$3<=5{if (onisland){island_end=\$2; if (island_end-island_start>=50){count=count+island_end-island_start}} onisland=0}END{if (onisland){island_end=\$2; if (island_end-island_start>=50){count=count+island_end-island_start}}  if (NR==0){NR=1} printf(\"%s,%0.02f\n\", fname, count*100/NR)}' \$f >> ${WORK_DIR}/stats.csv"  >> "$WORK_DIR"/collect_stats.sh 
echo "fi;	done "  >> "$WORK_DIR"/collect_stats.sh

jid=`sbatch < "$WORK_DIR"/collect_stats.sh`

dependcollectstats="afterok:${jid##* }"

######################################################################
######################################################################
########## Produce contigs - can happen in parallel
######################################################################
######################################################################

jid=`sbatch <<- CONTIG | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH -p $QUEUE_X86
	#SBATCH -o ${LOG_DIR}/contig-%j.out
	#SBATCH -e ${LOG_DIR}/contig-%j.err
	#SBATCH -t 600
	#SBATCH -n 1 
	#SBATCH -c 1
	#SBATCH --mem-per-cpu=10G
	#SBATCH --threads-per-core=1
	$slurm_depend_merge

	$LOAD_MEGAHIT
	echo "Running $MEGAHIT_CMD -1 $read1filescomma -2 $read2filescomma -o ${WORK_DIR}/contigs"
	rm -rf ${WORK_DIR}/contigs
	$MEGAHIT_CMD -1 $read1filescomma -2 $read2filescomma -o ${WORK_DIR}/contigs
	mv ${WORK_DIR}/contigs/final.contigs.fa ${FINAL_DIR}/.
CONTIG`


echo "CONTIG_LENGTH=\$(tail -n2 ${WORK_DIR}/contigs/log |grep -o 'total.*' | cut -f2 -d' ')" > ${WORK_DIR}/call_dotplot.sh
echo "$PYTHON_CMD ${PIPELINE_DIR}/remove_strays.py ${WORK_DIR}/${MATCH_NAME}/aligned/depth_per_base.txt ${WORK_DIR}/${MATCH_NAME}/aligned/depth_per_base_without_primers_islands.txt" >> ${WORK_DIR}/call_dotplot.sh
echo "$PYTHON_CMD ${PIPELINE_DIR}/dot_coverage.py ${WORK_DIR}/${MATCH_NAME}/aligned/depth_per_base_without_primers_islands.txt ${FINAL_DIR}/contig.paf ${WORK_DIR}/stats.csv  \$CONTIG_LENGTH ${FINAL_DIR}/report " >> ${WORK_DIR}/call_dotplot.sh

# need to wait for match alignment
dependcontig="${dependmatchdone}:$jid"

######################################################################
######################################################################
########## Minimap2 - after contig and alignment
######################################################################
######################################################################
jid=`sbatch <<- MINIMAP | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH --partition=$QUEUE_X86 
	#SBATCH -o ${LOG_DIR}/pairwise-%j.out
	#SBATCH -e ${LOG_DIR}/pairwise-%j.err
	#SBATCH -t 600
	#SBATCH -n 1 
	#SBATCH -c 2
	#SBATCH --mem-per-cpu=2G
	#SBATCH --threads-per-core=1 
	#SBATCH -d $dependcontig
	
	$LOAD_MINIMAP2
	$LOAD_SAMTOOLS
	$MINIMAP2_CMD -x asm5 $MATCH_REF ${FINAL_DIR}/final.contigs.fa > ${FINAL_DIR}/contig.paf
	if [ ! -s ${FINAL_DIR}/contig.paf ]
	then
    		echo "!*** Pairwise alignment by minimap2 failed."
		touch ${FINAL_DIR}/contig.paf
	fi
MINIMAP`
dependcollectstats="${dependcollectstats}:$jid"

######################################################################
######################################################################
########## Dot plot - after contig and alignment
######################################################################
######################################################################

jid=`sbatch <<- DOTPLOT | egrep -o -e "\b[0-9]+$"
	#!/bin/bash -l
	#SBATCH --partition=$QUEUE_X86
	#SBATCH -o ${LOG_DIR}/dotplot-%j.out
	#SBATCH -e ${LOG_DIR}/dotplot-%j.err
	#SBATCH -t 600
	#SBATCH -n 1 
	#SBATCH -c 2
	#SBATCH --mem-per-cpu=2G
	#SBATCH --threads-per-core=1 
	#SBATCH -d $dependcollectstats
	
	$LOAD_PYTHON
	source ${WORK_DIR}/call_dotplot.sh
DOTPLOT`


echo "(-: Finished adding all jobs... Now is a good time to get that cup of coffee... Last job id $jid"
