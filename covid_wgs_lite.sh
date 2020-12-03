#!/bin/bash
#Pipeline for whole genome sequence assembly for COVID19_WGS
#This is the "lite" version that doesn't use de novo assembly
#Nov 2020
#Pavitra Roychoudhury

PATH=$PATH:$HOME/.local/bin:
echo "Number of cores used: "$SLURM_CPUS_PER_TASK

while getopts ":1:2:u:s:ap:q" opt; do
	case $opt in
		1) in_fastq_r1="$OPTARG"
			paired="true"
		;;
		2) in_fastq_r2="$OPTARG"
			paired="true"
		;;
		u) in_fastq="$OPTARG"
			paired="false"
		;;
		a) adapter_trim="true"
		;;
		q) qual_trim="true"
		;;
		p) primer_trim="true"
			primerlist="$OPTARG"
		;;
		s) sampname="$OPTARG"
		;;
		\?) echo "Invalid option -$OPTARG" >&2
    	;;
	esac
done

printf "Input arguments:\n\n"
echo $@

ref_fasta='./refs/NC_045512.2.fasta'
ref_bowtie='NC_045512.2'

if [ -z $sampname ]
then
echo "Missing sample name."
exit 1
fi

##  PAIRED-END  ##

if [[ $paired == "true" ]]
then


if [ -z $in_fastq_r1 ] || [ -z $in_fastq_r2 ]
then
echo "Missing fastq name(s)."
exit 1
fi


#FastQC report on raw reads
printf "\n\nFastQC report on raw reads ... \n\n\n"
mkdir -p ./fastqc_reports_raw
fastqc $in_fastq_r1 $in_fastq_r2 -o ./fastqc_reports_raw -t $SLURM_CPUS_PER_TASK  


#Adapter trimming with bbduk
if [[ $adapter_trim == "true" ]]
then
printf "\n\nAdapter trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
tmp_fastq1='./preprocessed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz'
tmp_fastq2='./preprocessed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz'
processed_fastq1='./preprocessed_fastq/'$sampname'_trimmed_r1.fastq.gz'
processed_fastq2='./preprocessed_fastq/'$sampname'_trimmed_r2.fastq.gz'

bbduk.sh in1=$in_fastq_r1 in2=$in_fastq_r2 out1=$tmp_fastq1 out2=$tmp_fastq2 ref=adapters,artifacts k=21 ktrim=r mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 

bbduk.sh in1=$tmp_fastq1 in2=$tmp_fastq2 out1=$processed_fastq1 out2=$processed_fastq2 ref=adapters,artifacts k=21 ktrim=l mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
rm $tmp_fastq1 $tmp_fastq2

else
processed_fastq1=$in_fastq_r1 
processed_fastq2=$in_fastq_r2
fi


#Quality trimming
if [[ $qual_trim == "true" ]]
then
printf "\n\nQuality trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
processed_fastq_old1=$processed_fastq1
processed_fastq_old2=$processed_fastq2
processed_fastq1='./preprocessed_fastq/'$sampname'_preprocessed_r1.fastq.gz'
processed_fastq2='./preprocessed_fastq/'$sampname'_preprocessed_r2.fastq.gz'

bbduk.sh in1=$processed_fastq_old1 in2=$processed_fastq_old2 out1=$processed_fastq1 out2=$processed_fastq2 t=$SLURM_CPUS_PER_TASK qtrim=rl trimq=20 maq=10 overwrite=TRUE minlen=20

fi


#Primer trimming: settings based on discussions in SPHERES consortium-- this assumes longer reads, so setting minimum read length to 75 if using primer trimming
if [[ $primer_trim == "true" ]]
then
printf "\n\nPrimer trimming ... \n\n\n"

if [ -z $primerlist ]
then
printf "Missing primer list!"
exit
fi

mkdir -p ./preprocessed_fastq
tmp_fastq1=$processed_fastq1
tmp_fastq2=$processed_fastq2
processed_fastq1='./preprocessed_fastq/'$sampname'_trimmed2_r1.fastq.gz'
processed_fastq2='./preprocessed_fastq/'$sampname'_trimmed2_r2.fastq.gz'

bbduk.sh in1=$tmp_fastq1 in2=$tmp_fastq2 out1=$processed_fastq1 out2=$processed_fastq2  ref=$primerlist k=18 ktrim=l hdist=3 qhdist=1 mink=4 rcomp=f overwrite=TRUE restrictleft=30 t=$SLURM_CPUS_PER_TASK minlen=75

tmp_fastq1=$processed_fastq1
tmp_fastq2=$processed_fastq2
processed_fastq1='./preprocessed_fastq/'$sampname'_trimmed3_r1.fastq.gz'
processed_fastq2='./preprocessed_fastq/'$sampname'_trimmed3_r2.fastq.gz'
bbduk.sh in1=$tmp_fastq1 in2=$tmp_fastq2 out1=$processed_fastq1 out2=$processed_fastq2 ref=$primerlist k=18 ktrim=r hdist=3 qhdist=1 mink=4 rcomp=f overwrite=TRUE restrictright=30 t=$SLURM_CPUS_PER_TASK minlen=75
rm $tmp_fastq1 $tmp_fastq2

fi




#Map reads to reference
printf "\n\nMapping reads to reference ... \n\n\n"
mkdir -p ./mapped_reads
mappedtoref_bam='./mapped_reads/'$sampname'.bam'
bowtie2 -x ./refs/$ref_bowtie -1 $processed_fastq1 -2 $processed_fastq2 -p ${SLURM_CPUS_PER_TASK} | samtools view -bS -F 4 - > $mappedtoref_bam
samtools sort -@ ${SLURM_CPUS_PER_TASK} -o './mapped_reads/'$sampname'.sorted.bam' $mappedtoref_bam 
rm $mappedtoref_bam 
mv './mapped_reads/'$sampname'.sorted.bam' $mappedtoref_bam 



##  SINGLE-END  ## 
else 
if [[ $paired == "false" ]]
then

if [ -z $in_fastq ]
then
echo "Missing fastq name."
exit 1
fi

#FastQC report on raw reads
printf "\n\nFastQC report on raw reads ... \n\n\n"
mkdir -p ./fastqc_reports_raw
fastqc -o ./fastqc_reports_raw -t $SLURM_CPUS_PER_TASK $in_fastq 

#Adapter trimming with bbduk
if [[ $adapter_trim == "true" ]]
then
printf "\n\nAdapter trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
tmp_fastq='./preprocessed_fastq/'$sampname'_trimmed_tmp.fastq.gz'
processed_fastq='./preprocessed_fastq/'$sampname'_trimmed.fastq.gz'

bbduk.sh in=$in_fastq out=$tmp_fastq ref=adapters,artifacts k=21 ktrim=r mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
bbduk.sh in=$tmp_fastq out=$processed_fastq ref=adapters,artifacts k=21 ktrim=l mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
rm $tmp_fastq

else
processed_fastq=$in_fastq 
fi

  
#Quality trimming
if [[ $qual_trim == "true" ]]
then
printf "\n\nQuality trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
processed_fastq_old=$processed_fastq
processed_fastq='./preprocessed_fastq/'$sampname'_preprocessed.fastq.gz'

bbduk.sh in=$processed_fastq_old out=$processed_fastq t=$SLURM_CPUS_PER_TASK qtrim=rl trimq=20 maq=10 overwrite=TRUE minlen=20

fi


#Primer trimming -- this assumes longer reads, so setting minimum read length to 75 if using primer trimming
if [[ $primer_trim == "true" ]]
then
printf "\n\nPrimer trimming ... \n\n\n"

if [ -z $primerlist ]
then
printf "Missing primer list!"
exit
fi

mkdir -p ./preprocessed_fastq
tmp_fastq=$processed_fastq
processed_fastq='./preprocessed_fastq/'$sampname'_trimmed2.fastq.gz'

bbduk.sh in=$tmp_fastq out=$processed_fastq ref=$primerlist k=18 ktrim=l hdist=3 qhdist=1 mink=4 rcomp=f overwrite=TRUE restrictleft=30 t=$SLURM_CPUS_PER_TASK minlen=75

tmp_fastq=$processed_fastq
processed_fastq='./preprocessed_fastq/'$sampname'_trimmed3.fastq.gz'
bbduk.sh in=$tmp_fastq out=$processed_fastq ref=$primerlist k=18 ktrim=r hdist=3 qhdist=1 mink=4 rcomp=f overwrite=TRUE restrictright=30 t=$SLURM_CPUS_PER_TASK minlen=75
rm $tmp_fastq

fi


#Map reads to reference
printf "\n\nMapping reads to reference ... \n\n\n"
mkdir -p ./mapped_reads
mappedtoref_bam='./mapped_reads/'$sampname'.bam'
bowtie2 -x ./refs/$ref_bowtie -U $processed_fastq -p ${SLURM_CPUS_PER_TASK} | samtools view -bS -F 4 - > $mappedtoref_bam
samtools sort -@ ${SLURM_CPUS_PER_TASK} -o './mapped_reads/'$sampname'.sorted.bam' $mappedtoref_bam 
rm $mappedtoref_bam 
mv './mapped_reads/'$sampname'.sorted.bam' $mappedtoref_bam 



fi
fi


# R script to make the consensus sequence 
Rscript --vanilla hcov_make_seq_lite.R sampname=\"$sampname\" bamfname=\"$mappedtoref_bam\" 


printf "\n\n Done. Celebrate! \n\n\n"