#!/bin/bash
#Pipeline for whole genome RSV
#March 2019
#Pavitra Roychoudhury

#### Load required modules ####
#**versions are current as of 5-Mar-19**
# cd /fh/fast/jerome_k/RSV_WGS/
# module load bowtie2/2.2.5
# module load FastQC/0.11.8-Java-1.8
# module load R/3.5.1-foss-2016b-fh1 #3.5.2 has no Rsamtools
# module load prokka/1.13-foss-2016b-BioPerl-1.7.0
# module load parallel/20170222-foss-2016b
# module load Perl/5.28.0-foss-2016b
# module load SAMtools/1.8-foss-2016b #1.9 has errors involving cURL
# module load BBMap/36.62-foss-2016b-Java-1.8.0_121
# module load LAST/926-foss-2016b
# wget https://github.com/tseemann/prokka/raw/master/binaries/linux/tbl2asn -O $HOME/.local/bin/tbl2asn
# curl -o ./rsv_make_reference.R https://raw.githubusercontent.com/proychou/ViralWGS/master/make_ref_from_assembly.R

#### One time: First build reference for bowtie and make a copy of the ref seqs ####
# module load bowtie2
# bowtie2-build './refs/RSVA_KU950499.fasta' ./refs/rsvA_ref
# bowtie2-build './refs/RSVB_KU950676.fasta' ./refs/rsvB_ref
# cp './refs/RSVA_KU950499.fasta' ./refs/rsvA_ref.fasta
# cp './refs/RSVB_KU950676.fasta' ./refs/rsvB_ref.fasta
# cat ./refs/rsvA_ref.fasta ./refs/rsvB_ref.fasta > ./refs/rsv_refs.fasta
# prokka-genbank_to_fasta_db ./refs/RSVA_KU950499.gb ./refs/RSVB_KU950676.gb > ./refs/rsv_proteins.faa
# bwa index ./refs/rsvA_ref.fasta
# bwa index ./refs/rsvB_ref.fasta

#### Usage ####
#For paired-end library
#		hsv1_pipeline.sh -1 yourreads_r1.fastq.gz -2 yourreads_r2.fastq.gz
#For single-end library
#		hsv1_pipeline.sh -s yourreads.fastq.gz
#This is meant to be run on the cluster (typically through sbatch) so if run locally,
#first set the environment variable manually, e.g.
#		SLURM_CPUS_PER_TASK=8
#or whatever is the number of available processors 
 
#Load required tools
#Note that spades and last are all locally installed and need to be updated manually as required

#To do: 
# - add a restart option
# - move some of the common parts of this script (between HHV6 and HSV) to ViralWGS 


#Test run
# in_fastq_r1='/fh/fast/jerome_k/SR/ngs/illumina/proychou/190222_D00300_0682_BHTMLLBCX2/Unaligned/Project_jboonyar/Sample_GH100084/GH100084_GGAGATTC-GTTCAGAC_L002_R1_001.fastq.gz'
# in_fastq_r2='/fh/fast/jerome_k/SR/ngs/illumina/proychou/190222_D00300_0682_BHTMLLBCX2/Unaligned/Project_jboonyar/Sample_GH100084/GH100084_GGAGATTC-GTTCAGAC_L002_R2_001.fastq.gz'
# SLURM_CPUS_PER_TASK=8


PATH=$PATH:$HOME/.local/bin:$HOME/SPAdes-3.11.1-Linux/bin:
export PATH=$PATH:$EBROOTPROKKA/bin:$EBROOTPROKKA/db:
echo "Number of cores used: "$SLURM_CPUS_PER_TASK
# echo "Path: "$PATH

while getopts ":1:2:s:f" opt; do
	case $opt in
		1) in_fastq_r1="$OPTARG"
			paired="true"
		;;
		2) in_fastq_r2="$OPTARG"
			paired="true"
		;;
		s) in_fastq="$OPTARG"
			paired="false"
		;;
		f) filter="true"
		;;
		\?) echo "Invalid option -$OPTARG" >&2
    	;;
	esac
done

printf "Input arguments:\n\n"
echo $@


##  PAIRED-END  ##
if [[ $paired == "true" ]]
then
if [ -z $in_fastq_r1 ] || [ -z $in_fastq_r2 ]
then
echo "Missing input argument."
fi

sampname=$(basename ${in_fastq_r1%%_R1_001.fastq*})

#FastQC report on raw reads
printf "\n\nFastQC report on raw reads ... \n\n\n"
mkdir -p ./fastqc_reports_raw
fastqc $in_fastq_r1 $in_fastq_r2 -o ./fastqc_reports_raw -t $SLURM_CPUS_PER_TASK 

#Remove optical duplicates -- TO DO 
# printf "\n\nRemove optical duplicates from raw reads ... \n\n\n"
# clumpify.sh in=reads.fq.gz out=clumped.fq.gz dedupe optical dist=40 spantiles=f


#Adapter trimming with bbduk
printf "\n\nAdapter trimming ... \n\n\n"
mkdir -p ./trimmed_fastq
bbduk.sh in1=$in_fastq_r1 in2=$in_fastq_r2  out1='./trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' out2='./trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz' ref=~/bbmap/resources/adapters.fa k=21 ktrim=r mink=4 hdist=2 tpe tbo overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
bbduk.sh in1='./trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' in2='./trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz'  out1='./trimmed_fastq/'$sampname'_trimmed_r1.fastq.gz' out2='./trimmed_fastq/'$sampname'_trimmed_r2.fastq.gz' ref=~/bbmap/resources/adapters.fa k=21 ktrim=l mink=4 hdist=2 tpe tbo overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
rm './trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' './trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz'

#Quality trimming
printf "\n\nQuality trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
bbduk.sh in1='./trimmed_fastq/'$sampname'_trimmed_r1.fastq.gz' in2='./trimmed_fastq/'$sampname'_trimmed_r2.fastq.gz' out1='./preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' out2='./preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' t=$SLURM_CPUS_PER_TASK qtrim=rl trimq=20 maq=10 overwrite=TRUE minlen=20

#Use bbduk to filter reads that match RSV genomes -- ** not tested, needs edit **
# if [[ $filter == "true" ]]
# then
# printf "\n\nK-mer filtering using hsv_refs.fasta ... \n\n\n"
# mkdir -p ./filtered_fastq/
# 
# bbduk.sh in1='./preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' in2='./preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' out1='./filtered_fastq/'$sampname'_unmatched_r1.fastq.gz' out2='./filtered_fastq/'$sampname'_unmatched_r2.fastq.gz' outm1='./filtered_fastq/'$sampname'_matched_r1.fastq.gz' outm2='./filtered_fastq/'$sampname'_matched_r2.fastq.gz' ref='./refs/hsv_refs.fasta' k=31 hdist=2 stats='./filtered_fastq/'$sampname'_stats_hsv.txt' overwrite=TRUE t=$SLURM_CPUS_PER_TASK
# 
# rm './filtered_fastq/'$sampname'_unmatched_r1.fastq.gz' './filtered_fastq/'$sampname'_unmatched_r2.fastq.gz' 
# mv './filtered_fastq/'$sampname'_matched_r1.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz'
# mv './filtered_fastq/'$sampname'_matched_r2.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz'
# fi

#FastQC report on processed reads
mkdir -p ./fastqc_reports_preprocessed
printf "\n\nFastQC report on preprocessed reads ... \n\n\n"
fastqc './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' -o ./fastqc_reports_preprocessed -t $SLURM_CPUS_PER_TASK 


#Map reads to reference
printf "\n\nMapping reads to reference seqs rsvA_ref, rsvB_ref ... \n\n\n"
mkdir -p ./mapped_reads
for ref in rsvA_ref rsvB_ref; do
bowtie2 -x ./refs/$ref -1 './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' -2 './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' -p ${SLURM_CPUS_PER_TASK} -S './mapped_reads/'$sampname'_'$ref'.sam'
done

#Assemble with SPAdes 
printf "\n\nStarting de novo assembly ... \n\n\n"
mkdir -p './contigs/'$sampname
spades.py -1 './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' -2 './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' -o './contigs/'$sampname --careful -t ${SLURM_CPUS_PER_TASK}

#Delete some spades folders to free up space
rm -r './contigs/'$sampname'/corrected' 


#Filter contigs by length, throw out any shorter than 200bp -- already being done within R script, plus filtering by coverage so commenting this out
# printf "\n\nFiltering contigs by length ... \n\n\n"
# reformat.sh in='./contigs/'$sampname'/scaffolds.fasta' out='./contigs/'$sampname'/scaffolds_filtered.fasta' minlength=200 overwrite=TRUE


##  SINGLE-END  ## ** not tested, needs edits **
else 
if [[ $paired == "false" ]]
then

printf "Single-end runs not tested yet. Exiting."

# if [ -z $in_fastq ]
# then
# echo "Missing input argument."
# fi
# 
# sampname=$(basename ${in_fastq%%_R1_001.fastq*})
# 
# #FastQC report on raw reads
# printf "\n\nFastQC report on raw reads ... \n\n\n"
# mkdir -p ./fastqc_reports_raw
# fastqc -o ./fastqc_reports_raw -t $SLURM_CPUS_PER_TASK $in_fastq 
# 
# #Adapter trimming with bbduk
# printf "\n\nAdapter trimming ... \n\n\n"
# mkdir -p ./trimmed_fastq
# bbduk.sh in=$in_fastq out='./trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz' ref=~/bbmap/resources/adapters.fa k=21 ktrim=r mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
# bbduk.sh in='./trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz'  out='./trimmed_fastq/'$sampname'_trimmed.fastq.gz' ref=~/bbmap/resources/adapters.fa k=21 ktrim=l mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
# rm './trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz'
# 
# #Quality trimming
# printf "\n\nQuality trimming ... \n\n\n"
# mkdir -p ./preprocessed_fastq
# bbduk.sh in='./trimmed_fastq/'$sampname'_trimmed.fastq.gz' out='./preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' t=$SLURM_CPUS_PER_TASK qtrim=rl trimq=20 maq=10 overwrite=TRUE minlen=20
# 
# #Use bbduk to filter reads that match HSV genomes
# if [[ $filter == "true" ]]
# then
# printf "\n\nK-mer filtering using hsv_refs.fasta ... \n\n\n"
# mkdir -p ./filtered_fastq/
# bbduk.sh in='./preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' out='./filtered_fastq/'$sampname'_unmatched.fastq.gz' outm='./filtered_fastq/'$sampname'_matched.fastq.gz' ref='./refs/hsv_refs.fasta' k=31 hdist=2 stats='./filtered_fastq/'$sampname'_stats_hsv.txt' overwrite=TRUE t=$SLURM_CPUS_PER_TASK
# rm './filtered_fastq/'$sampname'_unmatched.fastq.gz' 
# mv './filtered_fastq/'$sampname'_matched.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz'
# fi
# 
# #FastQC report on processed reads
# printf "\n\nFastQC report on preprocessed reads ... \n\n\n"
# mkdir -p ./fastqc_reports_preprocessed
# fastqc -o ./fastqc_reports_preprocessed -t $SLURM_CPUS_PER_TASK './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' 
# 
# 
# #Map reads to reference
# printf "\n\nMapping reads to reference seqs hsv1_ref, hsv2_ref_hg52 and hsv2_sd90e ... \n\n\n"
# mkdir -p ./mapped_reads
# for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
# bowtie2 -x ./refs/$ref -U './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -p ${SLURM_CPUS_PER_TASK} -S './mapped_reads/'$sampname'_'$ref'.sam'
# done
# 
# #Assemble with SPAdes
# printf "\n\nStarting de novo assembly ... \n\n\n"
# mkdir -p './contigs/'$sampname
# spades.py -s './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -o './contigs/'$sampname --careful -t ${SLURM_CPUS_PER_TASK}
# 
fi
fi



#Generate sorted bams for mapped reads (this bam is used for typing: A vs B)
printf "\n\nMaking and sorting bam files ... \n\n\n"
for ref in rsvA_ref rsvB_ref; do
if [ -f './mapped_reads/'$sampname'_'$ref'.sam' ]
then
samtools view -bh -o './mapped_reads/'$sampname'_'$ref'.bam' './mapped_reads/'$sampname'_'$ref'.sam' -T ./refs/$ref'.fasta'  
rm './mapped_reads/'$sampname'_'$ref'.sam'
samtools sort -o './mapped_reads/'$sampname'_'$ref'.sorted.bam' './mapped_reads/'$sampname'_'$ref'.bam' 
rm './mapped_reads/'$sampname'_'$ref'.bam' 
else
echo 'Mapping to '$ref 'failed. No sam file found'
fi
done


# Now call an R script that merges assembly and mapping and ultimately makes the consensus sequence assuming either A or B as reference
Rscript --vanilla rsv_make_seq.R sampname=\"$sampname\"


#Remap reads to "new" reference
printf "\n\nRe-mapping reads to assembled sequence ... \n\n\n"
mkdir -p ./remapped_reads

# if [[ $paired == "false" ]] # ** single-end not tested **
# then
# for ref in rsvA_ref rsvB_ref; do
# bowtie2-build './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref'_consensus.fasta' './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref
# bowtie2 -x './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref -U './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -p ${SLURM_CPUS_PER_TASK} -S './remapped_reads/'$sampname'_'$ref'.sam'
# done
# fi
 
if [[ $paired == "true" ]]
then
for ref in rsvA_ref rsvB_ref; do
bowtie2-build -q './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref'_consensus.fasta' './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref
bowtie2 -x './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref -1 './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' -2 './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' -p ${SLURM_CPUS_PER_TASK} -S './remapped_reads/'$sampname'_'$ref'.sam'
done
fi

#Make sorted bam
for ref in rsvA_ref rsvB_ref; do
if [ -f './remapped_reads/'$sampname'_'$ref'.sam' ]
then
samtools view -bh -o './remapped_reads/'$sampname'_'$ref'.bam' './remapped_reads/'$sampname'_'$ref'.sam' -T './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref'_consensus.fasta'
rm './remapped_reads/'$sampname'_'$ref'.sam'
samtools sort -o './remapped_reads/'$sampname'_'$ref'.sorted.bam' './remapped_reads/'$sampname'_'$ref'.bam'
rm './remapped_reads/'$sampname'_'$ref'.bam'
mv './remapped_reads/'$sampname'_'$ref'.sorted.bam'  './remapped_reads/'$sampname'_'$ref'.bam' 
else
echo 'No sam file found'
fi
done

#Call R script to generate a consensus sequence
printf "\n\nGenerating consensus sequence ... \n\n\n"
mkdir -p ./consensus_seqs_all
mkdir -p ./stats
if [[ $paired == "true" ]]
then
Rscript --vanilla rsv_generate_consensus.R s1=\"$in_fastq_r1\"
else
if [[ $paired == "false" ]]
then
Rscript --vanilla rsv_generate_consensus.R s1=\"$in_fastq\"
fi
fi

#Annotate
printf "\n\nAnnotating with prokka ... \n\n\n"
mkdir -p ./annotations_prokka_rsvA
prokka --outdir './annotations_prokka_rsvA/'$sampname'/' --force --kingdom 'Viruses' --genus 'Human orthopneumovirus' --species '' --proteins ./refs/rsv_proteins.faa --locustag '' --strain $sampname --prefix $sampname --gcode 1 --evalue 1e-9 './annotations_prokka_rsvA/'$sampname/*.fa
mkdir -p ./annotations_prokka_rsvB
prokka --outdir './annotations_prokka_rsvB/'$sampname'/' --force --kingdom 'Viruses' --genus 'Human orthopneumovirus' --species '' --proteins ./refs/rsv_proteins.faa --locustag '' --strain $sampname --prefix $sampname --gcode 1 --evalue 1e-9 './annotations_prokka_rsvB/'$sampname/*.fa

#Clean up some files
rm './ref_for_remapping/'$sampname*'.fai'
rm './ref_for_remapping/'$sampname*'.bt2'
