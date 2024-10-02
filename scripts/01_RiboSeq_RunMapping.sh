#!/bin/bash
### Footprints analysis - based on previous scripts available in the lab
### Virginie Ricci - October 2024

DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE

######################################################################
######################################################################
### Variables to define ### 
Prefix='RiboSeq'
Run_nb=2371
Nb_max_seq=1000000 # 0 >> no limit; default >> 10mio

echo 'Defined variables:'
echo "Genome = ${GENOME}"
echo "Prefix = ${Prefix}"
echo "Run_nb = ${Run_nb}"
echo "Nb_max_seq = ${Nb_max_seq} (if not specified, 10mio by default)"

# and all necessary paths
mkdir -p config/
mkdir -p raw_data/
mkdir -p mapping_data/


### Required files ###
# config/${Prefix}_pipeline_conf-riboprof_umi.sh
## config file with all necessary parameters

# {Prefix}_samples.txt in mapping_data/
## list of sample names (one column, no header)

# RiboSeq .fastq files in raw_data/


######################################################################
######################################################################
### This script run the mapping and QC analysis on samples that were
### generated with some modifications on the rRNA purification
### steps. They are a test to check if the protocol works. Moreover
### this is a test run with very low throughput

# Add the modules:
module add samtools/1.19.2;
module add fastx_toolkit/0.0.14;
module add bowtie2/2.3.5;
module add STAR/2.7.11b;
module add cutadapt/3.5;
module add chip-seq/1.5.5; 
module add bedtools/2.28.0; 
module add UMI-tools/1.0.0git;
module add pipeline/0.3.2;
module add rpfTools/1.2.0;


######################################################################
echo "pipeline_get_barcodes.sh"
pipeline_get_barcodes.sh -s ${Nb_max_seq} \
			 -c config/${Prefix}_pipeline_conf-riboprof_umi_STARonly_sequentialmapping.sh \
			 -r raw_data/ \
			 ${Prefix}_samples.txt



echo "meta_pipeline_STAR_sequentialmapping.sh"
echo "sequential mapping to mouse rRNA, human rRNA, mouse tRNA, and mouse genome"
meta_pipeline_STAR_sequentialmapping.sh -o mapping_data/ \
		-d ${Prefix}_samples.db \
		-f config/${Prefix}_pipeline_conf-riboprof_umi_STARonly_sequentialmapping.sh \
		-s _${Prefix} \
		${Prefix}_samples.txt



echo 'Dedup and split .bam file of reads mapped to genome (not part of the pipeline)'
for infile in $(ls *genome.Aligned.out.bam); do
    filebase1=$(basename $infile .genome.Aligned.out.bam)
    filebase2=genome
    logfile=umi_tools.log
    
    umi_tools dedup -I $infile --method=directional --per-cell --read-length --log=$logfile | samtools view -h - | awk -v filebase1=$filebase1 -v filebase2=$filebase2'{split($1,barcode,"_"); print $0 > filebase1 "_split_" barcode[2] filebase2 ".sam"}'
done

for infile in $(ls *_split_*.sam); do
    outfile=$(basename .sam).bam

    samtools view -bS $infile > $outfile
    samtools index $outfile
done


echo "pipelineQC_STAR on reads mapped to genome"
pipelineQC_STAR -i mapping_data/ \
	   -f config/${Prefix}_pipeline_conf-riboprof_umi_STARonly_sequentialmapping.sh \
	   -s _${Prefix} \
	   -a ${Prefix}_samples.txt

Rscript PlotQC.R



######################################################################
######################################################################

DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE