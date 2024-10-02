#!/bin/bash
### RNASeq analysis
### Virginie Ricci - October 2024

DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'Start of the job' $DATE

######################################################################
######################################################################
### Variables to define ### 
Prefix='RNASeq'
Run_nb=2327

echo 'Defined variables:'
echo "Genome = ${GENOME}"
echo "Prefix = ${Prefix}"
echo "Run_nb = ${Run_nb}"

# and all necessary paths
mkdir -p raw_data/
mkdir -p mapping_data/

N_CORE=''
rRNA_DB=''
H_rRNA_DB=''
tRNA_DB=''
GENOME_DB=''



### Required files ###
# RNASeq .fastq files in raw_data/

######################################################################
######################################################################
# Add the modules:
module add samtools/1.9;
module add bowtie2/2.3.5;
module add cutadapt/3.5;
module load cufflinks/2.2.1
module load python3/3.7.3
module add pipeline/0.3.2_VR
module add rpfTools/1.2.0_VR;
module add STAR/2.7.11b


######################################################################
# Remove first 2 nt that are bad quality
for infile in $(ls *.fq); do
        outfile=$(basename $infile .fq).trimmed.fq

        cutadapt -u 2 -o $outfile $infile
done



######################################################################
# STAR mapping

# Mouse rRNA - Step 1
echo 'STAR mapping to Mus musculus rRNA'
for infile in $(ls *.trimmed.fq); do
        sample=$(basename $infile .trimmed.fq)
        echo $sample

        STAR \
        --runThreadN $N_CORE \
        --runMode alignReads \
        --genomeDir $rRNA_DB \
        --readFilesIn $infile \
        --outFileNamePrefix ${sample}.STAR_1. \
        --outSAMmultNmax 1 \
        --seedSearchStartLmax 50 \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate Unsorted \
        --limitBAMsortRAM 15000000000 \
        --outReadsUnmapped Fastx \
        --outTmpDir STAR_1/${sample} \
        >(tee /dev/null) -S /dev/null - 2> ${sample}.STAR_1.log
done



# Human rRNA - Step 2
echo 'STAR mapping to Human rRNA'
for infile in $(ls *.trimmed.fq); do
        sample=$(basename $infile .trimmed.fq)
        echo $sample

        STAR \
        --runThreadN $N_CORE \
        --runMode alignReads \
        --genomeDir $H_rRNA_DB \
        --readFilesIn ${sample}.STAR_1.Unmapped.out.mate1 \
        --outFileNamePrefix ${sample}.STAR_2. \
        --outSAMmultNmax 1 \
        --seedSearchStartLmax 50 \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate Unsorted \
        --limitBAMsortRAM 15000000000 \
        --outReadsUnmapped Fastx \
        --outTmpDir STAR_2/${sample} \
        >(tee /dev/null) -S /dev/null - 2> ${sample}.STAR_2.log
done



# Mouse tRNA - Step 3
echo 'STAR mapping to Mus musculus tRNA'
for infile in $(ls *.trimmed.fq); do
        sample=$(basename $infile .trimmed.fq)
        echo $sample

        STAR \
        --runThreadN $N_CORE \
        --runMode alignReads \
        --genomeDir $tRNA_DB \
        --readFilesIn ${sample}.STAR_2.Unmapped.out.mate1 \
        --outFileNamePrefix ${sample}.STAR_3. \
        --outSAMmultNmax 1 \
        --seedSearchStartLmax 50 \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate Unsorted \
        --limitBAMsortRAM 15000000000 \
        --outReadsUnmapped Fastx \
        --outTmpDir STAR_3/${sample} \
        >(tee /dev/null) -S /dev/null - 2> ${sample}.STAR_3.log
done



# Mouse genome - Step 4
echo 'STAR mapping to Mus musculus genome'
for infile in $(ls *.trimmed.fq); do
        sample=$(basename $infile .trimmed.fq)
        echo $sample

        STAR \
        --runThreadN $N_CORE \
        --runMode alignReads \
        --genomeDir $GENOME_DB \
        --readFilesIn ${sample}.STAR_3.Unmapped.out.mate1 \
        --outFileNamePrefix ${sample}. \
        --outSAMmultNmax 1 \
        --seedSearchStartLmax 50 \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate Unsorted \
        --limitBAMsortRAM 15000000000 \
        --outReadsUnmapped Fastx \
        --outTmpDir ${sample} \
        >(tee /dev/null) -S /dev/null - 2> ${sample}.log


        samtools view -bS -F 4,8,256 ${sample}.Aligned.sortedByCoord.out.bam > ${sample}.Aligned.sortedByCoord.out.final.bam

        samtools index ${sample}.Aligned.sortedByCoord.out.final.bam
done



######################################################################
######################################################################

DATE=`date '+%d-%m-%Y %H:%M:%S'`
echo 'End of the job' $DATE