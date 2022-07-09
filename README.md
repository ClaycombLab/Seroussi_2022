# A Comprehensive Survey of C. elegans Argonaute Proteins Reveals Organism-wide Gene Regulatory Networks and Functions

This repository holds scripts to perform small RNA alignment and counting as described in Seroussi et. al. 2022.

## Summary of scripts

### The pipeline

If you wish to follow the exact same process as described in Seroussi et. al. 2022, you can follow the instruction in main_pipline.sh

```
# QA and Counting of reads against the C. elegans genome as in Seroussi et. al. 2022

# Starting point: current working directory should hold fastq.gz files

# Run the commands below in order in terminal

# MD5 Checksum

# Check if any files were corrupted during downloading

ls \*fastq.gz > fastq.txt
$(pwd)/scripts/checkMD5.sh fastq.txt

# Run FastQC - QC analysis tool,reads a set of sequence files and produces from each one a quality control report.

mkdir qc_fastq
fastqc --outdir ./qc_fastq --threads 12 \*fastq.gz

# Cutadapt to remove adapters

# Also discard reads shorter than length=16 (-m), and discard reads longer than 30 (-M)

mkdir cutadapt
adapter3_nebnext=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
parallel -j 12 \
cutadapt -a $adapter3_nebnext -f fastq -m 16 -M 30 --discard-untrimmed -o cutadapt/{} {} ::: \*.fastq.gz

# Rerun fastQC post adapter removal

mkdir qc_cutadapt_fastq
fastqc --outdir qc_cutadapt_fastq --threads 12 cutadapt/\*.fastq.gz

# Create a text file with all the samples that need to be processed

find $(pwd)/cutadapt/\*.fastq > fastq.txt

# Generate indexes for your favourite genome reference as in https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# Directory of where the STAR indices of your genome are

GENOME=$(pwd)/genome_data/STAR_indices/YOUR_STAR_INDICES_OF_SELECTED_GENOME

# modify this directory for your current working directory

MAIN_DIR=$(pwd)

# Run STAR

while read fastq;
do
base=`basename ${fastq}`
mkdir ${base%.fastq.gz} && cd $\_
STAR \
 --runThreadN 12 \
 --genomeDir ${GENOME} \
 --readFilesIn ${fastq} \
 --readFilesCommand zcat \
 --outSAMtype BAM SortedByCoordinate \
 --quantMode GeneCounts \
 --limitBAMsortRAM 30000000000 \
 --outFilterMultimapNmax 50 \
 --outFilterMultimapScoreRange 0 \
 --outFilterMismatchNoverLmax 0.05 \
 --outFilterMatchNmin 16 \
 --outFilterScoreMinOverLread 0 \
 --outFilterMatchNminOverLread 0 \
 --alignIntronMax 1 \
 --genomeLoad LoadAndKeep
cd ${MAIN_DIR}

done < fastq.txt

# Move resulting folders to a new folder called 'STAR'

mkdir STAR;

# Make sure all your samples are in the STAR folder then run the counting script

for bam in `find STAR | grep .bam$`;
do
echo $bam;
Rscript $(pwd)/scripts/counting_script.R ${bam}
done >> output.txt
```

### Counting Script

If you wish to only use the counting script on a given BAM file, run `> Rscript counting_script.R file.bam`.
The reference directory holds genome annotations of genes and their parts.
The WormBase version WS276 PRJNA13758 ce11 canonical geneset annotations were used (excluding miRNAs, repeats and transposons). C. elegans miRNA annotations were obtained from miRBase (release 22.1). For repeats and transposons RepeatMasker + Dfam (ce10 - Oct 2010 - RepeatMasker open-4.0.6 - Dfam 2.0) annotations were used. The UCSC Lift Genome Annotations tool was used to convert ce10 coordinates to ce11 coordinates.
