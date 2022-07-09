# QA and Counting of reads against the C. elegans genome as in Seroussi et. al. 2022
# Starting point: current working directory should hold fastq.gz files
# Run the commands below in order in terminal

# MD5 Checksum
# Check if any files were corrupted during downloading
ls *fastq.gz > fastq.txt
$(pwd)/scripts/checkMD5.sh fastq.txt

# Run FastQC - QC analysis tool,reads a set of sequence files and produces from each one a quality control report.
mkdir qc_fastq
fastqc --outdir ./qc_fastq --threads 12 *fastq.gz

# Cutadapt to remove adapters
# Also discard reads shorter than length=16 (-m), and discard reads longer than 30 (-M)
mkdir cutadapt
adapter3_nebnext=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
parallel -j 12 \
cutadapt -a $adapter3_nebnext -f fastq -m 16 -M 30 --discard-untrimmed -o cutadapt/{} {} ::: *.fastq.gz 

# Rerun fastQC post adapter removal
mkdir qc_cutadapt_fastq
fastqc --outdir qc_cutadapt_fastq --threads 12 cutadapt/*.fastq.gz

# Create a text file with all the samples that need to be processed
find $(pwd)/cutadapt/*.fastq > fastq.txt

# Generate indexes for your favourite genome reference as in https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf 
# Directory of where the STAR indices of your genome are
GENOME=$(pwd)/genome_data/STAR_indices/YOUR_STAR_INDICES_OF_SELECTED_GENOME

# modify this directory for your current working directory
MAIN_DIR=$(pwd)

# Run STAR
while read fastq; 
do
   base=`basename ${fastq}`					
   mkdir ${base%.fastq.gz} && cd $_			
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
