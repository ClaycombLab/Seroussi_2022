#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

################################################################################################
################################################################################################
################################################################################################
# Seroussi et. al. 2022
# This script finds overlaps between the C. elegans genome and sequencing data in BAM format
# This counting script uses the findOverlaps function from the
# GenomicAlignments package to assign reads to genomic features. Multiple
# aligning reads or reads that align to more than one feature are dealt with by
# counting reads in a sequential manner to the different gene biotypes in the following
# order (AS stands for antisense): miRNA, piRNA, rRNA, snoRNA, snRNA, tRNA,
# ncRNA, lincRNA, repeats AS, protein coding AS, pseudogene AS, lincRNA AS,
# antisense RNA, rRNA AS, snoRNA AS, snRNA AS, tRNA AS, ncRNA AS, miRNA
# AS, piRNA AS, protein coding, pseudogene, antisense RNA AS, repeats. 
# For reads that align to more than one feature in the same biotype group, the read count is
# split between the features based on the fraction of uniquely aligned reads to each of
# those features (unique weighing).
################################################################################################
################################################################################################
################################################################################################

# load libraries
suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(GenomicRanges)
  library(stringr)
  library(seqinr)
  require(data.table)
})

# read in references
load("ce11_WS276_miRbase_22.1_repeat_masker_dfam_2.0.RData")
load("ce11_WS276_miRbase_22.1_repeat_masker_dfam_2.0_parts.RData")

bam_file <- args[1]
name <- gsub(".*/([^/]+)/.*", "\\1", bam_file, perl=TRUE)

#read in entire BAM file
param <- ScanBamParam(what = c("seq","strand", "flag", "qname", "cigar", "rname", "pos"), reverseComplement=TRUE)
bam <- readGAlignments(bam_file, param=param)

total.count <- length(unique(mcols(bam)$qname))

perfect <- grepl("(^\\d+M$)|(^\\d+M\\d+N\\d+M$)", cigar(bam)) # perfect and splice junction matches only
bam.unperfect <- bam[!perfect]
bam <- bam[perfect]
bam.unperfect <- bam.unperfect[!(mcols(bam.unperfect)$qname %in% mcols(bam)$qname)]  # if perfect exists it takes precedence

bam.unperfect.best <- as.data.table(mcols(bam.unperfect))
bam.unperfect.best$width <- width(bam.unperfect)

bam.unperfect.best <- bam.unperfect[bam.unperfect.best[, .I[width == max(width)], by=qname]$V1] #get only the best matches out of all unperfect matches per read

bam = c(bam,bam.unperfect.best)
rm(perfect)
rm(bam.unperfect)
rm(bam.unperfect.best)
rm(param)

#find overlaps
ov <- findOverlaps(bam, gene.whole, type = "any", ignore.strand = TRUE)
ov.parts <- findOverlaps(bam, gene.parts, type = "any", ignore.strand = TRUE)

bam <- as.data.table(mcols(bam)) #save only the metadata columns
bam[, pos := paste(rname,pos,flag, sep = ":")]
bam$flag <- NULL
bam$rname <- NULL

ov <- as.data.table(ov)
ov[, reads := bam[ov[, queryHits]]$qname]
ov[, read_strand := bam[ov[, queryHits]]$strand]
ov[, biotype := mcols(gene.whole[ov[, subjectHits]])$gene_biotype]
ov[, biotype_strand := as.data.table(granges(gene.whole))[ov[, subjectHits]]$strand]

ov.parts <- as.data.table(ov.parts)
ov.parts[, reads := bam[ov.parts[, queryHits]]$qname]
ov.parts[, read_strand := bam[ov.parts[, queryHits]]$strand]
ov.parts[, biotype := mcols(gene.parts[ov.parts[, subjectHits]])$gene_biotype]
ov.parts[, biotype_strand := as.data.table(granges(gene.parts))[ov.parts[, subjectHits]]$strand]

colnames(ov)[2] <- "gene"
ov.parts$gene <- gene.whole[gene.parts[ov.parts$subjectHits]$whole.ref]$parts.ref

ov.parts$subjectHits <- gene.parts[ov.parts$subjectHits]$type
colnames(ov.parts)[2] <- "feature"

# little sanity check
stopifnot(sum(gene.parts$gene_id != gene.whole[gene.whole[gene.parts$whole.ref]$parts.ref]$gene_id) == 0)
rm(gene.parts) # remove gene.parts as not needed anymore

# update whether read is antisense to feature and remove no longer needed columns
ov[(read_strand == "-" & biotype_strand == "+") | (read_strand == "+" & biotype_strand == "-"), biotype := paste0(biotype,"_AS")]
ov.parts[(read_strand == "-" & biotype_strand == "+") | (read_strand == "+" & biotype_strand == "-"), biotype := paste0(biotype,"_AS")]
ov$read_strand <- NULL
ov.parts$read_strand <- NULL
ov$biotype_strand <- NULL
ov.parts$biotype_strand <- NULL

# merge whole and parts
ov <- merge(ov, ov.parts, by=c("gene", "queryHits", "reads", "biotype"))
rm(ov.parts)

# set gene names
ov[, gene := gene.whole[gene]$gene_id]
rm(gene.whole)

# sequential counting of biotypes
biotypes = c("miRNA","piRNA","rRNA","snoRNA","snRNA","tRNA","ncRNA","lincRNA",
             "dfam_repeats_AS","protein_coding_AS","pseudogene_AS","lincRNA_AS",
             "antisense_RNA","rRNA_AS","snoRNA_AS","snRNA_AS","tRNA_AS","ncRNA_AS","miRNA_AS","piRNA_AS",
             "protein_coding","pseudogene","antisense_RNA_AS","dfam_repeats")

ov.counts = list()

for (btype in biotypes){
  # determine which reads map to more than one feature
  ov[biotype == btype, shared := .N, by=reads]
  # determine distribution of unique reads
  ov[biotype == btype, prior := sum(shared == 1), by=list(gene, feature)]
  ov[biotype == btype, sum.prior := sum(prior), by=reads]
  
  ov.counts[[btype]] <- ov[biotype == btype]
  
  seen <- ov[biotype == btype]$reads
  ov <- ov[!(reads %in% seen) & is.na(shared), ]
  
  rm(seen)
}

ov <- rbindlist(ov.counts)
rm(ov.counts)

# partition reads between genes
ov[, fraction := prior / sum.prior ]
ov[is.nan(fraction), fraction := 1 / shared]  # can this be optimized
# remove unnecessary columns
ov$shared <- NULL
ov$prior <- NULL
ov$sum.prior <- NULL
# add sequence properties of reads
ov[, seq := as.character(bam[queryHits]$seq)]
ov[, cigar := as.character(bam[queryHits]$cigar)]
ov[, pos := as.character(bam[queryHits]$pos)]
# remove unnecessary columns
ov$queryHits <- NULL

# seen reads
seen <- unique(ov[, reads])

# what didn't map to a feature => deal with this first so can remove bam and then count those mapped to feature
no.features <- data.table(
  gene = "no_feature",
  biotype = "no_feature",
  feature = 0,
  seq = as.character(bam[!(bam$qname %in% seen)]$seq),
  cigar = as.character(bam[!(bam$qname %in% seen)]$cigar),
  pos = as.character(bam[!(bam$qname %in% seen)]$pos),
  count = 1 / as.numeric(table(bam[!(bam$qname %in% seen)]$qname)[bam[!(bam$qname %in% seen)]$qname])
)[,.(count = sum(count)),by=.(gene,biotype,feature,seq,cigar,pos)]

rm(bam)
rm(seen)

# count reads
ov <- ov[, .(count = sum(fraction)), by=.(gene,biotype,feature,seq,cigar,pos)]

# merge everything
ov <- rbind(ov, no.features)

# place in lists and save
results <- list()
results.sizes <- list()

results[[name]] <- ov
results.sizes[[name]] <- c(mapped=total.count)

cat("DONE!\n",
    "total mapped: ",total.count,"\n",
    "total counted: ",round(sum(ov$count)),"\n")

# create results directory and save
dir.create("sequential_count_results")
save(list=c("results", "results.sizes"), file=paste("sequential_count_results/", name, ".RData", sep=""))

