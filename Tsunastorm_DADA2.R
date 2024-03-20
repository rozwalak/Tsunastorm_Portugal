#################################################################################
# This is a DADA2 script to analyse sequencing data from the project Tsunastorm #
# The entire script is based on DADA2 Pipeline Tutorial:                        #
# https://benjjneb.github.io/dada2/tutorial.html                                #
#                                                                               # 
# Script is only an example and is slightly different than original analysis,   #
# performed using SLIM web application:                                         #
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2663-2#
#################################################################################

#DADA2 installation, here I installed software in R 4.3 version
#Skip if you have installed DADA2

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.18")

#Data preprocessing

library(dada2)

path <- "./data/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_Martinhal_V9_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_Martinhal_V9_R2.fastq", full.names = TRUE))
#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality profiles

#Visualise quality of forward reads
plotQualityProfile(fnFs[1:2])

#Visualise quality of reverse reads
plotQualityProfile(fnRs[1:2])

#Filter and trim 

#Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
#Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

#If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Code below is different than in the paper. Here taxonomy assignment is based on naive Bayes classifier incorporated in DADA2
#In publication we replaced it with more accurate vsearch, which is not possible to run in R studio easily without installing external software in linux

taxa <- assignTaxonomy(seqtab.nochim, "./pr2_version_4.14.0_SSU_dada2.fasta", multithread=TRUE, taxLevels = c("kingdom", "supergroup", "division", "class", "order", "family", "genus", "species"))
taxa

#Convert sequence table to a dataframe
seqtab_df <- as.data.frame(seqtab.nochim)
seqtab_df_transposed <- t(seqtab_df)

ASV_table <- merge(seqtab_df_transposed, taxa, by = "row.names", all.x = TRUE)
colnames(ASV_table )[1] <- "Sequence"

#final result:
ASV_table
