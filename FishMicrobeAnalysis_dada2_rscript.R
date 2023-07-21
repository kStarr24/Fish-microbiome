# 2022-09-21
# Krista Starr
# Fish Microbiome DADA2 script to run on HCC crane
# This code is designed to be run in R/4.1 in Crane.

# Beging DAD2 pipeline
fastq_files = "/work/samodha/krista/fastqFiles_final2"
list.files(fastq_files)

mapping = read.csv("fishProject_meta.csv")
head(mapping)

# import may be an issue in batch jobs. try commenting out if the job fails
# library("import")
library("knitr")
library("BiocManager")
library("BiocStyle")
library("ggplot2")
library("gridExtra")
library("Rcpp")
library("dada2")
library("phyloseq")
library("DECIPHER")
library("ape")
library("phangorn")
library("ShortRead")

# list forward and reverse fastq files
fnFs = sort(list.files(fastq_files, pattern="_F_sampleNames.fastq"))
fnRs = sort(list.files(fastq_files, pattern="_R_sampleNames.fastq"))

# check file list
head(fnFs)
head(fnRs)

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)

#Specify full file path
fnFs <- file.path(fastq_files, fnFs)
fnRs <- file.path(fastq_files, fnRs)
fnFs[1:3]
fnRs[1:3]

# Plot quality profiles
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])

# Create a pathway to store filtered fastq files
filt_path <- file.path(fastq_files, "filtered") 
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

# Trimming and filtering the F/R reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE, matchIDs = T)
# view result
out
# save it!
save(out, file = "out.rds")

# Statistics after trimming
#total reads in
sum(out[,1]) 
#total reads out
sum(out[,2]) 
#reads lost
sum(out[,1]) - sum(out[,2])
# percentage data retained
sum(out[,2])/sum(out[,1]) 

#to avoid error, due to very low reads
exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames
# save
save(derepFs, derepRs, file = "dereps.rds")

# Learn the error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF)
plotErrors(errR)
# save
save(errF, errR, file = "error_rates.rds")

# Use the core dada2 inference method
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]] #summary of first sample
# save
save(dadaFs, dadaRs, file = "dada_function.rds")

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=T)
# save
save(mergers, file = "mergers.rds")

# Create a sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # returns #samples, #numASVs
save(seqtab, file = "seqtab.rds")

# Get a table with the distribution of sequence lengths
table = table(nchar(getSequences(seqtab)))
save(table, file = "seqLength_dist.rds")

# Remove chimera sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
save(seqtab.nochim, file = "seqtab.nochim.rds")

# Assign taxonomy
reference_file = paste0(fastq_files, "/", "silva_nr99_v138.1_wSpecies_train_set.fa.gz")
taxa <- assignTaxonomy(seqtab.nochim, reference_file , multithread=TRUE)
save(taxa, file = "taxa.rds")

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sampleNames
head(track)
write.csv(track, file = "track_reads.csv")


dir.create("Analysis")

# Give our sequence headers more manageble names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# Making and writing out a fasta of our final ASV seqs:
  asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "Analysis/ASVs.fa")

# Count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "Analysis/ASVs_counts.txt", sep="\t", quote=F)

# Tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "Analysis/ASVs_taxonomy.txt", sep="\t", quote=F)
write.table(taxa, "Analysis/taxonomy.txt", sep="\t", quote=F)


# # MOTHUR CODE
# module load mothur
# mkdir mothur
# cd mothur
# wget https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138.tgz
# tar -zxvf silva.nr_v138.tgz
# 
# # copy your .fa file into the mothur folder
# cp /work/samodha/krista/fastqFiles_final2/Analysis/ASVs.fa .
# 
# mothur
# 
# system(mv silva.nr_v138.align silva.v4.fasta)
# align.seqs(fasta=ASVs.fa, reference=silva.v4.fasta)
# 
# quit()
# 
# sed -i -e 's/>/>AAAAAAAAAA/g' ASVs.align
# sed -i -e 's/\./-/g' ASVs.align
# 
# module load mothur
# mothur
# 
# dist.seqs(fasta=ASVs.align, processors=16, cutoff=.10, output=phylip)
# clearcut(phylip=ASVs.phylip.dist)
# 
# quit()
# 
# sed -i -e 's/AAAAAAAAAA//g' ASVs.phylip.tre
# 
# # set to directory of your project!
# cp ASVs.phylip.tre /work/samodha/krista/fastqFiles_final2/Analysis/





