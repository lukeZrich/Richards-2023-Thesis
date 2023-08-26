rm(list = ls())
library(ShortRead)
library(dada2)
#Set working directory
setwd("~/ArabidopsisWheat/Fungi/")
#Set working directory for trimmed reads.
rawDir <- "./Raw"

#check files list and make lists of the forward sequences ONLY
# get sample names
list.files(rawDir)
fnFs <- sort(list.files(rawDir, pattern="_L001_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(rawDir, pattern="_L001_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "-110220"), `[`, 1)

#Filter and trim
#create filepaths for results
filtFs <- file.path("./Filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("./Filtered", paste0(sample.names, "_R_filt.fastq.gz"))
#trim and filter
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(25, 25),
                     maxN=0, maxEE = 2, rm.phix=TRUE, truncLen = c(280, 250), 
                     compress=TRUE, multithread=10, verbose = TRUE)
perc.reads.retained <- (out[,2] / out[,1]) * 100
out <- cbind(out, perc.reads.retained)
#write out a summary of filtering results
save(out, file = "./Filtered/out.R")


#learn error rates for algorithm and plot
errF <- learnErrors(filtFs, multithread=10, verbose = 1)
errR <- learnErrors(filtRs, multithread=10, verbose = 1)
errFplot <- plotErrors(errF)
errRplot <- plotErrors(errR)
dir.create("./ErrorPlots")
save(errFplot, file = "./ErrorPlots/errFplot.R")
save(errFplot, file = "./ErrorPlots/errRplot.R")
save(errF, file = "./ErrorPlots/errF.R")
save(errR, file = "./ErrorPlots/errR.R")

#dereplicate reads
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dir.create("./Dereplicated")
save(derepFs, file = "./Dereplicated/derepFs.R")
save(derepRs, file = "./Dereplicated/derepRs.R")

#apply core sample inference algorithm 
#lots of diagnostics within this "dada" class object
dadaFs <- dada(derepFs, err=errF, multithread=10)
dadaRs <- dada(derepRs, err=errR, multithread=10)
dir.create("./DadaObjects")
save(dadaFs, file = "./DadaObjects/dadaFs.R")
save(dadaRs, file = "./DadaObjects/dadaRs.R")

#and merge
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
save(mergers, file = "./DadaObjects/mergers.R")

#make amplicon sequence variant table
seqtab <- makeSequenceTable(mergers)
dir.create("./SequenceTables")
save(seqtab, file = "./SequenceTables/seqtab.R")

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=10, verbose=TRUE)
save(seqtab.nochim, file = "./SequenceTables/seqtab.nochim.R")

#reload tables back into R 
setwd("~/ArabidopsisWheat/Fungi/")
load("./SequenceTables/seqtab.R")
load("./DadaObjects/dadaFs.R")
load("./DadaObjects/dadaRs.R")
load("./DadaObjects/mergers.R")
load("./Filtered/out.R")
load("./SequenceTables/seqtab.nochim.R")

#Final track of reads throughout the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out[,1:2], sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
save(track, file = "./track.R")


#Fungal taxonomy assignment
rm(list = ls())
library(ShortRead)
library(dada2)
setwd("~/ArabidopsisWheat/Fungi")
load("./SequenceTables/seqtab.nochim.R")
taxa <- assignTaxonomy(seqtab.nochim,
                       "../../References/UNITE_ITS_v8",
                       multithread = 10, verbose = TRUE)

#inspect assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
#if happy save file out for use
save(taxa, file  = "./taxa.R")


#Fungal sequences
#rarefaction and rarefication.
rm(list = ls())
setwd("")
library(dada2)
library(RColorBrewer)
library(phyloseq)
library(ggplot2)
library(vegan)
library(ape)
library(Biostrings)
load("./AnalysisFungiTax002/taxa.R")
load("./AnalysisFungiTax002/seqtab.nochim.R")
sample.info <- read.csv("AnalysisFungiTax002/sample.info.csv")
rownames(sample.info) <- sample.info$Name

#there are reads in the control samples.
#need to subtract a mean of control reads from all the samples
#Must be careful here so as to account for each sample's (and control's) total read number
#calculate the proportion of reads from each OTU (in controls)
#remove that proportion of reads from each OTU in every sample.
controls <- seqtab.nochim[49, ]
prop.controls <- controls / sum(controls)
seqtab.nochim <- seqtab.nochim[1:48, ]
for(i in 1:length(rownames(seqtab.nochim))){
  seqtab.nochim[i, ] <- seqtab.nochim[i, ] - (sum(seqtab.nochim[i, ]) * prop.controls)
  print(i)
}
sample.info <- sample.info[which(rownames(sample.info) %in% rownames(seqtab.nochim)), ]
seqtab.nochim <- round(seqtab.nochim)
seqtab.nochim[which(seqtab.nochim < 0)] <- 0

#Identifying the sequences that only classify to the phylum level and remove them all.
#Removes anything spurious or anything likely to be a plant read.
potential.plant <-  rownames(taxa)[is.na(taxa[, "Class"])]
potential.plant <- DNAStringSet(potential.plant)
taxa <- taxa[-which(rownames(taxa) %in% potential.plant), ]
seqtab.nochim <- seqtab.nochim[, which(colnames(seqtab.nochim) %in% rownames(taxa))]

#remove reads that occur once or in only one sample
seqtab.nochim <- seqtab.nochim[, which(colSums(seqtab.nochim) > 1)]
logical.test <- seqtab.nochim != 0
seqtab.nochim <- seqtab.nochim[, colSums(logical.test) > 1]

save(taxa, file = "./AnalysisFungiTax002/taxa.noplant.R")
save(seqtab.nochim, file = "./AnalysisFungiTax002/seqtab.nochim.noplant.R")
load("./AnalysisFungiTax002/seqtab.nochim.noplant.R")

#rarefaction curves
depths <- c(10, 100, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,
            10000, 12000, 14000, 16000, 18000, 20000,
            25000, 30000, 35000, 40000, 45000, 50000)
number.samples <- length(seqtab.nochim[, 1])
number.depths <- length(depths)
result.matrix <- matrix(0, number.samples, number.depths)
colnames(result.matrix) <- depths
rownames(result.matrix) <- rownames(seqtab.nochim)
for(i in 1:number.depths){
  for(q in 1:number.samples){
    this.sample <- seqtab.nochim[q, ]
    abundance <- 0
    if(sum(this.sample) >= depths[i]){
      rarefied <- rrarefy(this.sample, depths[i])
      abundance <- sum(rarefied > 0)
    }
    result.matrix[q, i] <- abundance
  }
}
result.matrix[which(result.matrix == 0)] <- NA
plot(depths, result.matrix[1, ], type = "l", lwd = 0.25, ylim = c(0, max(result.matrix, na.rm = TRUE)*1.1),
     ylab = "Observed species abundance",
     xlab = "Rarefication depth", 
     main = "Fungal Rarefaction Curve")
for(i in 2:number.samples){
  lines(depths, result.matrix[i, ], lwd = 0.25)
}
abline(v = 4726, col = "red")
#rarefy at 4726 keeps all but one sample which has 666 reads which is clearly a very bad sample

rowSums(seqtab.nochim)

set.seed(10041991)
totals <- rowSums(seqtab.nochim)
seqtab.nochim <- seqtab.nochim[which(totals >= 4203), ]
sample.info <- sample.info[which(rownames(sample.info) %in% rownames(seqtab.nochim)), ]
orig.seqtab.nochim <- seqtab.nochim
seqtab.nochim <- rrarefy(orig.seqtab.nochim, sample = 4203)
#remove ASV's with no counts associated now
totals2 <- colSums(seqtab.nochim)
seqtab.nochim <- seqtab.nochim[, which(totals2 != 0)]
#saveout sequence table and sample info for use later now rarefied and other alterations made
save(seqtab.nochim, file = "./AnalysisFungiTax002/seqtab.nochim.rarefied.R")
save(sample.info, file = "./AnalysisFungiTax002/sample.info.rarefied.R")
save(taxa, file = "./AnalysisFungiTax002/taxa.noplant.R")


##############################################################################################################################
#Alpha diversity

rm(list=ls())
setwd("")
library(dada2)
library(RColorBrewer)
library(phyloseq)
library(ggplot2)
library(vegan)
load("./AnalysisFungiTax002/taxa.noplant.R")
load("./AnalysisFungiTax002/seqtab.nochim.rarefied.R")
load("./AnalysisFungiTax002/sample.info.rarefied.R")

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(sample.info))

#Remove rep 1 looks very bad on the PCoA clearly skewing things. (Consider for fungi too)
ps <- subset_samples(ps, sample_data(ps)$BioRep != "1")

#calculate various alpha diversity measures included
alpha_diversities <- estimate_richness(ps, split = TRUE)

#add sample information to the data frame
sample.info$Plant <- as.character(sample.info$Plant)
sample.info$Soil <- as.character(sample.info$Soil)
sample.info$BioRep <- as.character(sample.info$BioRep)
sample.info$Fraction <- as.character(sample.info$Fraction)
rownames(alpha_diversities) <- gsub("X", "", rownames(alpha_diversities))
rownames(alpha_diversities) <- gsub("\\.", "-", rownames(alpha_diversities))
alpha_diversities$Plant <- ""
alpha_diversities$Soil <- ""
alpha_diversities$BioRep <- ""
alpha_diversities$Fraction <- ""
for(i in 1:length(rownames(alpha_diversities))){
  alpha_diversities[i, 10:13] <- as.character(sample.info[which(rownames(sample.info) == rownames(alpha_diversities)[i]), 1:4])
}
#remove rhizosphere unplanted
alpha_diversities <- alpha_diversities[-c(which(alpha_diversities$Fraction == "Rhizosphere " & alpha_diversities$Plant == "Unplanted")), ]
alpha_diversities$Plant <- factor(alpha_diversities$Plant, levels = c("Unplanted", "A. thaliana", "T. aestivum"), labels = c("Unplanted", "A. thaliana", "T. aestivum"))
alpha_diversities$Fraction <- factor(alpha_diversities$Fraction, levels = c("Soil", "Rhizosphere ", "Endosphere"), labels = c("Soil", "Rhizosphere", "Endosphere"))

#pretty graphs for abundance values
ggplot(data = alpha_diversities, aes(x = Fraction, y = Observed, fill = Soil)) + 
  geom_boxplot(outlier.colour = NA) + 
  facet_wrap("Plant", scales = "fixed") +
  geom_point(aes(y = Observed, x = Fraction, fill = Soil, col = BioRep), 
             pch = 19, position = position_dodge(width = 0.75)) + 
  scale_color_grey(name = "Replicate", labels = c("1", "2", "3", "4")) +
  scale_fill_manual(name = "Soil", values = c("darkgoldenrod", "navy")) +
  xlab("Plant") +
  ylab("Observed species richness") +
  theme(axis.text.x = element_blank(), axis.ticks.x.bottom = element_blank()) + 
  ggtitle("Alpha diversities")

fitAll <- aov(Observed ~ Fraction * Soil * Plant, data = alpha_diversities)
summary(fitAll)
TukeyHSD(fitAll)

fitAll <- aov(Observed ~ Soil * Plant, data = alpha_diversities[which(alpha_diversities$Fraction == "Endosphere"), ])
summary(fitAll)
TukeyHSD(fitAll)

fitAll <- aov(Observed ~ Soil * Plant, data = alpha_diversities[which(alpha_diversities$Fraction == "Rhizosphere"), ])
summary(fitAll)
TukeyHSD(fitAll)

fitAll <- aov(Observed ~ Soil * Plant, data = alpha_diversities[which(alpha_diversities$Fraction == "Soil"), ])
summary(fitAll)
TukeyHSD(fitAll)

w1 <- pairwise.wilcox.test(alpha_diversities$Shannon, alpha_diversities$Fraction, p.adjust.method = "BH")
w2 <- pairwise.wilcox.test(alpha_diversities$Shannon, alpha_diversities$Plant, p.adjust.method = "BH")
w3 <- pairwise.wilcox.test(alpha_diversities$Shannon, alpha_diversities$Soil, p.adjust.method = "BH")
w4 <- pairwise.wilcox.test(alpha_diversities$Shannon, paste0(alpha_diversities$Fraction, alpha_diversities$Plant), 
                     p.adjust.method = "BH")


capture.output(print("Alpha statistics (No interaction effect found)"), file = "./AnalysisFungiTax002/alpha.statistics.txt")
capture.output(summary(fitAll), file = "./AnalysisFungiTax002/alpha.statistics.txt", append = TRUE)
capture.output(TukeyHSD(fitAll), file = "./AnalysisFungiTax002/alpha.statistics.txt", append = TRUE)
capture.output(w1, file = "./AnalysisFungiTax002/alpha.statistics.txt", append = TRUE)
capture.output(w2, file = "./AnalysisFungiTax002/alpha.statistics.txt", append = TRUE)
capture.output(w3, file = "./AnalysisFungiTax002/alpha.statistics.txt", append = TRUE)
capture.output(w4, file = "./AnalysisFungiTax002/alpha.statistics.txt", append = TRUE)

##############################################################################################
#Beta diversity

rm(list=ls())
setwd("")
library(dada2)
library(RColorBrewer)
library(phyloseq)
library(ggplot2)
library(vegan)
load("./AnalysisFungiTax002/taxa.noplant.R")
load("./AnalysisFungiTax002/seqtab.nochim.rarefied.R")
load("./AnalysisFungiTax002/sample.info.rarefied.R")

capture.output(print("Beta statistics"), file = "./AnalysisFungiTax002/beta.statistics.txt")

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(sample.info))


#Remove rep 1 looks very bad on the PCoA clearly skewing things. (Consider for fungi too)
ps <- subset_samples(ps, sample_data(ps)$BioRep != "1")

#Remove unplanted rhizosphere samples becsue it is skewing the statistics becasue it clusters so well with soil samples
ps <- subset_samples(ps, (sample_data(ps)$Fraction == "Rhizosphere " & sample_data(ps)$Plant == "Unplanted") == FALSE)

#subset by endosphere
psendo <- subset_samples(ps, sample_data(ps)$Fraction %in% c("Endosphere"))
#beta diversities and permanovas
test <- adonis(phyloseq::distance(psendo, method = "bray") ~ Plant + Soil, 
               data = data.frame(sample_data(psendo)), method  = "bray", permutations = 1000)
test
capture.output(test, file = "./AnalysisFungiTax002/beta.statistics.txt", append = TRUE)
#transform sample counts
ps.prop <- transform_sample_counts(psendo, function(otu) otu/sum(otu))
#ordinate
ord.nmds.bray <- ordinate(ps.prop, method="MDS", distance="bray")
#percentage explained variation for each axis is in
ord.nmds.bray$values[1:2, "Relative_eig"]
#make an object for ggplot
vectors <- as.data.frame(ord.nmds.bray$vectors)
vectors <- cbind(vectors, sample_data(ps.prop))
#plot
ggplot(data = vectors) +
  geom_point(aes(x = Axis.1, y = Axis.2,
                 shape = factor(Plant),
                 col = factor(Soil)), size = 4)+
  scale_color_manual(name = "Soil Type", values = c("darkgoldenrod", "navy")) +
  scale_shape_discrete(name = "Plant") +
  theme_bw()


#subset by rhizosphere
psrhiz <- subset_samples(ps, sample_data(ps)$Fraction %in% c("Rhizosphere ") & sample_data(ps)$Plant != "Unplanted")
#beta diversities and permanovas
test <- adonis(phyloseq::distance(psrhiz, method = "bray") ~ Soil + Plant, 
               data = data.frame(sample_data(psrhiz)), method  = "bray", permutations = 1000)
test
capture.output(test, file = "./AnalysisFungiTax002/beta.statistics.txt", append = TRUE)
#transform sample counts
ps.prop <- transform_sample_counts(psrhiz, function(otu) otu/sum(otu))
#ordinate
ord.nmds.bray <- ordinate(ps.prop, method="MDS", distance="bray")
#percentage explained variation for each axis is in
ord.nmds.bray$values[1:2, "Relative_eig"]
#make an object for ggplot
vectors <- as.data.frame(ord.nmds.bray$vectors)
vectors <- cbind(vectors, sample_data(ps.prop))
#plot
ggplot(data = vectors) +
  geom_point(aes(x = Axis.1, y = Axis.2,
                 shape = factor(Plant),
                 col = factor(Soil)), size = 4)+
  scale_color_manual(name = "Soil Type", values = c("darkgoldenrod", "navy", "red")) +
  scale_shape_discrete(name = "Plant") +
  theme_bw()



#subset by bulk soil
pssoil <- subset_samples(ps, sample_data(ps)$Fraction %in% c("Soil"))
#beta diversities and permanovas
test <- adonis(phyloseq::distance(pssoil, method = "bray") ~ Plant + Soil, 
               data = data.frame(sample_data(pssoil)), method  = "bray", permutations = 1000)
test
capture.output(test, file = "./AnalysisFungiTax002/beta.statistics.txt", append = TRUE)
#transform sample counts
ps.prop <- transform_sample_counts(pssoil, function(otu) otu/sum(otu))
#ordinate
ord.nmds.bray <- ordinate(ps.prop, method="MDS", distance="bray")
#percentage explained variation for each axis is in
ord.nmds.bray$values[1:2, "Relative_eig"]
#make an object for ggplot
vectors <- as.data.frame(ord.nmds.bray$vectors)
vectors <- cbind(vectors, sample_data(ps.prop))
#plot
ggplot(data = vectors) +
  geom_point(aes(x = Axis.1, y = Axis.2,
                 shape = factor(Plant),
                 col = factor(Soil)), size = 4)+
  scale_color_manual(name = "Soil Type", values = c("darkgoldenrod", "navy", "red")) +
  scale_shape_discrete(name = "Plant") +
  theme_bw()



ord.cap.bray <- ordinate(ps.prop, method="CAP", distance="bray", formula = ~ Plant + Soil)
ord.cap.bray
vectors <- ord.cap.bray$CCA$wa
vectors <- cbind(vectors, sample_data(ps.prop))

ggplot(data = vectors) +
  geom_point(aes(x = CAP3, y = CAP2,
                 shape = factor(Plant),
                 col = factor(Soil)), size = 4)+
  scale_color_manual(name = "Soil Type", values = c("darkgoldenrod", "navy")) +
  scale_shape_discrete(name = "Plant") +
  xlab("CAP3 [53.07% constrained variation]") +
  ylab("CAP2 [25.41% constrained variation]") +
  ggtitle("24.92% total variation") +
  theme_bw()


capture.output(ord.cap.bray, file = "./AnalysisFungiTax002/beta.statistics.txt", append = TRUE)

#quick permanova to test the plant effect without the unplanted
pssoil <- subset_samples(ps, sample_data(ps)$Fraction %in% c("Soil") & sample_data(ps)$Plant != "Unplanted")
#beta diversities and permanovas
test <- adonis(phyloseq::distance(pssoil, method = "bray") ~ Plant + Soil, 
               data = data.frame(sample_data(pssoil)), method  = "bray", permutations = 100)
test
capture.output(test, file = "./AnalysisFungiTax002/beta.statistics.txt", append = TRUE)

#############################################################################################
#Stacked bar charts
rm(list=ls())
setwd("")
library(dada2)
library(RColorBrewer)
library(phyloseq)
library(ggplot2)
library(vegan)
load("./AnalysisFungiTax002/taxa.noplant.R")
load("./AnalysisFungiTax002/seqtab.nochim.rarefied.R")
load("./AnalysisFungiTax002/sample.info.rarefied.R")

#remove biorep1 first
sample.info <- sample.info[which(sample.info$BioRep != 1), ]
seqtab.nochim <- seqtab.nochim[which(rownames(seqtab.nochim) %in% rownames(sample.info)), ]

#filter for highest abundance reads
total.abundance <- colSums(seqtab.nochim)
rel.abundance <- (total.abundance / sum(total.abundance)) * 100
seqtab.small <- seqtab.nochim[, which(rel.abundance > 0.1)]
#seqtab.small <- seqtab.nochim

ps <- phyloseq(otu_table(seqtab.small, taxa_are_rows=FALSE), 
               sample_data(sample.info), tax_table(taxa))

sample_data(ps)$SampleType <- paste(sample_data(ps)$Plant, sample_data(ps)$Soil, sample_data(ps)$Fraction)
sample_data(ps)$SampleType2 <- paste(sample_data(ps)$Plant, sample_data(ps)$Soil, sample_data(ps)$Fraction, sample_data(ps)$BioRep)
psall <- merge_samples(ps, "SampleType")
psall <- transform_sample_counts(psall, function(x) x / sum(x))
sample_data(psall)$Plant <- factor(sample_data(psall)$Plant, levels = c("3", "1", "2"), labels = c("Unplanted", "A. thaliana", "T. aestivum"))
sample_data(psall)$Soil <- factor(sample_data(psall)$Soil, levels = c("1", "2"), labels = c("WG", "WG+"))
sample_data(psall)$Fraction <- factor(sample_data(psall)$Fraction, levels = c("3", "2", "1"), labels = c("Soil", "Rhizosphere", "Endosphere"))
sample_data(psall)$SampleType <- paste(sample_data(psall)$Plant, sample_data(psall)$Fraction)
psall <- subset_samples(psall, sample_data(psall)$SampleType != "Unplanted Rhizosphere")
sample_data(psall)$SampleType <- paste(sample_data(psall)$Soil, sample_data(psall)$Fraction)
sample_data(psall)$SampleType <- factor(sample_data(psall)$SampleType, 
                                        levels = c("WG Soil", "WG+ Soil", "WG Rhizosphere", "WG+ Rhizosphere", "WG Endosphere", "WG+ Endosphere"), 
                                        labels = c("WG Soil", "WG+ Soil", "WG Rhizosphere", "WG+ Rhizosphere", "WG Endosphere", "WG+ Endosphere"))
ps <- transform_sample_counts(ps, function(x) x / sum(x))

#barplot
plot_bar(psall, x="SampleType", fill="Class") +
  facet_wrap(~Plant, scales="free_x") + 
  geom_bar(stat = "identity") +
  xlab("Treatment") + 
  ylab("Relative abundance")


plot_bar(ps, x="SampleType2", fill="Class") +
  facet_wrap(~Plant, scales="free_x") + 
  geom_bar(stat = "identity") +
  xlab("Treatment") + 
  ylab("Relative abundance")



####################################################################################################
#Differential abundance analysis 
######## New method do all soil comparisons at once
rm(list=ls())
setwd("")
library(dada2)
library(RColorBrewer)
library(phyloseq)
library(ggplot2)
library(vegan)
library(DESeq2)
load("./AnalysisFungiTax002/taxa.noplant.R")
load("./AnalysisFungiTax002/seqtab.nochim.rarefied.R")
load("./AnalysisFungiTax002/sample.info.rarefied.R")

ASV1 <- colnames(seqtab.nochim)
Names <- paste0("ASV", 1:length(ASV1))
ASV1 <- data.frame(ASV1, Names)
save(ASV1, file = "./AnalysisFungiTax002/ASV1.R")

#Need to subset for each compartment and remove the ASV's that are low-abundance to avoid the sparse matrix problem
#Also remove rep 1.
sample.info <- sample.info[which(sample.info$BioRep != "1"), ]
sample.info$combined  <- paste0(sample.info$Plant, sample.info$Fraction, sample.info$Soil)
sample.info.endo <- sample.info[which(sample.info$Fraction == "Endosphere"), ]
sample.info.rhiz <- sample.info[which(sample.info$Fraction == "Rhizosphere "), ]
sample.info.soil <- sample.info[which(sample.info$Fraction == "Soil"), ]
seqtab.endo <- seqtab.nochim[which(rownames(seqtab.nochim) %in% rownames(sample.info.endo)), ]
seqtab.rhiz <- seqtab.nochim[which(rownames(seqtab.nochim) %in% rownames(sample.info.rhiz)), ]
seqtab.soil <- seqtab.nochim[which(rownames(seqtab.nochim) %in% rownames(sample.info.soil)), ]
seqtab.endo <- seqtab.endo[, which(colSums(seqtab.endo) > 0)]
seqtab.rhiz <- seqtab.rhiz[, which(colSums(seqtab.rhiz) > 0)]
seqtab.soil <- seqtab.soil[, which(colSums(seqtab.soil) > 0)]



#make phyloseq objects
psendo <- phyloseq(otu_table(seqtab.endo, taxa_are_rows=FALSE), 
                   sample_data(sample.info.endo), 
                   tax_table(taxa))
psrhiz <- phyloseq(otu_table(seqtab.rhiz, taxa_are_rows=FALSE), 
                   sample_data(sample.info.rhiz), 
                   tax_table(taxa))
pssoil <- phyloseq(otu_table(seqtab.soil, taxa_are_rows=FALSE), 
                   sample_data(sample.info.soil), 
                   tax_table(taxa))

#Make DESeq object
ds1 <- phyloseq_to_deseq2(psendo, ~ combined)
dds1 <- estimateSizeFactors(ds1, type = "poscount")
dds1 <- estimateDispersions(dds1, fitType = "local")
dds1 <- nbinomWaldTest(dds1)

comparisons <- list(c("combined", "A. thalianaEndosphereWG", "A. thalianaEndosphereWG+"),
                    c("combined", "T. aestivumEndosphereWG", "T. aestivumEndosphereWG+"),
                    c("combined", "A. thalianaEndosphereWG",  "T. aestivumEndosphereWG"),
                    c("combined", "A. thalianaEndosphereWG+", "T. aestivumEndosphereWG+"))


############ Deseq fucntion to pull out results for all different comparisons
make.result.df <- function(comparison){
  res1 <- DESeq2::results(dds1, contrast = comparison)
  df1 <- as.data.frame(res1)
  df1$taxon <- rownames(df1)
  rownames(df1) <- paste("ArbOTU", 1:length(rownames(df1)), sep = "")
  df <- df1[which(df1$padj < 0.05), ]
  if(nrow(df) > 0){
    df$Phylum <- "phylum"
    df$Class <- "class"
    df$Order <- "order"
    df$Family <- "family"
    df$Genus <- "genus"
    for(x in 1:length(rownames(df))){
      this.data <- taxa[which(rownames(taxa) == df$taxon[x]), ]
      df$Phylum[x] <- this.data[[2]]
      df$Class[x] <- this.data[[3]]
      df$Order[x] <- this.data[[4]]
      df$Family[x] <- this.data[[5]]
      df$Genus[x] <- this.data[[6]]
      print(x)
    }
    df$Comparison <- paste(comparison, collapse = " ")
  }
  return(df)
}
#loops for every comparison in comparisons list
final.df <- c()
for(i in 1:length(comparisons)){
  this.result <- make.result.df(comparisons[[i]])
  if(nrow(this.result) > 0){
    final.df <- rbind(final.df, this.result)
  }
  print(i)
}
final.df

load("./AnalysisFungiTax002/ASV1.R")
#rename with standardise ASV names
final.df$taxon <- as.character(final.df$taxon)
final.df$ASV <- ""
for(i in 1:length(final.df$taxon)){
  final.df$ASV[i] <- as.character(ASV1$Names[which(as.character(ASV1$ASV1) == as.character(final.df$taxon[i]))])
}


write.csv(final.df, file = "./AnalysisFungiTax002/final.df.fungi.endo.csv")
final.df <- read.csv("./AnalysisFungiTax002/final.df.fungi.endo.csv")


#rhizosphere

#Make DESeq object
ds1 <- phyloseq_to_deseq2(psrhiz, ~ combined)
dds1 <- estimateSizeFactors(ds1, type = "poscount")
dds1 <- estimateDispersions(dds1, fitType = "local")
dds1 <- nbinomWaldTest(dds1)


comparisons <- list(c("combined", "UnplantedRhizosphere WG", "UnplantedRhizosphere WG+"),
                    c("combined", "A. thalianaRhizosphere WG", "A. thalianaRhizosphere WG+"),
                    c("combined", "T. aestivumRhizosphere WG", "T. aestivumRhizosphere WG+"),
                    c("combined", "UnplantedRhizosphere WG", "A. thalianaRhizosphere WG"),
                    c("combined", "UnplantedRhizosphere WG", "T. aestivumRhizosphere WG"),
                    c("combined", "A. thalianaRhizosphere WG", "T. aestivumRhizosphere WG"),
                    c("combined", "UnplantedRhizosphere WG+", "A. thalianaRhizosphere WG+"),
                    c("combined", "UnplantedRhizosphere WG+", "T. aestivumRhizosphere WG+"),
                    c("combined", "A. thalianaRhizosphere WG+", "T. aestivumRhizosphere WG+"))


############ Deseq fucntion to pull out results for all different comparisons
make.result.df <- function(comparison){
  res1 <- DESeq2::results(dds1, contrast = comparison)
  df1 <- as.data.frame(res1)
  df1$taxon <- rownames(df1)
  rownames(df1) <- paste("ArbOTU", 1:length(rownames(df1)), sep = "")
  df <- df1[which(df1$padj < 0.05), ]
  if(nrow(df) > 0){
    df$Phylum <- "phylum"
    df$Class <- "class"
    df$Order <- "order"
    df$Family <- "family"
    df$Genus <- "genus"
    for(x in 1:length(rownames(df))){
      this.data <- taxa[which(rownames(taxa) == df$taxon[x]), ]
      df$Phylum[x] <- this.data[[2]]
      df$Class[x] <- this.data[[3]]
      df$Order[x] <- this.data[[4]]
      df$Family[x] <- this.data[[5]]
      df$Genus[x] <- this.data[[6]]
      print(x)
    }
    df$Comparison <- paste(comparison, collapse = " ")
  }
  return(df)
}
#loops for every comparison in comparisons list
final.df <- c()
for(i in 1:length(comparisons)){
  this.result <- make.result.df(comparisons[[i]])
  if(nrow(this.result) > 0){
    final.df <- rbind(final.df, this.result)
  }
  print(i)
}
final.df

load("./AnalysisFungiTax002/ASV1.R")
#rename with standardise ASV names
final.df$taxon <- as.character(final.df$taxon)
final.df$ASV <- ""
for(i in 1:length(final.df$taxon)){
  final.df$ASV[i] <- as.character(ASV1$Names[which(as.character(ASV1$ASV1) == as.character(final.df$taxon[i]))])
}

write.csv(final.df, file = "./AnalysisFungiTax002/final.df.fungi.rhiz.csv")
final.df <- read.csv("./AnalysisFungiTax002/final.df.fungi.rhiz.csv")


#bulk soil

#Make DESeq object
ds1 <- phyloseq_to_deseq2(pssoil, ~ combined)
dds1 <- estimateSizeFactors(ds1, type = "poscount")
dds1 <- estimateDispersions(dds1, fitType = "local")
dds1 <- nbinomWaldTest(dds1)


comparisons <- list(c("combined", "UnplantedSoilWG", "UnplantedSoilWG+"),
                    c("combined", "A. thalianaSoilWG", "A. thalianaSoilWG+"),
                    c("combined", "T. aestivumSoilWG", "T. aestivumSoilWG+"),
                    c("combined", "UnplantedSoilWG", "A. thalianaSoilWG"),
                    c("combined", "UnplantedSoilWG", "T. aestivumSoilWG"),
                    c("combined", "A. thalianaSoilWG", "T. aestivumSoilWG"),
                    c("combined", "UnplantedSoilWG+", "A. thalianaSoilWG+"),
                    c("combined", "UnplantedSoilWG+", "T. aestivumSoilWG+"),
                    c("combined", "A. thalianaSoilWG+", "T. aestivumSoilWG+"))



############ Deseq fucntion to pull out results for all different comparisons
make.result.df <- function(comparison){
  res1 <- DESeq2::results(dds1, contrast = comparison)
  df1 <- as.data.frame(res1)
  df1$taxon <- rownames(df1)
  rownames(df1) <- paste("ArbOTU", 1:length(rownames(df1)), sep = "")
  df <- df1[which(df1$padj < 0.05), ]
  if(nrow(df) > 0){
    df$Phylum <- "phylum"
    df$Class <- "class"
    df$Order <- "order"
    df$Family <- "family"
    df$Genus <- "genus"
    for(x in 1:length(rownames(df))){
      this.data <- taxa[which(rownames(taxa) == df$taxon[x]), ]
      df$Phylum[x] <- this.data[[2]]
      df$Class[x] <- this.data[[3]]
      df$Order[x] <- this.data[[4]]
      df$Family[x] <- this.data[[5]]
      df$Genus[x] <- this.data[[6]]
      print(x)
    }
    df$Comparison <- paste(comparison, collapse = " ")
  }
  return(df)
}
#loops for every comparison in comparisons list
final.df <- c()
for(i in 1:length(comparisons)){
  this.result <- make.result.df(comparisons[[i]])
  if(nrow(this.result) > 0){
    final.df <- rbind(final.df, this.result)
  }
  print(i)
}
final.df


load("./AnalysisFungiTax002/ASV1.R")
#rename with standardise ASV names
final.df$taxon <- as.character(final.df$taxon)
final.df$ASV <- ""
for(i in 1:length(final.df$taxon)){
  final.df$ASV[i] <- as.character(ASV1$Names[which(as.character(ASV1$ASV1) == as.character(final.df$taxon[i]))])
}

write.csv(final.df, file = "./AnalysisFungiTax002/final.df.fungi.soil.csv")
final.df <- read.csv("./AnalysisFungiTax002/final.df.fungi.soil.csv")


###################Heatmaps for everything together, endosphere differentially abundant only. (With pattern in all samples)
rm(list=ls())
setwd("")
library(dada2)
library(RColorBrewer)
library(phyloseq)
library(ggplot2)
library(vegan)
library(DESeq2)
library(factoextra)
library(gplots)
load("./AnalysisFungiTax002/taxa.noplant.R")
load("./AnalysisFungiTax002/seqtab.nochim.rarefied.R")
load("./AnalysisFungiTax002/sample.info.rarefied.R")
load("./AnalysisFungiTax002/ASV1.R")
final.df <- read.csv("./AnalysisFungiTax002/final.df.fungi.endo.csv")
final.df2 <- read.csv("./AnalysisFungiTax002/final.df.fungi.rhiz.csv")
final.df3 <- read.csv("./AnalysisFungiTax002/final.df.fungi.soil.csv")
final.df <- rbind(final.df, final.df2, final.df3)
#Get rid of the unplanted rhizosphere comparisons
final.df <- final.df[-c(grep("UnplantedRhizosphere", final.df$Comparison)), ]
to.plot <- unique(final.df$taxon)

seqtab.nochim <- seqtab.nochim[, which(colnames(seqtab.nochim) %in% to.plot)]
#remove rep 1 and remove unplanted rhizsphere
sample.info <- sample.info[which(sample.info$BioRep != "1"), ]

eLNwheatSamples <- as.character(rownames(sample.info[which(sample.info$Soil == "WG" & sample.info$Plant == "T. aestivum" & sample.info$Fraction == "Endosphere"),]))
eLNarabiSamples <- as.character(rownames(sample.info[which(sample.info$Soil == "WG" & sample.info$Plant == "A. thaliana" & sample.info$Fraction == "Endosphere"),]))
eHNwheatSamples <- as.character(rownames(sample.info[which(sample.info$Soil == "WG+" & sample.info$Plant == "T. aestivum" & sample.info$Fraction == "Endosphere"),]))
eHNarabiSamples <- as.character(rownames(sample.info[which(sample.info$Soil == "WG+" & sample.info$Plant == "A. thaliana" & sample.info$Fraction == "Endosphere"),]))
rLNwheatSamples <- as.character(rownames(sample.info[which(sample.info$Soil == "WG" & sample.info$Plant == "T. aestivum" & sample.info$Fraction == "Rhizosphere "),]))
rLNarabiSamples <- as.character(rownames(sample.info[which(sample.info$Soil == "WG" & sample.info$Plant == "A. thaliana" & sample.info$Fraction == "Rhizosphere "),]))
rHNwheatSamples <- as.character(rownames(sample.info[which(sample.info$Soil == "WG+" & sample.info$Plant == "T. aestivum" & sample.info$Fraction == "Rhizosphere "),]))
rHNarabiSamples <- as.character(rownames(sample.info[which(sample.info$Soil == "WG+" & sample.info$Plant == "A. thaliana" & sample.info$Fraction == "Rhizosphere "),]))
sLNwheatSamples <- as.character(rownames(sample.info[which(sample.info$Soil == "WG" & sample.info$Plant == "T. aestivum" & sample.info$Fraction == "Soil"),]))
sLNarabiSamples <- as.character(rownames(sample.info[which(sample.info$Soil == "WG" & sample.info$Plant == "A. thaliana" & sample.info$Fraction == "Soil"),]))
sHNwheatSamples <- as.character(rownames(sample.info[which(sample.info$Soil == "WG+" & sample.info$Plant == "T. aestivum" & sample.info$Fraction == "Soil"),]))
sHNarabiSamples <- as.character(rownames(sample.info[which(sample.info$Soil == "WG+" & sample.info$Plant == "A. thaliana" & sample.info$Fraction == "Soil"),]))
sHNunplantedSamples <- as.character(rownames(sample.info[which(sample.info$Soil == "WG+" & sample.info$Plant == "Unplanted" & sample.info$Fraction == "Soil"),]))
sLNunplantedSamples <- as.character(rownames(sample.info[which(sample.info$Soil == "WG" & sample.info$Plant == "Unplanted" & sample.info$Fraction == "Soil"),]))
eLNwheatSums <- colSums(seqtab.nochim[which(rownames(seqtab.nochim) %in% eLNwheatSamples), ]) / length(eLNwheatSamples)
eLNarabiSums <- colSums(seqtab.nochim[which(rownames(seqtab.nochim) %in% eLNarabiSamples), ]) / length(eLNarabiSamples)
eHNwheatSums <- colSums(seqtab.nochim[which(rownames(seqtab.nochim) %in% eHNwheatSamples), ]) / length(eHNwheatSamples)
eHNarabiSums <- colSums(seqtab.nochim[which(rownames(seqtab.nochim) %in% eHNarabiSamples), ]) / length(eHNarabiSamples)
rLNwheatSums <- colSums(seqtab.nochim[which(rownames(seqtab.nochim) %in% rLNwheatSamples), ]) / length(rLNwheatSamples)
rLNarabiSums <- colSums(seqtab.nochim[which(rownames(seqtab.nochim) %in% rLNarabiSamples), ]) / length(rLNarabiSamples)
rHNwheatSums <- colSums(seqtab.nochim[which(rownames(seqtab.nochim) %in% rHNwheatSamples), ]) / length(rHNwheatSamples)
rHNarabiSums <- colSums(seqtab.nochim[which(rownames(seqtab.nochim) %in% rHNarabiSamples), ]) / length(rHNarabiSamples)
sLNwheatSums <- colSums(seqtab.nochim[which(rownames(seqtab.nochim) %in% sLNwheatSamples), ]) / length(sLNwheatSamples)
sLNarabiSums <- colSums(seqtab.nochim[which(rownames(seqtab.nochim) %in% sLNarabiSamples), ]) / length(sLNarabiSamples)
sHNwheatSums <- colSums(seqtab.nochim[which(rownames(seqtab.nochim) %in% sHNwheatSamples), ]) / length(sHNwheatSamples)
sHNarabiSums <- colSums(seqtab.nochim[which(rownames(seqtab.nochim) %in% sHNarabiSamples), ]) / length(sHNarabiSamples)
sHNunplantedSums <- colSums(seqtab.nochim[which(rownames(seqtab.nochim) %in% sHNunplantedSamples), ]) / length(sHNunplantedSamples)
sLNunplantedSums <- colSums(seqtab.nochim[which(rownames(seqtab.nochim) %in% sLNunplantedSamples), ]) / length(sLNunplantedSamples)

#need to add in a line to remove all teh rown that sum to 0
my.matrix.s <- cbind(sLNunplantedSums, sHNunplantedSums, sLNwheatSums,  sHNwheatSums, sLNarabiSums, sHNarabiSums)
my.matrix.r <- cbind(rLNwheatSums, rHNwheatSums, rLNarabiSums,  rHNarabiSums)
my.matrix.e <- cbind(eLNwheatSums, eHNwheatSums, eLNarabiSums, eHNarabiSums)
mat.e <- t(scale(x = t(my.matrix.e), center = TRUE, scale = TRUE))
mat.r <- t(scale(x = t(my.matrix.r), center = TRUE, scale = TRUE))
mat.s <- t(scale(x = t(my.matrix.s), center = TRUE, scale = TRUE))
mat <- cbind(mat.s, mat.r, mat.e)

#need to get rid of Na's for the clustering to work
mat.e[which(is.na(mat.e))] <- 0
mat.r[which(is.na(mat.r))] <- 0
mat.s[which(is.na(mat.s))] <- 0
#for plotting
mat[which(is.na(mat))] <- 0
#for clustering
#weighting each of the compartments for the clustering
mat2 <- mat
mat2[, 1:6] <- mat2[, 1:6] / 3
mat2[, 7:10] <- mat2[, 7:10] / 2



#elbow method
fviz_nbclust(x = mat2, FUNcluster = hcut, method = "wss", k.max = 25) +
  labs(subtitle = "Hierarchical clustering - Elbow method")
hr <- hclust(d = dist(x = mat2, method = "euclidean"))
my_k = 3
mycl <- cutree(hr, k = my_k)
mycolhc_4 <- c("white", "black", "grey56", "red", "darkgoldenrod", "navy", "brown", "green")
mycolhc <- mycolhc_4[as.vector(mycl)]
#Need to read in other databases to get labelling right
load("./AnalysisFungiTax002/ASV1.R")
labels <- ASV1[which(ASV1[, "ASV1"] %in% hr$labels), "Names"]
hr$labels2 <- labels
## Plot heatmap 
mycol <- colorRampPalette(c("yellow", "blue"))(255)
heatmap.2(as.matrix(mat), col = mycol, dendrogram = 'row',
          Rowv = as.dendrogram(hr), Colv = FALSE,
          scale ="none", RowSideColors = mycolhc, 
          density.info="none", trace = "none",
          labRow = hr$labels2, cexCol = 0.9, srtCol = 45)



