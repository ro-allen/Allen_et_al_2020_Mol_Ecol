#Packages ----
library(ggplot2)
library(phyloseq)
library(reshape2)
library(phangorn)
library(ape)
library(gridExtra)
library(vegan)
library(cowplot)
library(knitr)
library(dada2)
library(vegan)
library(picante)

#ggrare source code ----
#' Make a rarefaction curve using ggplot2
#' @param physeq_object A phyloseq class object, from which abundance data are extracted
#' @param step Step Size for sample size in rarefaction curves
#' @param label Default `NULL`. Character string. The name of the variable to map to text labels on the plot. Similar to color option but for plotting text.
#' @param color Default `NULL`. Character string. The name of the variable to map to the colors in the plot. This can be a sample variables among the set returned by sample_variables(physeq_object) or taxonomic rank, among the set returned by rank_names(physeq_object)
#' @param plot default `TRUE`. Logical. Should the graph be plotted
#' @param parallel default `FALSE`. Logical. Should rarefaction be parallelized
#' @param se default `TRUE`. Logical. Should standard errors be calculated.
#' @export

ggrare <- function(physeq_object, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  
  x <- methods::as(phyloseq::otu_table(physeq_object), "matrix")
  if (phyloseq::taxa_are_rows(physeq_object)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  # Get sample data
  if (!is.null(phyloseq::sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(phyloseq::sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  # Add, any custom-supplied plot-mapped variables
  if ( length(color) > 1 ) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  
  if ( length(label) > 1 ) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot2::ggplot(data = data,
                       ggplot2::aes_string(x = "Size",
                                           y = ".S",
                                           group = "Sample",
                                           color = color))
  
  p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")
  
  if (!is.null(label)) {
    p <- p + ggplot2::geom_text(data = labels,
                                ggplot2::aes_string(x = "x",
                                                    y = "y",
                                                    label = label,
                                                    color = color),
                                size = 4, hjust = 0)
  }
  
  p <- p + ggplot2::geom_line()
  if (se) { ## add standard error if available
    p <- p +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                               ymax = ".S + .se",
                                               color = NULL,
                                               fill = color),
                           alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

#Data import----
path = "~/Desktop/PhD/Data/Sequencing Data/geotraces_transect_run/OG3896/OG3896_fastq"
list.files(path)

#Bioconductor pipeline----

#Extract sample names (used in Filtering and trimming)
f.names = as.vector(list.files(path, pattern = "_R1_001.fastq", 
                               full.names = F))
r.names = as.vector(list.files(path, pattern = "_R2_001.fastq", 
                               full.names = F))

##Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

##Extract sample names. Filenames have format: SAMPLENAME_XXX.fastq
setwd("~/Desktop/PhD/Data/Sequencing Data/geotraces_transect_run/OG3896/OG3896_fastq")
list.files()
sn <- read.csv('geotraces_transect_sample_names.csv', header = FALSE)
sample.names = sn[,1]
sample.names = as.vector(sample.names)

#Plot quality scores
qpf = plotQualityProfile(fnFs[1:20]) 
qpr = plotQualityProfile(fnRs[1:20])
ggsave("geotraces_transect_quality_F.jpeg", qpf, width = 60, height = 40, units = "cm", device = "jpeg")
ggsave("geotraces_transect_quality_R.jpeg", qpr, width = 60, height = 40, units = "cm", device = "jpeg")

##Filtering and trimming
filt_path = file.path(path, "filtered")
filtFs = file.path(filt_path, paste0(f.names, "_F_filt.fastq.gz"))
filtRs = file.path(filt_path, paste0(r.names, "_R_filt.fastq.gz"))

##Primers removed at this stage
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,230), trimLeft=c(20,18),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

#Calculating the error model
errF = learnErrors(filtFs, multithread = TRUE)
errR = learnErrors(filtRs, multithread = TRUE)

#Plotting errors
plot.errF = plotErrors(errF, nominalQ = TRUE) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot.errR = plotErrors(errR, nominalQ = TRUE) + theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Dereplicate sequences
derepFs = derepFastq(filtFs, verbose = TRUE)
derepRs = derepFastq(filtRs, verbose = TRUE)
names(derepFs) = f.names
names(derepRs) = r.names

#DADA2 algorithm on dereplicated sequences
dadaFs = dada(derepFs, err = errF, multithread = TRUE)
dadaRs = dada(derepRs, err = errR, multithread = TRUE)
dadaFs[[1]]
dadaRs[[1]]

#Merge paired end reads
mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Check merger success by sample ([[1]] refers to sample 1)
head(mergers[[1]])

#Generate first sequencing table 
seqtab = makeSequenceTable(mergers)
seqtab.ex = seqtab[, nchar(colnames(seqtab)) %in% seq(252, 256)]

#Remove chimeras 
seqtab.ex.chi = removeBimeraDenovo(seqtab.ex, method = "consensus", 
                                   multithread = T, verbose = T)

#Track sequence loss 
getN = function(x) sum(getUniques(x))
track = cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), 
              rowSums(seqtab.ex), rowSums(seqtab.ex.chi))
colnames(track) = c("input", "filtered", "denoised", "merged", 
                    "tabled", "no chim")
rownames(track) = f.names

###Assign taxonomy 
silva.taxa = assignTaxonomy(seqtab.ex.chi, "~/Desktop/PhD/Data/Sequencing Data/geotraces_transect_run/OG3896/OG3896_fastq/silva_nr_v132_train_set.fa", 
                         multithread = T)
##assign species level taxonomy via SILVA
silva.taxa = addSpecies(silva.taxa, "~/Desktop/PhD/Data/Sequencing Data/geotraces_transect_run/OG3896/OG3896_fastq/silva_species_assignment_v132.fa")

#Phylogenetic tree construction 
source("http://bioconductor.org/biocLite.R")
biocLite("DECIPHER")
install.packages("phangorn")
library(DECIPHER)
library(phangorn)

seqs = getSequences(seqtab.ex.chi)
names(seqs) = seqs 
alignment = AlignSeqs(DNAStringSet(seqs), anchor = NA)

phang.align = phyDat(as(alignment, "matrix"), type = "DNA")
dm = dist.ml(phang.align)
treeNJ = NJ(dm)
fit = pml(treeNJ, data = phang.align)

fitGTR = update(fit, k=4, inv=0.2)
fitGTR = optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace=0))
detach("package:phangorn", unload=TRUE)

#Upload metadata
sd.geo.trans = read.csv("~/Desktop/PhD/Data/Sequencing Data/geotraces_transect_run/OG3896/OG3896_fastq/geotraces_transect_metadata_full.csv", header = T)
dim(sd.geo.trans)
head(sd.geo.trans)

# create a vector for sample names
s.vec = as.vector(1:110)  #number should reflect your total number of samples (note that samples from multiple projects were included here and processed simultaneously before subsetting)
s.nam = cbind("sample_", s.vec)
s.nam = as.data.frame(s.nam)
s.names = paste0(s.nam$V1, s.nam$s.vec)
s.names = as.data.frame(s.names)

# apply sample names to metadata
row.names(sd.geo.trans) = s.names$s.names
sd.geo.trans = as.data.frame(sd.geo.trans)
head(sd.geo.trans)

# apply sample names to sequence table
row.names(seqtab.ex.chi) = s.names$s.names

#Construct Phyloseq Object
geo.trans = phyloseq(tax_table(silva.taxa), otu_table(seqtab.ex.chi, taxa_are_rows = FALSE), sample_data(sd.geo.trans), phy_tree(fitGTR$tree))

geo.trans = subset_samples(geo.trans, Type == "Geo") #select only geotraces transect samples (exclude Munida)
geo.trans

munida.trans = phyloseq(tax_table(silva.taxa), otu_table(seqtab.ex.chi, taxa_are_rows = FALSE), sample_data(sd.geo.trans), phy_tree(fitGTR$tree))
munida.trans = subset_samples(munida.trans, Type == "Munida") #select only munida transect samples


# create vector for ASV names
dim(otu_table(geo.trans))
dim(tax_table(geo.trans))
a.vec = as.vector(1:6258)  #number should reflect your total ASVs
a.nam = cbind("asv_", a.vec)
a.nam = as.data.frame(a.nam)
asv.names = paste0(a.nam$V1, a.nam$a.vec)
asv.names = as.data.frame(asv.names)

taxa_names(geo.trans) = asv.names$asv.names

# Root phylogenetic tree 
pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <- 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>% 
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}

pick_new_outgroup(phy_tree(geo.trans)) #asv_891
rootedTree = ape::root(phy_tree(geo.trans), outgroup="asv_891", resolve.root=TRUE) 
geo.trans = phyloseq(tax_table(geo.trans), otu_table(geo.trans), sample_data(geo.trans), phy_tree(rootedTree))
is.rooted(phy_tree(geo.trans))
phy_tree(geo.trans)


#Analysis ----

#reads remaining per sample
rowSums(otu_table(geo.trans))
mean(rowSums(otu_table(geo.trans))) #57,698.41
min(rowSums(otu_table(geo.trans))) #31,178
max(rowSums(otu_table(geo.trans))) #12,5497
View(tax_table(geo.trans))

# Restructuring tax_table to include 'best' classification
library(zoo)
library(tibble)
bc.t = t(as.data.frame(tax_table(geo.trans)))
bc.t[bc.t=="k__"] <- ""
bc.t[bc.t=="p__"] <- ""
bc.t[bc.t=="c__"] <- ""
bc.t[bc.t=="o__"] <- ""
bc.t[bc.t=="f__"] <- ""
bc.t[bc.t=="g__"] <- ""
bc.t[bc.t=="s__"] <- ""
bc.t[bc.t==""] <- NA
bc.fill = na.locf(bc.t, na.rm = TRUE)
t.bc.fill = as.data.frame(t(bc.fill))
head(t.bc.fill)
rnc.bc = rownames_to_column(t.bc.fill, "ASV")

## Creates a column with the best classification and the ASV
rnc.bc$taxa_ASV = paste(rnc.bc$Species,rnc.bc$ASV)

## Bind this column back onto the original tax_table 
safe.bc = as.data.frame(tax_table(geo.trans))
safe.bc$taxa_ASV = paste(rnc.bc$taxa_ASV)
View(safe.bc)

# Setup object as tax_table
bc.tax = tax_table(safe.bc)
colnames(bc.tax) = colnames(safe.bc)
rownames(bc.tax) = rownames(safe.bc)
View(bc.tax)

## Update phyloseq object with new table
identical(bc.tax[1:6258,1:7], tax_table(geo.trans)) #should be true
tax_table(geo.trans) = bc.tax
View(tax_table(geo.trans))

###Subset to remove chloroplast sequences, and unassigned sequnces (only keep Bacteria and Archaea)
dim(tax_table(geo.trans))
gt.x = subset_taxa(geo.trans, Kingdom=="Bacteria" | Kingdom=="Archaea")
dim(tax_table(gt.x))
gt.x = subset_taxa(gt.x, (Order!="Chloroplast") | is.na(Order))
dim(tax_table(gt.x))
gt.x = subset_taxa(gt.x, (Family!="Mitochondria") | is.na(Family))
dim(tax_table(gt.x))

###Rarefaction Curve
gt.curve = ggrare(gt.x, step = 1000, color = "Type", se = FALSE)
gt.curve = gt.curve + facet_wrap(~Type) + theme_bw()
quartz()
gt.curve
ggsave("geotraces_transect_rarefaction_curve.jpeg", gt.curve, width = 15, height = 7.5, units = "cm", device = "jpeg")

###Rarefaction
rowSums(otu_table(gt.x))
mean(rowSums(otu_table(gt.x))) #53,254.1
min(rowSums(otu_table(gt.x))) #29,280
max(rowSums(otu_table(gt.x))) #114,880


set.seed(711)
gt.rare = rarefy_even_depth(gt.x, sample.size = 30462, trimOTUs = TRUE) 

