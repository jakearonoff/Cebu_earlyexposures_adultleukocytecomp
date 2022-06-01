library(dplyr)
library(ggplot2)
library(ewastools)
library(stringi)
library(data.table)
library(magrittr)   
library(purrr)
library(svd)
library(wateRmelon)
library(minfi)


################################################################################
## Flag/remove samples with low overall intensities 
################################################################################

# load in sample sheet and IDATS
targets <- read.metharray.sheet("path_to_sample_sheet")
RGset <- read.metharray.exp(targets = targets, verbose = T, extended = T,
                            force = TRUE)

# extract raw values
rawMSet <- preprocessRaw(RGset)
rawMSet

Meth <- getMeth(rawMSet)
Unmeth <- getUnmeth(rawMSet)

# Get overall intensities for Meth and Unmeth in each sample to find outliers with low overall intensities 
targets$MQC <- log2(colMedians(Meth))
targets$UQC <- log2(colMedians(Unmeth))

# Visualize intensities 
# look for datapoints that stand out for low overall intensity 
ggplot(targets, aes(UQC, MQC)) + geom_point() + 
  coord_cartesian(xlim = c(11,13), ylim = c(11,15)) + 
  labs(x = "Log2 Median Unmethylated Intensity", y = "Log2 Median Methylated Intensity") + 
  theme_bw()


################################################################################
### check if samples passed control metrics with ewastools 
################################################################################

meth <- read_idats(targets$Basename,quiet = FALSE) %>% detectionP %>% correct_dye_bias()
ctrls <- control_metrics(meth)
targets$ctrlfail <- ewastools::sample_failure(ctrls)
table(targets$ctrlfail) # check if any samples failed control metrics 

################################################################################
### checking for mismatches between reported sex and estimated chromosomal sex with ewastools (a way to identify potential samples swaps) 
#      before beginning, make sure the sample sheet (i.e., targets) has a FACTOR variable named "sex" and coded values as "m" and "f" 
################################################################################
targets <- as.data.table(targets)
targets[,c("X","Y") := check_sex(meth)]
targets[,predicted_sex:=predict_sex(X,Y,which(sex=="m"),which(sex=="f"))]

# see the mismatching samples 
targets[sex!=predicted_sex,.(Sample_Name,sex,predicted_sex)]


###################################
# detection P check, a bad sample is identified if more than 10% of its probes have detection P > 0.01 and/or a bead count < 5
###################################
# get beta matrix for later reference
beta <- meth %>% dont_normalize
# get detection p values 
detp <- meth$detP
# give detection p matrix same row/column names as beta matrix
dimnames(detp) <- dimnames(beta)

# get number beads 
mani <- ewastools:::manifest_epic
nbeads <- wateRmelon::beadcount(RGset)

# exclude SNP probes from consideration
nbeads <- nbeads[match(mani[mani$probe_type != "rs", probe_id], rownames(nbeads)),]

# visualize average bead counts per probe 
beadmean <- data.frame(bmean = rowMeans(nbeads, na.rm = T))
ggplot(beadmean, aes(x = bmean)) + 
  geom_histogram(bins = 50, color = "black", fill = "grey") + 
  coord_cartesian(xlim = c(0,25)) + 
  labs(x = "Probe mean beadcount") + 
  theme_bw()

# see % of probes with bead counts under 5
(table(nbeads < 5) / length(unlist(nbeads))) %>% round(digits = 3)

detp <- detp[!grepl("rs", meth$manifest$probe_id),] # removing SNP probes from detection p matrix

## specify bead number min. cutoff of 5, and detection p max cutoff of 0.01
beadmin <- 5
detpmax <- 0.01

reliable <- (detp > detpmax) | (nbeads < beadmin)  # creates a matrix of TRUE of FALSE results, in which TRUE indicates a bad probe 
reliable[is.na(reliable)] <- F  # fill in NA's with FALSE


# filtering samples 
# this is the sample threshold, specifying that a sample cannot have more than 10% of failed probes 
sthresh <- 0.1

# columns are samples in the matrix, so the column means reflect the % of failed probes in each sample 
smeans <- colMeans(reliable)
sum(smeans > sthresh) # look at % of bad samples 

# selecting bad samples from the detection p matrix 
detp <- detp[, smeans > sthresh] 
# shows the names of failed samples for later removal 
colnames(detp)




###################################
# Contamination Check 
###################################

# get names of SNP probes 
snps = meth$manifest[probe_type=="rs",index]

# extract betas for SNP probes 
snps = beta[snps,]

# call genotypes for SNP probes 
genotypes = call_genotypes(snps,learn=FALSE)

# create outlier values and match up with sample names in targets 
targets$outlier = snp_outliers(genotypes)

# identify contaminated samples using outlier values 
contam <- targets[targets$outlier > -4,c("Sample_Name")]
contam # samples that were likely contaminated 


###################################
# Additional step in ewastools that uses SNP probes to create unique genotype ID's. If samples have same genotype ID, this can flag unidentified technical replicates 
###################################
targets$sample_geno_id = enumerate_sample_donors(genotypes)

# List duplicates
targets[,n:=.N,by=sample_geno_id]
targets[n>1,.(sample_name,sample_geno_id)]




###################################
# Cell deconvolution using samples that passed QC
###################################

# here I have already stored the samples that failed QC into the object called "badsamples" to filter them out prior to cell deconvolution 
targets <- filter(targets, !(Sample_Name %in% badsamples))
RGset <- read.metharray.exp(targets = targets, verbose = T,
                            force = TRUE)

# cell deconvolution using Salas et al. 2018 reference, also noob normalizing 
noob <- preprocessNoob(RGset, dyeCorr = TRUE, verbose = TRUE,
               dyeMethod="single")
beta <- getBeta(noob)
##  all 450 IDOL probes in noob beta matrix 
cells <- estimateLC(beta, ref = "Salas", constrained = T)












