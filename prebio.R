# Data analysis for BMT303 trial
# Fiona Tamburini

# required packages
library(ggplot2)
library(genefilter)
library(RColorBrewer)
library(plyr)
library(dplyr)
library(reshape2)
library(scales)
library(MASS)
library(gtools)
library(vegan)
library(ggpubr)
library(cowplot)

######################################################################
### Setup ############################################################
######################################################################

### set this to /your/path/to/prebio2
setwd("/Users/Fiona/scg4_fiona/prebio2/prebio")

# color palette
# FOS, Control
my_pal <- c("#D55E00", "#0072B2")
names(my_pal) <- c("FOS", "Control")

######################################################################
### Read in data and metadate files for prebiotic project analysis ###
######################################################################

# TO DO: change filepaths/organize for portability
# TO DO: remove P83 and re-save

### Read sample metadata -- which stools were collected/sequenced
prebio_meta_all <- read.table("metadata/prebio_meta.tsv", sep = '\t', header = T, quote="\"")

# set FOS/Control grouping
prebio_meta_all$group <- ifelse(startsWith(as.character(prebio_meta_all$patient_id), '303'), "FOS", "Control")
prebio_meta_all$group <- factor(prebio_meta_all$group, levels = c("FOS", "Control"))

# format columns as date
prebio_meta_all$date <- as.Date(prebio_meta_all$date)
prebio_meta_all$trx <- as.Date(prebio_meta_all$trx)

# set factor levels for downstream plots
prebio_meta_all$patient_id <- factor(prebio_meta_all$patient_id, levels = mixedsort(unique(prebio_meta_all$patient_id)))

# metadata for sequenced samples only
prebio_meta <- filter(prebio_meta_all, sequenced_status == T)
prebio_meta <- prebio_meta[mixedorder(unique(prebio_meta$sequencing_id)), ]

### Read taxonomic classification data
## bracken species read counts including unclassifed
brack_sp_reads <- read.table("input_data/bracken_species_reads.txt", sep = '\t', header = T, quote = "")
brack_g_reads <- read.table("input_data/bracken_genus_reads.txt", sep = '\t', header = T, quote = "")

## bracken species percentage -- classified only
brack_sp_perc <- read.table("input_data/bracken_species_perc.txt", sep = '\t', header = T, quote = "")
brack_g_perc <- read.table("input_data/bracken_genus_perc.txt", sep = '\t', header = T, quote = "")

## Read short chain fatty acid measurements
# repeated measurements may 2019
scfa2_f <- "input_data/prebio_scfa_may19.txt"
scfa2 <- read.table(scfa2_f, sep = '\t', header = T)
scfa2[is.na(scfa2)] <- 0

######################################################################
### Summary statistics ###############################################
######################################################################

# n patients, controls
print("FOS")
length(unique(filter(prebio_meta_all, group == "FOS")$patient_id))

print("Controls")
length(unique(filter(prebio_meta_all, group == "Control")$patient_id))

# n samples collected
length(prebio_meta_all$sequencing_id[!is.na(prebio_meta_all$sequencing_id)])

# n samples sequenced
length(prebio_meta_all$sequencing_id[prebio_meta_all$sequenced_status])

# samples collected but not sequenced
not_seqd <- filter(prebio_meta_all, !sequenced_status)

# samples collected per patient
all_freq <- plyr::count(prebio_meta_all[!is.na(prebio_meta_all$sequencing_id),], "patient_id")
fos_freq <- plyr::count(filter(prebio_meta_all[!is.na(prebio_meta_all$sequencing_id),], group == "FOS"), "patient_id")
ctrl_freq <- plyr::count(filter(prebio_meta_all[!is.na(prebio_meta_all$sequencing_id),], group == "Control"), "patient_id")

# median samples collected per patient
median(all_freq$freq)
median(fos_freq$freq)
median(ctrl_freq$freq)

# mean samples collected per patient
mean(all_freq$freq)
mean(fos_freq$freq)
mean(ctrl_freq$freq)

# range
range(all_freq$freq)
range(fos_freq$freq)
range(ctrl_freq$freq)

# samples not sequenced
filter(prebio_meta_all, sequenced_status == F & !is.na(date))

# samples not collected
filter(prebio_meta_all, sequenced_status == F & is.na(date))


######################################################################
### Readcount plots ##################################################
######################################################################

# readcounts file from preprocessing pipeline
readcounts_f <- "input_data/readcounts.tsv"
readcounts <- read.table(readcounts_f, sep = '\t', header = T)
counts <- readcounts[, c(1:3, 5, 7)]
colnames(counts) <- c("Sample", "Raw reads", "Trimmed reads", "Deduplicated reads", "Non-human reads")
counts_long <- melt(counts, id.vars = "Sample", variable.name = "step", value.name = "reads")
counts_long$reads_m <- (counts_long$reads / 1e6)

# plot readcounts
readcount_plot <- ggplot(counts_long, aes(x=reads_m, fill=step)) +
  geom_histogram(binwidth = 1) +
  scale_x_continuous(labels = comma, breaks = seq(0, 100, 10)) +
  facet_grid(step ~ ., scales = "free_y") +
  theme_cowplot(12) +
  labs(
    x = "\nReads (M)",
    y = "Count\n",
    fill = ""
  ) +
  background_grid()

ggsave("plots/readcounts_preproccessing.png", readcount_plot, device = "png", height = 6, width = 7)


######################################################################
### Sample collection plot ###########################################
######################################################################

# plot relative to date of transplant
samples <- prebio_meta_all
samples$sample_day <- (samples$date - samples$trx)

# create patient labels, set order
fos <- filter(samples, group == "FOS")
control <- filter(samples, group == "Control")
labels <- data.frame(patient_id = sort(unique(fos$patient_id)), label = paste0("F", seq(unique(fos$patient_id))))
labels <- rbind(labels, data.frame(patient_id = mixedsort(as.character(unique(control$patient_id))), label = paste0("C", seq(unique(control$patient_id)))))

samples <- merge(samples, labels, by = "patient_id", all = T)
samples$label <- factor(samples$label, levels = rev(labels$label))

# set sequenced vs no
samples$sequenced_status <- ifelse(samples$sequenced_status, "Sequenced", "Not sequenced")
samples$sequenced_status <- ifelse(is.na(samples$sequencing_id), "Not collected", samples$sequenced_status)
samples$sequenced_status <- factor(samples$sequenced_status, levels = c("Sequenced", "Not sequenced", "Not collected"))

# if the sample wasn't collected and the day is NA, change sample day to actual day
samples$sample_day <- ifelse(is.na(samples$sample_day), samples$day, samples$sample_day)

# remove samples > day 100
# maybe change this so that samples >100 are included and axis is >100 ?
samples <- filter(samples, sample_day <= 100)

# plot collected samples
sample_plot <- ggplot(samples, aes(x=sample_day, y=label, shape=sequenced_status)) +
  geom_point(size = 2, color = "black") +
  scale_shape_manual(values = c(16, 1, 4)) +
  facet_wrap(~ group, ncol = 1, strip.position = "top", scales = "free_y") +
  theme_cowplot() +
  labs(
    x = "\nDay relative to transplant",
    y = "Patient\n",
    shape = "Status"
  ) +
  scale_x_continuous(labels = comma, breaks = c(-5, 0, 7, 14, 28, 60, 100))

ggsave("plots/stool_sampling.png", sample_plot, device = "png", height = 6, width = 6)

# color by timepoint
sample_plot2 <- ggplot(samples, aes(x=sample_day, y=label, shape=sequenced_status)) +
  geom_point(size = 2, aes(color = factor(day, levels = c(-5, 0, 7, 14, 28, 60, 100)))) +
  scale_shape_manual(values = c(16, 1, 4)) +
  facet_wrap(~ group, ncol = 1, strip.position = "top", scales = "free_y") +
  theme_cowplot(12) +
  labs(
    x = "\nDay relative to transplant",
    y = "Patient\n",
    shape = "Status",
    color = "Timepoint"
  ) +
  scale_x_continuous(labels = comma, breaks = c(-5, 0, 7, 14, 28, 60, 100))

ggsave("plots/stool_sampling_colored.png", sample_plot2, device = "png", height = 6, width = 6)


######################################################################
### SCFA measurements ################################################
######################################################################

## repeated measurements may 2019
scfa_long2 <- melt(scfa2, id.vars = c("sample", "patient_id", "sequencing_id", "group"), variable.name = "scfa")
scfa_long2$scfa <- gsub("\\.A", " a", scfa_long2$scfa)

# set factor level for group
scfa_long2$group <- factor(scfa_long2$group, levels = c("FOS", "Control"))

## plot without log transformation, free y axis
pvals <- compare_means(value ~ group, data = scfa_long2, group.by = "scfa", method = "wilcox.test", p.adjust.method = "fdr")
pvals$p.signif <- ifelse(pvals$p.adj < 0.05, "*", "ns")
pvals$p.signif <- ifelse(pvals$p.adj < 0.01 & pvals$p.adj >= 0.001, "**", pvals$p.signif)
pvals$p.signif <- ifelse(pvals$p.adj < 0.001, "***", pvals$p.signif)

# set y position of signif for each plot
maxs <- aggregate(value ~ scfa,scfa_long2, FUN = max)
pvals$y.position <- maxs[match(pvals$scfa, maxs$scfa), "value"] * 1.10

scfa_plot <- ggplot(scfa_long2, aes(x = group, y = value)) + 
  geom_violin(aes(fill = group)) +
  geom_point() +
  facet_wrap(. ~ scfa, scales = "free_y") +
  # pseudo_log_trans() +
  labs(
    x = "Short-chain fatty acid",
    y = "Concentration (umol/g stool)",
    fill="") +
  scale_fill_manual(values = my_pal) +
  stat_pvalue_manual(pvals, label = "p.signif") +
  theme_cowplot(12)

ggsave("plots/scfa_may19_facet.png", scfa_plot, device = "png", height = 9, width = 8)


######################################################################
### Classified reads #################################################
######################################################################

## Plot histogram of classified reads
classified <- (1 - sweep(brack_sp_reads, 2, colSums(brack_sp_reads), "/")["Unclassified",]) * 100

read_plot <- ggplot(melt(classified), aes(x=value)) +
  geom_histogram(binwidth = 1, fill = "cornflowerblue", color = "white") +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  theme_cowplot(12) +
  scale_fill_manual(values = my_pal) +
  labs(
    x = "Percentage of reads classified",
    y = "Count"
  )

ggsave("plots/readcounts_classified_histo.png", read_plot, device = "png", height = 4, width = 5)

######################################################################
### Diversity plots ##################################################
######################################################################

# find shannon diversity with vegdist
shannon_div <- diversity(t(brack_sp_perc), index = "shannon")
div <- data.frame("shannon_div" = shannon_div, "sequencing_id" = names(shannon_div))
div_meta <- merge(div, prebio_meta, by = "sequencing_id")

## stat smooth shannon diversity over time
shannon_plot_smooth <- ggplot(div_meta, aes(day, shannon_div, color = group)) +
  geom_point() +
  stat_smooth() +
  labs(
       x = "Day",
       y = "Shannon Diversity",
       color="") +
  theme_cowplot(12) +
  scale_color_manual(values = my_pal) +
  scale_x_continuous(labels = comma, breaks = c(-5, 0, 7, 14, 28, 60, 100))

ggsave("plots/shannon_line_smooth.png", shannon_plot_smooth, device = "png", height = 4, width = 6)


## violin plot -- alpha diversity at each timepoint

## compare means
pvals <- compare_means(shannon_div ~ group, data = div_meta, group.by = "day", method = "wilcox.test", p.adjust.method = "fdr")
pvals$p.signif <- ifelse(pvals$p.adj < 0.05, "*", "ns")
pvals$p.signif <- ifelse(pvals$p.adj < 0.01 & pvals$p.adj >= 0.001, "**", pvals$p.signif)
pvals$p.signif <- ifelse(pvals$p.adj < 0.001, "***", pvals$p.signif)
pvals$y.position <- 8

# plot
shannon_plot <- ggplot(div_meta, aes(x=group, y=shannon_div)) + 
  geom_violin(aes(fill = group), position=position_dodge(.9), trim = F) +
  stat_summary(fun.data=mean_sdl, aes(group=group), position=position_dodge(.9), geom="pointrange", color="black") +
  facet_grid(. ~ day, scales = "free") +
  labs(
       x = "\nTreatment",
       y = "Shannon Diversity\n",
       fill="") +
  theme_cowplot(12) +
  scale_fill_manual(values = my_pal) +
  stat_pvalue_manual(pvals, label = "p.signif")

ggsave("plots/shannon_div.png", shannon_plot, device = "png", height = 4, width = 10)


######################################################################
### NMDS ordination ##################################################
######################################################################

### ordinate species-level classifications

### find pairwise bray-curtis distances with vegdist
vare_dis <- vegdist(t(brack_sp_perc), method = "bray")

### nmds ordinate
vare_mds0 <- isoMDS(vare_dis)
mds <- data.frame(vare_mds0$points)
mds$sequencing_id <- row.names(mds)

### merge pheno data
mds_meta <- merge(mds, prebio_meta, by = "sequencing_id")

### function to create scatterplot
nmds_plot <- ggplot(mds_meta, aes(x = X1, y = X2, color = group)) +
  geom_point(size = 2) +
  theme_cowplot(12) +
  scale_color_manual(values = my_pal) +
  labs(
       x = "NMDS1",
       y = "NMDS2",
       color = ""
       )

# add 95% confidence ellipse
nmds_plot_ci <- nmds_plot + stat_ellipse(type = 't', size = 1)

ggsave("plots/nmds_by_treatment_ci.png", nmds_plot_ci, device = "png", height = 5, width = 6)

# test group differences
# beta dispersions -- are assumptions for PERMANOVA met?
dispersion <- betadisper(vare_dis, group = prebio_meta$group)
permutest(dispersion)

adonis(vare_dis ~ group, data = filter(prebio_meta, sequenced_status == T))


######################################################################
### Input tables for lefse ###########################################
######################################################################

dir.create("lefse")

# function to keep features that are at least A relative abundance and p prevalence
subset_lefse <- function(bracken_data, filt_day, relab, prop, rank){
  
  # filter metadata
  lefse_meta <- filter(prebio_meta, sequenced_status == T, day == filt_day)[, c("sequencing_id", "group")]
  lefse_meta <- lefse_meta[order(as.character(lefse_meta$sequencing_id)), ]
  lefse_meta_t <- t(lefse_meta)
  
  # filter taxa
  tax <- bracken_data[, sort(as.character(lefse_meta$sequencing_id))]
  rownames(tax) <- gsub(' ', '_', rownames(tax))
  
  # remove rows that sum to zero
  tax <- tax[rowSums(tax) > 0, ]
  
  keep <- data.frame(genefilter(tax, pOverA(p=prop, A=relab * 1e6)))
  colnames(keep) <- "taxon"
  keep$tax <- row.names(keep)
  keep <- filter(keep, taxon == T)$tax
  tax_filt <- tax[keep, ]
  
  fname <- paste0("lefse/lefse_input_relab", relab, "_p", prop, "_", rank, ".txt")
  write.table(lefse_meta_t, fname, sep = '\t', row.names = T, col.names = F, quote = F)
  write.table(tax_filt, fname, sep = '\t', row.names = T, col.names = F, quote = F, append = T)
  
  # print(F %in% (colnames(tax) == lefse_meta_t[1,]))
}

# 0.01% relative abundance, 10%
subset_lefse(brack_sp_perc * 1e6, 14, 0.01, 0.10, "sp")
subset_lefse(brack_g_perc * 1e6, 14, 0.01, 0.10, "g")

# no filtering
subset_lefse(brack_sp_perc * 1e6, 14, 0, 0, "sp")
subset_lefse(brack_g_perc * 1e6, 14, 0, 0, "g")

# next, run lefse on the Huttenhower lab galaxy server (https://huttenhower.sph.harvard.edu/galaxy/)
# or on the command line

######################################################################
### Taxonomy area plots ##############################################
######################################################################

dir.create("plots/area_plots", showWarnings = F)

## species
sp_data <- brack_sp_perc
sp_data$taxon <- row.names(sp_data)
sp_long <- melt(sp_data, id.vars = "taxon", variable.name = "sequencing_id", value.name = "rel_abundance")
sp_long_meta <- merge(sp_long, prebio_meta, by = "sequencing_id")
patient_list <- unique(sp_long_meta$patient_id)

for (patient in patient_list) {
  plot_data <- filter(sp_long_meta, patient_id == patient)
  
  # plot only n top taxa
  n_taxa <- 20
  
  # color palette for n taxa
  myCols <- colorRampPalette(brewer.pal(12, "Paired"))
  my_pal <- myCols(n_taxa)
  my_pal <- sample(my_pal)
  
  tax <- aggregate(rel_abundance ~ taxon, data = plot_data, sum)
  tax <- tax[rev(order(tax$rel_abundance)), ]
  top_taxa <- tax[1:n_taxa, "taxon"]
  
  plot_filt <- filter(plot_data, taxon %in% top_taxa)
  
  area_plot <- ggplot(plot_filt, aes(day, rel_abundance * 100, group = taxon)) +
    geom_area(aes(fill = taxon)) +
    labs(
      title=paste("Patient", patient),
      x = "Day",
      y = "Species Relative Abundance",
      fill="Species") +
    scale_fill_manual(values=my_pal, guide = guide_legend(ncol = 1)) +
    scale_x_continuous(breaks = c(-5, 0, 7, 14, 28, 60, 100)) +
    scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
    theme_cowplot(12)
  
  ggsave(paste0("plots/area_plots/", patient, "_species.png"), area_plot, device = "png", height = 6, width = 10)
}

## genus
g_data <- brack_g_perc
g_data$taxon <- row.names(g_data)
g_long <- melt(g_data, id.vars = "taxon", variable.name = "sequencing_id", value.name = "rel_abundance")
g_long_meta <- merge(g_long, prebio_meta, by = "sequencing_id")
patient_list <- unique(g_long_meta$patient_id)

for (patient in patient_list) {
  plot_data <- filter(g_long_meta, patient_id == patient)
  
  # plot only n top taxa
  n_taxa <- 20
  
  # color palette for n taxa
  myCols <- colorRampPalette(brewer.pal(12, "Paired"))
  my_pal <- myCols(n_taxa)
  my_pal <- sample(my_pal)
  
  tax <- aggregate(rel_abundance ~ taxon, data = plot_data, sum)
  tax <- tax[rev(order(tax$rel_abundance)), ]
  top_taxa <- tax[1:n_taxa, "taxon"]
  
  plot_filt <- filter(plot_data, taxon %in% top_taxa)
  
  area_plot <- ggplot(plot_filt, aes(day, rel_abundance * 100, group = taxon)) +
    geom_area(aes(fill = taxon)) +
    labs(
      title=paste("Patient", patient),
      x = "Day",
      y = "Species Relative Abundance",
      fill="Species") +
    scale_fill_manual(values=my_pal, guide = guide_legend(ncol = 1)) +
    scale_x_continuous(breaks = c(-5, 0, 7, 14, 28, 60, 100)) +
    scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100)) +
    theme_cowplot(12)
  
  ggsave(paste0("plots/area_plots/", patient, "_genus.png"), area_plot, device = "png", height = 6, width = 10)
}

######################################################################
### Boxplots of specific features ####################################
######################################################################

# plot_data <- brack_g_perc
# plot_data$taxon <- row.names(plot_data)
# data_long <- melt(plot_data, id.vars = "taxon", variable.name = "sequencing_id", value.name = "rel_abundance")
# data_long_meta <- merge(data_long, prebio_meta, by = "sequencing_id")
# 
# taxa <- c("Lactobacillus", "Blautia")
# data_filt <- filter(data_long_meta, taxon %in% taxa, day == 14)
# 
# tax_boxplot <- ggplot(data_filt, aes(x=taxon, y=rel_abundance)) + 
#   geom_boxplot(aes(fill = group), position=position_dodge(.9)) +
#   # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, aes(fill = Treatment), position=position_dodge(.9)) +
#   # stat_summary(fun.data=mean_sdl, mult=1, aes(group=group), position=position_dodge(.9), geom="pointrange", color="black") +
#   facet_wrap(. ~ taxon, scales = "free") +
#   # scale_y_log10() +
#   labs(title='',
#        x = "\nGenus",
#        y = "Relative abundance (%)\n",
#        fill="") +
#   theme_cowplot(12)
# 
# ggsave("plots/lefse_g_boxplot.png", tax_boxplot, device = "png", height = 4, width = 2.5 * length(taxa))
