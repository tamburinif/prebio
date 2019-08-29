# Data analysis for BMT303 trial
# Fiona Tamburini

# required packages
library(ggplot2)
# # library(ggrepel)
# # library(RColorBrewer)
# library(plyr)
library(dplyr)
# library(zoo)
library(reshape2)
# library(scales)
library(gtools)
# library(vegan)
# library(MASS)
# library(ggpubr)
# library(ggbeeswarm)
# library(cowplot)

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
brack_sp <- read.table("input_data/bracken_species_reads.txt", sep = '\t', header = T, quote = "")
brack_g <- read.table("input_data/bracken_genus_reads.txt", sep = '\t', header = T, quote = "")

# relative abundance normalize -- classified/no human only
remove <- c("Unclassified", "Homo sapiens")
brack_sp_filt <- brack_sp[!rownames(brack_sp) %in% remove,]
brack_g_filt <- brack_g[!rownames(brack_g) %in% remove,]
brack_sp_perc <- sweep(brack_sp_filt, 2, colSums(brack_sp_filt), "/")
brack_g_perc <- sweep(brack_g_filt, 2, colSums(brack_g_filt), "/")

# remove P83 sample from all taxonomy tables
brack_sp_perc[, "P83"] <- NULL
brack_g_perc[, "P83"] <- NULL

# set order
brack_sp_perc <- brack_sp_perc[, mixedorder(unique(prebio_meta$sequencing_id))]
brack_g_perc <- brack_g_perc[, mixedorder(unique(prebio_meta$sequencing_id))]

write.table(brack_sp_perc, "input_data/bracken_species_perc.txt", sep = "\t", row.names = F, quote = F)
write.table(brack_g_perc, "input_data/bracken_genus_perc.txt", sep = "\t", row.names = F, quote = F)

## Read short chain fatty acid measurements
# initial
scfa1_f <- "/Users/Fiona/Google Drive/Bhatt lab/Projects/Prebiotic/Prebio2/prebio_scfa_initial.txt"
scfa1 <- read.table(scfa1_f, sep = '\t', header = T)
scfa1[is.na(scfa1)] <- 0

# repeated may 2019
scfa2_f <- "/Users/Fiona/Google Drive/Bhatt lab/Projects/Prebiotic/Prebio2/prebio_scfa_may19.txt"
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
# write.table(not_seqd, "/Users/Fiona/Desktop/prebio_not_sequenced.tsv", sep = '\t', quote = F, row.names = F)

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
readcounts_f <- "01_processing/readcounts.tsv"
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
pvals$y.position <- maxs[match(pvals$scfa, maxs$scfa), "value"]*1.10

scfa_plot <- ggplot(scfa_long2, aes(x = group, y = value)) + 
  # geom_violin(aes(fill = group), position=position_dodge(.9)) +
  geom_violin(aes(fill = group)) +
  geom_point() +
  facet_wrap(. ~ scfa, scales = "free_y") +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, aes(fill = group), position=position_dodge(.9)) +
  # stat_summary(fun.data=mean_sdl, mult=1, aes(group=group), position=position_dodge(.9), geom="pointrange", color="black") +
  # scale_y_log10() +
  # pseudo_log_trans() +
  labs(
    x = "Short-chain fatty acid",
    y = "Concentration (umol/g stool)",
    fill="") +
  theme_cowplot(12) +
  scale_fill_manual(values = my_pal) +
  # my_thm +
  # theme(
  #   axis.text.x = element_text(size = 18, angle = 20, hjust = 1),
  #   strip.text = element_text(size = 16)
  # ) +
  stat_pvalue_manual(pvals, label = "p.signif")

ggsave("plots/scfa_may19_facet.png", scfa_plot, device = "png", height = 9, width = 8)


######################################################################
### Classified reads #################################################
######################################################################

## Plot histogram of classified reads
classified <- (1 - sweep(brack_sp, 2, colSums(brack_sp), "/")["Unclassified",]) * 100

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

## shannon diversity over time

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


## Violin plot -- alpha diversity each day

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
### Taxonomy stream plots ############################################
######################################################################

## species
plot_data <- brack_sp_perc
plot_data$taxon <- row.names(plot_data)
data_long <- melt(plot_data, id.vars = "taxon", variable.name = "sequencing_id", value.name = "rel_abundance")
data_long_meta <- merge(data_long, prebio_meta, by = "sequencing_id")

patient_list <- unique(data_long_meta$patient_id)

for (patient in patient_list) {
  plot_data <- filter(data_long_meta, patient_id == patient)
  
  # plot only n top taxa
  n_taxa <- 20
  
  # color palette for n taxa
  myCols <- colorRampPalette(brewer.pal(12, "Paired"))
  my_pal <- myCols(n_taxa)
  my_pal <- sample(my_pal)
  # my_pal[n_taxa + 1] <- "gray"
  
  # find top n taxa
  # from Ben
  # abundance_threshold <- sort(rowSums(plot_data), decreasing = T)[n_taxa]
  # bracken_plot <- plot_data[rowSums(plot_data) >= abundance_threshold,]
  
  tax <- aggregate(rel_abundance ~ taxon, data = plot_data, sum)
  tax <- tax[rev(order(tax$rel_abundance)), ]
  top_taxa <- tax[1:n_taxa, "taxon"]
  
  plot_filt <- filter(plot_data, taxon %in% top_taxa)
  
  area_plot <- ggplot(plot_filt, aes(day, rel_abundance, group = taxon)) +
    geom_area(aes(fill = taxon)) +
    labs(title=paste("Patient", patient),
         subtitle='',
         x = "Day",
         y = "Species Relative Abundance",
         fill="Species") +
    scale_fill_manual(values=my_pal, guide = guide_legend(ncol = 1)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 18, hjust = 1),
      axis.text.y = element_text(size = 18),
      plot.title = element_text(size = 22)
    ) +
    scale_x_continuous(breaks = c(-5, 0, 7, 14, 28, 60, 100)) +
    scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100))
  
  ggsave(paste0("plots/area plots/", patient, ".png"), area_plot, device = "png", height = 8, width = 12)
}

## genus
plot_data <- brack_g_perc
plot_data$taxon <- row.names(plot_data)
data_long <- melt(plot_data, id.vars = "taxon", variable.name = "sequencing_id", value.name = "rel_abundance")
data_long_meta <- merge(data_long, prebio_meta, by = "sequencing_id")

patient_list <- unique(data_long_meta$patient_id)

for (patient in patient_list) {
  plot_data <- filter(data_long_meta, patient_id == patient)
  
  # plot only n top taxa
  n_taxa <- 20
  
  # color palette for n taxa
  myCols <- colorRampPalette(brewer.pal(12, "Paired"))
  my_pal <- myCols(n_taxa)
  my_pal <- sample(my_pal)
  # my_pal[n_taxa + 1] <- "gray"
  
  # find top n taxa
  # from Ben
  # abundance_threshold <- sort(rowSums(plot_data), decreasing = T)[n_taxa]
  # bracken_plot <- plot_data[rowSums(plot_data) >= abundance_threshold,]
  
  tax <- aggregate(rel_abundance ~ taxon, data = plot_data, sum)
  tax <- tax[rev(order(tax$rel_abundance)), ]
  top_taxa <- tax[1:n_taxa, "taxon"]
  
  plot_filt <- filter(plot_data, taxon %in% top_taxa)
  
  area_plot <- ggplot(plot_filt, aes(day, rel_abundance, group = taxon)) +
    geom_area(aes(fill = taxon)) +
    labs(title=paste("Patient", patient),
         subtitle='',
         x = "Day",
         y = "Species Relative Abundance",
         fill="Species") +
    scale_fill_manual(values=my_pal, guide = guide_legend(ncol = 1)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 18, hjust = 1),
      axis.text.y = element_text(size = 18),
      plot.title = element_text(size = 22)
    ) +
    scale_x_continuous(breaks = c(-5, 0, 7, 14, 28, 60, 100)) +
    scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100))
  
  ggsave(paste0("plots/area plots/", patient, "_genus.png"), area_plot, device = "png", height = 8, width = 12)
}


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
### Boxplots of specific features ####################################
######################################################################

plot_data <- brack_g_perc
plot_data$taxon <- row.names(plot_data)
data_long <- melt(plot_data, id.vars = "taxon", variable.name = "sequencing_id", value.name = "rel_abundance")
data_long_meta <- merge(data_long, prebio_meta, by = "sequencing_id")

# taxa <- c("Bacteroides", "Enterococcus")
taxa <- c("Lactobacillus", "Blautia", "Ralstonia")
data_filt <- filter(data_long_meta, taxon %in% taxa, day == 14)

g <- ggplot(data_filt, aes(x=taxon, y=rel_abundance)) + 
  geom_boxplot(aes(fill = group), position=position_dodge(.9)) +
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, aes(fill = Treatment), position=position_dodge(.9)) +
  # stat_summary(fun.data=mean_sdl, mult=1, aes(group=group), position=position_dodge(.9), geom="pointrange", color="black") +
  facet_wrap(. ~ taxon, scales = "free") +
  # scale_y_log10() +
  labs(title='',
       x = "\nGenus",
       y = "Relative abundance (%)\n",
       fill="") +
  theme_bw() +
  my_thm +
  theme(
    axis.text.x = element_text(size = 18, angle = 20, hjust = 1),
    strip.text = element_text(size = 18)
  )
# ggsave("plots/lefse_g_boxplot.png", g, device = "png", height = 10, width = 10)

######################################################################
### Input tables for lefse ###########################################
######################################################################

# function to keep features that are at least x relative abundance in at least n% samples
subset_lefse <- function(tax_data, this_day, relab, prop, fname){
  # filter metadata
  lefse_meta <- filter(prebio_meta, sequenced_status == T, day == this_day)[, c("sequencing_id", "group")]
  lefse_meta <- lefse_meta[order(as.character(lefse_meta$sequencing_id)), ]
  lefse_meta_t <- t(lefse_meta)
  
  # filter taxa
  tax <- tax_data[, sort(as.character(lefse_meta$sequencing_id))]
  rownames(tax) <- gsub(' ', '_', rownames(tax))
  
  # remove rows that sum to zero
  tax <- tax[rowSums(tax) > 0, ]
  
  # n_samples <- length(unique(lefse_meta$sequencing_id)) * prev
  # tax_long <- melt(data.frame(tax, row.names(tax)), id.vars = "row.names.tax.", variable.name = "sample", value.name = "perc")
  # tax_counts <- plyr::count(filter(tax_long, perc >= relab), "row.names.tax.")
  # taxa_names <- filter(tax_counts, freq >= n_samples)$row.names.tax.
  # tax <- subset(tax, rownames(tax) %in% taxa_names)
  
  keep <- data.frame(genefilter(tax, pOverA(p=prop, A=relab)))
  colnames(keep) <- "taxon"
  keep$tax <- row.names(keep)
  keep <- filter(keep, taxon == T)$tax
  tax_filt <- tax[keep, ]
  
  write.table(lefse_meta_t, fname, sep = '\t', row.names = T, col.names = F, quote = F)
  write.table(tax_filt, fname, sep = '\t', row.names = T, col.names = F, quote = F, append = T)
  
  print(F %in% (colnames(tax) == lefse_meta_t[1,]))
}

# no filtering
# subset_lefse(brack_sp_perc, 14, 0, 0, "03_lefse/lefse_sp_unfiltered_day14.tsv")
# subset_lefse(brack_g_perc, 14, 0, 0, "03_lefse/lefse_g_unfiltered_day14.tsv")

# for (day in c(-5, 0, 7, 14, 28, 60, 100)){
#   print(day)
#   subset_lefse(brack_sp_perc, day, 0, 0, paste0("03_lefse/lefse_sp_unfiltered_day", day, ".tsv"))
#   subset_lefse(brack_g_perc, day, 0, 0, paste0("03_lefse/lefse_g_unfiltered_day", day, ".tsv"))
# }

# filter 0.1% relative abundance, at least 1 sample
# subset_lefse(brack_sp_perc, 14, 0.1, 1, "03_lefse/lefse_sp_relab0.1_n1_day14.tsv")
# subset_lefse(brack_g_perc, 14, 0.1, 1, "03_lefse/lefse_g_relab0.1_n1_day14.tsv")

# filter 1% relative abundance, at least 1 sample
# subset_lefse(brack_sp_perc, 14, 1, 1, "03_lefse/lefse_sp_relab1_n1_day14.tsv")
# subset_lefse(brack_g_perc, 14, 1, 1, "03_lefse/lefse_g_relab1_n1_day14.tsv")

# filter 0.1% relative abundance, at least 10% of samples
# subset_lefse(brack_sp_perc, 14, 0.1, 0.1, "03_lefse/lefse_sp_relab1_prev0.1_day14.tsv")
# subset_lefse(brack_g_perc, 14, 0.1, 0.1, "03_lefse/lefse_g_relab1_prev0.1_day14.tsv")

# filter 0.1% relative abundance, at least 2 samples
subset_lefse(brack_sp_perc, 14, 0.1, 2/29, "03_lefse/lefse_sp_relab0.1_prev2_day14.tsv")
subset_lefse(brack_g_perc, 14, 0.1, 2/29, "03_lefse/lefse_g_relab0.1_prev2_day14.tsv")

######################################################################
### Process/plot lefse results #######################################
######################################################################

## plot lefse barplot results
cols <- c("feature", "log_class_avg", "class", "lda_effect_size", "p_value")
lda <- 2

# function to plot lefse results
read_lefse <- function(fname){
  lefse_res <- read.table(fname, sep = '\t')
  colnames(lefse_res) <- cols
  
  lefse_res <- filter(lefse_res, lda_effect_size > lda)
  lefse_res$feature <- gsub('_', ' ', lefse_res$feature)
  lefse_res$feature <- gsub(' sp ', ' sp. ', lefse_res$feature)
  
  # flip control axis
  lefse_res$lda_effect_size <- ifelse(lefse_res$class == "Control", -1 * lefse_res$lda_effect_size, lefse_res$lda_effect_size)
  lefse_res <- lefse_res[order(lefse_res$lda_effect_size), ]
  lefse_res$feature <- factor(lefse_res$feature, levels = lefse_res$feature)
  lefse_res
}

plot_lefse <- function(lefse_res) {
  lefse_plot <- ggplot(lefse_res, aes(x=feature, y=lda_effect_size, fill=class)) + 
    geom_bar(stat='identity', aes(fill=class), width = 0.75)  +
    # scale_fill_manual(name="Mileage", 
    #                   labels = c("Above Average", "Below Average"), 
    #                   values = c("above"="#00ba38", "below"="#f8766d")) + 
    labs(
      y = "LDA Effect Size",
      x = "Feature"
    ) + 
    coord_flip() +
    theme_bw() +
    scale_fill_manual(values = my_pal) +
    scale_y_continuous(breaks = seq(-4, 4, 1)) +
    my_thm +
    theme(
      axis.text.y = element_text(size = 16, face = "italic"),
      axis.text.x = element_text(size = 16)
    )
  lefse_plot
}

ggsave("plots/lefse_sp_results_filtered.png", plot_lefse(read_lefse("03_lefse/lefse_sp_relab0.1_prev2_day14/lefse_sp_relab0.1_prev2_day14.res")), device = "png", height = 6, width = 12)
ggsave("plots/lefse_g_results_filtered.png", plot_lefse(read_lefse("03_lefse/lefse_g_relab0.1_prev2_day14/lefse_g_relab0.1_prev2_day14.res")), device = "png", height = 3, width = 12)

# ggsave("plots/lefse_sp_results_unfiltered.png", plot_lefse(read_lefse("sp", "unfiltered")), device = "png", height = 6, width = 12)
# ggsave("plots/lefse_g_results_unfiltered.png", plot_lefse(read_lefse("g", "unfiltered")), device = "png", height = 14, width = 12)

## plot significant features from lefse

# plot function
lefse_boxplot <- function(level, taxa_df, filter){
  plot_data <- taxa_df
  plot_data$taxon <- row.names(plot_data)
  data_long <- melt(plot_data, id.vars = "taxon", variable.name = "sequencing_id", value.name = "rel_abundance")
  data_long_meta <- merge(data_long, prebio_meta, by = "sequencing_id")
  
  taxa <- read_lefse(level, filter)$feature
  data_filt <- filter(data_long_meta, taxon %in% taxa, day == 14)
  
  lefse_boxplot <- ggplot(data_filt, aes(x=group, y=rel_abundance)) + 
    geom_boxplot(aes(fill = group), position=position_dodge(.9)) +
    # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, aes(fill = Treatment), position=position_dodge(.9)) +
    # stat_summary(fun.data=mean_sdl, mult=1, aes(group=group), position=position_dodge(.9), geom="pointrange", color="black") +
    facet_wrap(. ~ taxon, scales = "free") +
    # scale_y_log10() +
    labs(title='',
         x = "\nGroup",
         y = "Relative abundance (%)\n",
         fill="") +
    theme_bw() +
    # my_thm +
    theme(
      axis.text.x = element_text(size = 12, angle = 20, hjust = 1),
      strip.text = element_text(size = 12)
    )
  lefse_boxplot
}

ggsave("plots/lefse_boxplot_g_filtered.png", lefse_boxplot("g", brack_g_perc, "filtered"), device = "png", height = 10, width = 10)
ggsave("plots/lefse_boxplot_sp_filtered.png", lefse_boxplot("sp", brack_sp_perc, "filtered"), device = "png", height = 10, width = 10)

ggsave("plots/lefse_boxplot_g_unfiltered.png", lefse_boxplot("g", brack_g_perc, "unfiltered"), device = "png", height = 10, width = 10)
ggsave("plots/lefse_boxplot_sp_unfiltered.png", lefse_boxplot("sp", brack_sp_perc, "unfiltered"), device = "png", height = 30, width = 30)

