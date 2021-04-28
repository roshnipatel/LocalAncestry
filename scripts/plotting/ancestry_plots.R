library(tidyverse)
library(ggpubr)

############################## For MVP ##############################
setwd("~/Google Drive/stanford/writing/admixed PGS paper/figures/schematic/")
MVP_anc = read_tsv("mvp_aa_global_afr_anc.txt")
pdf("mvp_global_ancestry.pdf")
ggplot(MVP_anc, aes(AFR)) + geom_histogram(bins=50, fill="#e48671") + theme_pubr() +
  theme(text = element_text(size=20)) + xlim(0, 1) + 
  xlab("African (YRI) global ancestry fraction")
dev.off()

############################## For MESA ##############################
setwd("~/sherlock/oak/LocalAncestry/results/v2")
dir.create("plots", showWarnings = FALSE)
comb_frac = read_tsv("combined_global_anc_frac.txt") %>% rename(nwd_id = ID)
eur_tracts = read_tsv("european_tract_lengths.txt")

# Remove excluded individuals from dataset
exclusion_idv = read_tsv("individuals_to_exclude.txt", col_names=c("nwd_id"))
comb_frac = comb_frac %>% anti_join(exclusion_idv)
eur_tracts = eur_tracts %>% anti_join(exclusion_idv)

# Plot correlation between RFMix and ADMIXTURE global ancestry estimates
pdf("plots/admixture_rfmix_correlation.pdf")
ggplot(comb_frac, aes(YRI_rfmix, YRI_admixture)) + geom_point() + 
  theme_pubr() + xlim(0.5, 1) + ylim(0.5, 1) + geom_abline(slope = 1) + 
  xlab("YRI fraction (RFMix)") + ylab("YRI fraction (ADMIXTURE)") +
  theme(text = element_text(size=20))
dev.off()

# Plot global ancestry fraction for all individuals using ADMIXTURE estimates
# Histogram
pdf("plots/global_anc_histogram.pdf")
ggplot(comb_frac, aes(YRI_admixture)) + geom_histogram(bins=50, fill="#e48671") + theme_pubr() +
  theme(text = element_text(size=20)) + xlim(0, 1) + 
  xlab("African (YRI) global ancestry fraction")
dev.off()

# Plot global ancestry fraction for all individuals using ADMIXTURE estimates
# STRUCTURE-type plot
pdf("plots/admixture_plot.pdf")
tib = comb_frac %>% transmute(id = row_number(), YRI = YRI_admixture, CEU = CEU_admixture) %>%
  gather("population", "prob", YRI, CEU) %>% group_by(id) %>%
  mutate(YRIprob = prob[which(population == "YRI")]) %>%
  arrange(desc(YRIprob)) %>% ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
ggplot(tib, aes(id, prob, fill = population)) + geom_col() + theme_pubr() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.title.x = element_blank(), axis.line.x = element_blank()) +
  ylab("Global ancestry fraction") + ylim(0, 1) + 
  scale_fill_manual(values=c("#455A8E", "#e48671")) +
  theme(text = element_text(size=20), legend.position = "bottom")
dev.off()

pdf("plots/eur_tract_dist.pdf")
ggplot(eur_tracts, aes(length)) + geom_histogram(bins=100, fill="#455A8E") + 
  theme_pubr(x.text.angle = 90) + xlab("length of European tracts (bp)") + 
  theme(text = element_text(size=20))
dev.off()

eur_tracts %>% summarise(max(length))
