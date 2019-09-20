library(tidyverse)

anc = read_tsv("data/chr1/chr1.combined_global_anc_frac.txt")
anc = anc %>% mutate(chr='1')
# pdf("plots/chr1.correlate_rfmix_admixture.pdf")
# print(ggplot(anc, aes(CEU_rfmix, CEU_admixture)) + geom_point())
# dev.off()

for (i in 2:22) {
  tmp = read_tsv(paste("data/chr", i, "/chr", i, ".combined_global_anc_frac.txt", sep=''))
  tmp = tmp %>% mutate(chr=toString(i))
  anc = anc %>% bind_rows(tmp)
  # pdf(paste("plots/chr", i, ".correlate_rfmix_admixture.pdf", sep=''))
  # print(ggplot(tmp, aes(CEU_rfmix, CEU_admixture)) + geom_point())
  # dev.off()
}

# pdf("plots/all.correlate_rfmix_admixture.pdf")
# ggplot(anc, aes(CEU_rfmix, CEU_admixture)) + geom_point()
# dev.off()

chr_len = read_tsv("data/maps/chrom_lengths.tsv", col_names=c('chrom', 'length'))

tib = anc %>% left_join(chr_len, by=c('chr'='chrom')) %>% mutate(YRI_admix_bp = YRI_admixture * length, CEU_admix_bp = CEU_admixture * length, CEU_rfmix_bp = CEU_rfmix * length, YRI_rfmix_bp = YRI_rfmix * length) %>% group_by(X1) %>% summarise(YRI_admix_sum = sum(YRI_admix_bp), CEU_admix_sum = sum(CEU_admix_bp), YRI_rfmix_sum = sum(YRI_rfmix_bp), CEU_rfmix_sum = sum(CEU_rfmix_bp)) %>% mutate(tot_admix = YRI_admix_sum + CEU_admix_sum, tot_rfmix = YRI_rfmix_sum + CEU_rfmix_sum) %>% mutate(YRI_admix = YRI_admix_sum/tot_admix, CEU_admix = CEU_admix_sum/tot_admix, YRI_rfmix = YRI_rfmix_sum/tot_rfmix, CEU_rfmix = CEU_rfmix_sum/tot_rfmix) %>% select(YRI_admix, CEU_admix, YRI_rfmix, CEU_rfmix)
# pdf("plots/all.correlate_rfmix_admixture.pdf")
# ggplot(tib, aes(CEU_rfmix, CEU_admix)) + geom_point()
# dev.off()

write_delim(tib, "data/all.combined_global_anc_frac.txt", delim='\t')

pdf("plots/rfmix.eur_prop.hist.pdf")
ggplot(tib, aes(CEU_rfmix)) + geom_histogram()
dev.off()

pdf("plots/admixture.eur_prop.hist.pdf")
ggplot(tib, aes(CEU_admix)) + geom_histogram()
dev.off()