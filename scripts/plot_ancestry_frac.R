library(tidyverse)

bed_files = list.files("data/bed/")

# Initialize large empty vector to store tracts - length is unknown, but we want to avoid spending time/space extending it in the for loop
init_length = 4 * length(bed_files)
eur_tract_lengths = rep(0, init_length)

cur_index = 1
buffer_full = FALSE # Boolean variable to keep track of whether we ran out of room in the inital vector and need to extend it
for (bf in bed_files) {
  bf_tib = read_tsv(paste("data/bed/", bf, sep=''), col_names=c('chrom', 'gpos_start', 'gpos_stop', 'anc', 'pos_start', 'pos_stop'))
  tracts = bf_tib %>% filter(anc == 'CEU') %>% mutate(size = pos_stop - pos_start) %>% pull(size)
  if (length(tracts) > 0) {
    if (buffer_full) {
      eur_tract_lengths = c(eur_tract_lengths, tracts)
    } else if (cur_index + length(tracts) - 1 > init_length) {
      eur_tract_lengths = eur_tract_lengths[eur_tract_lengths != 0]
      eur_tract_lengths = c(eur_tract_lengths, tracts)
      buffer_full = TRUE
    } else {
      eur_tract_lengths[cur_index:cur_index + length(tracts) - 1] = tracts 
      cur_index = cur_index + length(tracts)
    }
  }
}

eur_tract_lengths = eur_tract_lengths[eur_tract_lengths != 0]

tib = tibble(eur_tract_lengths)

pdf("plots/eur_tract_lengths.hist.pdf")
ggplot(tib, aes(eur_tract_lengths)) + geom_histogram()
dev.off()
