library(dplyr)
library(data.table)

out = '/n/holystore01/LABS/xlin/Lab/tonychen/L0LearnSum/Blocks/'

# block positions
blocks = fread('EUR_LDBlocks.txt') %>%
  mutate(chr = readr::parse_number(chr))

# split for LD blocks (extend first/last blocks)
blocks_split = blocks %>%
  group_by(chr) %>%
  mutate(block = row_number()) %>%
  mutate(start = ifelse(block == 1, 1, start)) %>%
  mutate(stop = ifelse(block == length(block), Inf, stop)) %>%
  ungroup() %>%
  data.table() %>%
  split(by=c('chr', 'block'))

n_blocks = blocks %>%
  group_by(chr) %>%
  summarize(n_blocks = n()) %>%
  pull(n_blocks)

for (chrom in 1:22) {
  cat(chrom, '..')
  for (block in 1:n_blocks[chrom]) {
    range = blocks_split[[paste0(chrom, '.', block)]]
    write.table(range, paste0(out, 'chr_', chrom, '_block_', block, '.range'),
                row.names=F, col.names=F, quote=F)
  }
}
cat('\n')

# construct LD blocks in plink
file = '/n/holystore01/LABS/xlin/Lab/tonychen/EUR/all_chr_tuning_hm3'
file = '/n/holyscratch01/xlin/tonychen/UKB/Genotypes/ukb_common'

system(paste0('sbatch block.sh ', file))

# compile into list
ld_list = list()
for (chr in 1:22) {
  cat(chrom, '..')
  for (block in 1:n_blocks[chrom]) {
    ld_list = append(ld_list, data.matrix(fread(paste0(out, 'LD/chr_', chrom, 
                                                       '_block_', block, '.ld'))))
  }
}
cat('\n')

ld_list = unlist(ld_list, recursive=F)
ld_list = Filter(f=function(x) length(x) > 0, ld_list)
length(ld_list)
summary(unlist(lapply(ld_list, length)))
sum(unlist(lapply(ld_list, nrow)))

saveRDS(ld_list, paste0(file, '_ld.RDS'))



