# ALL-Sum (Aggregated L0Learn using Summary Statistics)
Please send any suggestions and feedback on our software!

## Code
- L0LearnSum.cpp: main optimization function for L0Learn on summary data
- allsum.R: full analysis workflow for ALL-Sum

## Data
Please access simulated data and LD references from the Harvard Dataverse at https://dataverse.harvard.edu/dataverse/allsum.

## Simulated Data on Chromosome 21
- tuning.[bed/bim/fam/frq]: tuning genotypes and minor allele frequencies
- validation.[bed/bim/fam]: validation genotypes
- sumstat.txt: GWAS summary statistics
- pheno.txt: continuous phenotype
- ref.map: mapping file to join GWAS, LD, and genotype data
- ref_ld.RDS: block LD reference

## Reference Data with ~1.5 million SNPs from HapMap3 + MEGA chips (Harvard Dataverse)
- 1000G_EUR_hm3_mega_chr[1-22].map / UKB_EUR_hm3_mega_chr[1-22].map : mapping filele
- 1000G_EUR_hm3_mega_chr[1-22]_ld.RDS / UKB_EUR_hm3_mega_chr[1-22]_ld.RDS : block LD reference
Combine per-chromosome LD files with the following script:
```
ld_list = list()
for (chrom in 1:2) {
  ld_list[[chrom]] = readRDS(paste0('1000G_EUR_hm3_mega_chr', chrom, '_ld.RDS'))
  print(sum(sqrt(unlist(lapply(big_ld_list[[chrom]], length)))))
}
ld_list = unlist(ld_list, recursive=F)
length(ld_list) # should be 1701
sum(sqrt(unlist(lapply(ld_list, length)))) # should be 1494152
saveRDS(ld_list, '1000G_EUR_hm3_mega_ld.RDS')
```


# Tutorial
## Download package
Need R and gcc on system to run. Download plink2 from https://www.cog-genomics.org/plink/2.0/. 
```
git clone https://github.com/chen-tony/ALL-sum.git

# necessary R packages
Rscript -e 'install.packages(c('optparse', 'Rcpp','Rcpp','RcppArmadillo', 'dplyr', 'glmnet', 'RISCA'))'
```

## Test ALL-Sum
Change the `--plink2` argument to wherever it is installed.
```
# help
Rscript allsum.R -h 

# run on test data
Rscript allsum.R \
--out test \
--sumstat sumstat.txt \
--ref ref \
--plink2 ~/plink2 \
--tun tuning \
--val validation \
--pheno pheno.txt
```

## Full-genome analysis
Analysis of ~1.5 million SNPs should use around 20GB of memory and 45 minutes of runtime. Note that binary traits will likely take a little longer than continuous traits. Below is an example script 
```
Rscript allsum.R \
--out Trait \
--sumstat Trait_gwas.txt \
--sumstat-name id,chr,pos,ref,alt,stat,n \
--ref Reference/UKB_EUR_hm3_mega \
--plink2 ~/plink2 \
--tun tuning \
--val validation \
--pheno phenotypes.pheno \
--pheno-name FID,IID,Trait \
--cov covariates.cov \
--cov-name FID,IID,age,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10

```

## Creating new reference data
While UKB or 1000 Genomes reference data may work well, it may be preferable to use available tuning (or other reference) data to compute LD. Below are steps to compute new LD blocks to feed into the ALL-Sum workflow. 

### Create .map file
New LD will be computed based on an existing plink file (example name here 'REFERENCE'). In ALL-Sum analysis, replace '--ref Reference/UKB_EUR_hm3_mega' with '--ref REFERENCE'. 

```{r}
library(dplyr)
library(data.table)

blocks = fread('Range/EUR_LDBlocks.txt')

bim = fread('REFERENCE.bim', col.names=c('chr', 'rsid', 'posg', 'pos', 'alt', 'ref'))

full_table = NULL
for (chrom in 1:22) {
  cat(chrom, '..')
  # subset to chromosome
  bim_chr = filter(bim, chr==chrom)
  
  block_chr = filter(blocks, chr == chrom)
  
  block_split = split(block_chr, by=c('chr', 'block'))
  
  # split by blocks
  ix_list = lapply(block_split, FUN=function(X) {
    with(bim_chr, which(pos >= X$start & 
                      pos <= X$stop))
  })
  
  # append to bim
  ix_table = rbindlist(lapply(1:length(ix_list), FUN=function(i) 
    if (length(ix_list[[i]] > 0)) data.frame(bim_chr[ix_list[[i]], ], block=i)))
  
  full_table = rbind(full_table, ix_table)
}
cat('\n')

full_table = full_table %>%
  group_by(chr, block) %>%
  mutate(ix=row_number()) %>% 
  ungroup()

fwrite(full_table, 'REFERENCE.map', sep='\t')
```

### Compute LD blocks using plink
For 1000 Genomes EUR (~500 samples, 1.5 SNPs), this should take about 10 minutes and 100 MB. For larger data such as UKB EUR (~300k samples, 1.5 SNPs), this can take up to 10 hours and 750 MB. This procedure can also be done in parallel by running for separate chromosomes. In practice, a small subset with a few thousand samples should be enough. 

```
# download Ranges for block positions
unzip Range/EUR.zip

mkdir LD

for chrom in {1..22}; do

# number of LD blocks 
n_blocks=$(ls Range/EUR/chr_${chrom}_* | wc -l)
echo $chrom $n_blocks

for ((block=1; block<=$n_blocks; block++)); do

echo -n $block ..

# calculate LD for each block
plink --silent --bfile REFERENCE \
--extract range Range/EUR/chr_${chrom}_block_${block}.range \
--r square \
--allow-no-sex \
--silent \
--write-snplist \
--out LD/chr_${chrom}_block_${block}

done
echo 

done
```
                                
### Check alignment of SNPs and compile LD blocks into list
Saving LD for 1.5 million SNPs should take about 10 GB of memory. 

```{r}
library(dplyr)
library(data.table)

# get number of blocks within map file
map = fread('REFERENCE.map', col.names=c('chr', 'rsid', 'posg', 'pos', 'alt', 'ref', 'block', 'ix'))

map_blocks = map %>% 
  group_by(chr, block) %>% 
  summarize(n=n()) %>%
  ungroup() %>%
  mutate(name = paste0('LD/chr_', chr, '_block_', block))

# match with reference block information
blocks = fread('EUR_LDBlocks.txt')

block_file = blocks %>%
  group_by(chr) %>%
  mutate(block = row_number()) %>%
  ungroup() %>%
  mutate(name = paste0('LD/chr_', chr, '_block_', block)) %>%
  filter(name %in% map_blocks$name) %>%
  pull(name)

# check SNP alignment
snp_list = lapply(block_file, FUN=function(x) fread(paste0(x, '.snplist'), header=F)$V1)

all.equal(unlist(snp_list), map$chr)

# compile LD blocks into list
ld_list = lapply(block_file, FUN=function(x) data.matrix(fread(paste0(x, '.ld'))))

ld_list = Filter(f=function(x) length(x) > 0, ld_list) # remove any potentially empty blocks
length(ld_list) # total number of blocks
sum(unlist(lapply(ld_list, nrow))) # verify correct number of SNPs

saveRDS(ld_list, 'REFERENCE_ld.RDS') 
```
