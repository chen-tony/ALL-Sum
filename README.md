# ALL-Sum (Aggregated L0Learn using Summary Statistics)

## Data
- [EUR/AFR/ASN]_LDBlocks.txt: LD block information based on Berisa and Pickrell (2016)

- Range/[EUR/AFR/ASN]/: directory of set-range formatted files for computing LD blocks in plink

- Test/: directory of example data on chromosome 21

## Code
- L0LearnSum.cpp: main optimization function for L0Learn on summary data

- ALL_Sum_pipeline.R: full analysis pipeline for ALL-Sum

## Reference Data with ~1.5 million SNPs from HapMap3 + MEGA chips (Dropbox links?)
- 1000G_EUR_hm3_mega[.map/_ld.RDS]: based on 253 European samples in 1000 Genomes Project (Phase 3) 

- UKB_EUR_hm3_mega[.map/_ld.RDS]: based on 20,000 European samples in UK Biobank 

# Tutorial
## Download package
Need R and gcc on system to run. Download plink2 from https://www.cog-genomics.org/plink/2.0/. 
```
git clone https://github.com/chen-tony/ALL-sum.git

# necessary R packages
Rscript -e 'install.packages(c('optparse', 'Rcpp','Rcpp','RcppArmadillo', 'dplyr', 'glmnet', 'RISCA'))'

# download Reference data
wget Reference.zip
```

## Run ALL-Sum
Change the `--plink2` argument to wherever it is installed.
### Test data
```
# download test data
wget Test.zip
unzip Test.zip

# help
Rscript allsum.R -h 

# run on test data
Rscript allsum.R \
--out Test/test \
--sumstat Test/sumstat.txt \
--ref Test/ref \
--plink2 ~/plink2 \
--tun Test/tuning \
--val Test/validation \
--pheno Test/pheno.txt
```

### Full-genome analysis
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
While UKB or 1000 Genomes reference data may work well, it may be preferable to use available tuning (or other reference) data to compute LD. Below are steps to compute new LD blocks to feed into the ALL-Sum pipeline. 

### Create .map file
New LD will be computed based on an existing plink file (example name here 'REFERENCE'). In ALL-Sum analysis, replace '--ref Reference/UKB_EUR_hm3_mega' with '--ref REFERENCE'. 

```{r}
library(dplyr)
library(data.table)

blocks = fread(paste0(out, 'EUR_LDBlocks.txt'))

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
wget Range.zip

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
