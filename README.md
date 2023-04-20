# ALL-Sum (Aggregated L0Learn using Summary Statistics)

## Data
- [EUR/AFR/ASN]_LDBlocks.txt: LD block information based on Berisa and Pickrell (2016)

- Range/[EUR/AFR/ASN]/: set-range formatted files for computing LD blocks in plink

- Test/: example data on chromosome 21 for testing

## Code
- L0LearnSum.cpp: main optimization function

- construct_blocks.html / block.sh: LD block computation

- ALL_Sum_pipeline.R: full analysis pipeline for ALL-sum

## Reference Data with ~1.5 million SNPs from HapMap3 + MEGA chips (Dropbox links?)
- 1000G_EUR_hm3_mega[.map/_ld.RDS]: based on 253 European samples in 1000 Genomes Project (Phase 3) 
- UKB_EUR_hm3_mega[.map/_ld.RDS]: based on 20,000 European samples in UK Biobank 

## Test Data
```
unzip Test.zip
```

## Range files for LD block construction
```
unzip Range.zip
```

# Tutorial
## Download GitHub
```
git clone https://github.com/chen-tony/ALL-sum.git
```

## Load dependencies
Load the following R libraries (must have gcc and R on system).

```
Rscript -e 'install.packages(c('optparse', 'Rcpp','Rcpp','RcppArmadillo', 'dplyr', 'glmnet', 'RISCA'))'
```

Download plink2: https://www.cog-genomics.org/plink/2.0/

## Run ALL-Sum
Change the `--plink2` argument to wherever it is installed 

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



## Running ALL-Sum from scratch
First, create a ".map" file - appending LD block information to the relevant ".bim" file. For this example, we will be constructing EUR-based LD using a plink file called "ref". 

```{r}
blocks = fread(paste0(out, 'EUR_LDBlocks.txt'))

bim = fread('ref.bim', col.names=c('chr', 'rsid', 'posg', 'pos', 'alt', 'ref'))

full_table = NULL
for (chrom in 1:22) {
  # subset to chromosome
  bim_chr = filter(bim, chr==chrom)
  
  block_chr = filter(blocks, chr == chrom)
  
  block_split = split(block_chr, by=c('chr', 'block'))
  
  # split by blocks
  ix_list = lapply(block_split, FUN=function(X) {
    with(bim_chr, which(pos >= X$start & 
                      pos <= X$stop))
  })
  
  ix_list = Filter(f=function(x) length(x) > 0, ix_list)
  
  # append to bim
  ix_table = rbindlist(lapply(1:length(ix_list), FUN=function(i) 
    data.frame(bim_chr[ix_list[[i]], ], block=i)))
  
  full_table = rbind(full_table, ix_table)
}

fwrite(full_table, 'ref.map')

```

Then, we can use plink to compute LD blocks based on the positions in the Range directory. It may be fastest to separate this by chromosome. 
```

mkdir LD

for chrom in {1..22}; do

# number of LD blocks 
n_blocks=$(ls Range/EUR/chr_${chrom}_* | wc -l)
echo $chrom $n_blocks

for ((block=1; block<=$n_blocks; block++)); do

echo -n $block ..

# calculate LD for each block
plink --silent --bfile ref \
--extract range ${block_dir}/Range/chr_${chrom}_block_${block}.range \
--r square \
--allow-no-sex \
--silent \
--write-snplist \
--out LD/chr_${chrom}_block_${block}

done
echo 

done

```

Check that everything has been aligned properly
```{r}
map = fread('ref.map', ol.names=c('chr', 'rsid', 'posg', 'pos', 'alt', 'ref', 'block'))

map_blocks = map %>% 
  group_by(chr, block) %>% 
  summarize(n=n()) %>%
  ungroup() %>%
  mutate(name = paste0('LD/chr_', chr, '_block_', block))

blocks = fread('EUR_LDBlocks.txt')

block_file = blocks %>%
  group_by(chr) %>%
  mutate(block = row_number()) %>%
  ungroup() %>%
  mutate(name = paste0('LD/chr_', chr, '_block_', block)) %>%
  filter(name %in% map_blocks$name) %>%
  pull(name)

snp_list = lapply(block_file, FUN=function(x) fread(paste0(x, '.snplist'), header=F)$V1)

all.equal(unlist(snp_list), map$chr)

```

Finally, compile the blocks into an R list.
```{r}
ld_list = lapply(block_file, FUN=function(x) data.matrix(fread(paste0(x, '.ld'))))

ld_list = Filter(f=function(x) length(x) > 0, ld_list) # remove any potentially empty blocks
length(ld_list) # total number of blocks
sum(unlist(lapply(ld_list, nrow))) # check total number of SNPs

saveRDS(ld_list, 'ref_ld.RDS') # save

```

Analysis of ~1.5 million SNPs should use around 20GB of memory and 45 minutes of runtime. Note that binary traits will likely take a little longer than continuous traits. 
```
# download LD reference data
wget Reference.zip
unzip Reference.zip

# run on real data (example syntax)
Rscript allsum.R \
--out trait \
--sumstat trait_gwas.txt \
--sumstat-name rsid,chr,pos,a0,a1,stat,n_eff \
--ref ref \
--plink2 ~/plink2 \
--tun tuning \
--val validation \
--pheno phenotypes.pheno \
--pheno-name FID,IID,trait \
--cov covariates.cov \
--cov-name FID,IID,age,sex,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10




```
