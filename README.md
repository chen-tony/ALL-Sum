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



# Tutorial
## Download GitHub
```
git something
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

For full-genome analysis (~1.5 million SNPs), allow around 20GB of memory and 45 minutes of runtime. Note that analysis of binary traits will likely take a little longer than continuous traits. 
