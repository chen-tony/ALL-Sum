# ALL-Sum
Aggregated L0Learn using Summary Statistics



## Data
- [EUR/AFR/ASN]_LDBlocks.txt: LD block information based on Berisa and Pickrell (2016)

- Range/[EUR/AFR/ASN]/: set-range formatted files for computing LD blocks in plink

- Test/: example data on chromosome 21 for testing

## Reference Data with 1.5 SNPs from HapMap3 + MEGA chips (Dropbox links?)
- 1000G_EUR_hm3_mega[.map/_ld.RDS]: based on 253 European samples in 1000 Genomes Project (Phase 3) 
- UKB_EUR_hm3_mega[.map/_ld.RDS]: based on 20,000 European samples in UK Biobank 

## Code
- L0LearnSum.cpp: main optimization function

- construct_blocks.html / block.sh: LD block computation

- ALL_Sum_pipeline.R: full analysis pipeline for ALL-sum
