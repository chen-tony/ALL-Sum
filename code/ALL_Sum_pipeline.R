library(Rcpp)
library(RcppArmadillo)
library(data.table)
library(dplyr)

invisible(sourceCpp('L0LearnSum.cpp'))

setwd(utils::getSrcDirectory(function(){}))

# out_file: root of output file

# sumstat_file: summary statistics (header: id, chr, pos, ref, alt, stat, n)

# ld_file: LD reference (extensions: .bed, .bim., .fam, _ld.RDS)
# ld_snp_file: SNPs in LD reference

# tun_file: tuning genotype (extensions: .bed, .bim, .fam, .frq)
# val_file: validation genotype (extensions: .bed, .bim, .fam)

# pheno_file: phenotype data
# pheno_name: column name/number in pheno_file (include FID/IID)

# cov_file: covariate data
# cov_name: column name/number in cov_file (include FID/IID)

# options: modify options (parameter grid, partial sorting, ensembling step)

ALL_Sum_pipeline = function(sumstat_file, sumstat_names,
                            ld_file, ld_snp_file,
                            out_file = 'test',
                            tun_file = NULL, val_file = NULL, frq_file = NULL, 
                            pheno_file = NULL, pheno_name = NULL,
                            cov_file = NULL, cov_name = NULL,
                            options=list()) {
  #####################
  ## Data Processing ##
  #####################
  # read in summary statistics
  cat('summary statistics \n')
  sumstat = fread(sumstat_file) 
  
  # match SNPs between sumstats and LD reference
  ld_map = fread(paste0(ld_file, '.bim'),
                 col.names=c('chr', 'id', 'posg', 'pos', 'alt', 'ref'))
  
  combined = ld_map %>% 
    left_join(map, by=c('id'), suffix=c('.gen', '.sum')) %>%
    
    # remove incorrectly specified SNPs
    filter((ref.gen == ref.sum) & (alt.gen == alt.sum) | (ref.gen == alt.sum) & (alt.gen == ref.sum)) 
    
  # check for flipped effects
  flipped = (ref.gen != ref.sum)
  combined$stat[flipped] = -combined$stat[flipped]
  
  cat(nrow(sumstat), 'variants,', 
      nrow(combined), 'kept,', 
      sum(flipped), 'effects flipped \n')
  
  # convert summary statistics
  r = with(combined, stat / sqrt(n - 2 + stat^2))
  
  # hold onto map for plink output
  map0 = combined[,c('id', 'alt')]
  
  rm(sumstat, ld_map); invisible(gc())
  
  # read in LD blocks
  cat('LD blocks \n')
  
  ld_list = readRDS(paste0(ld_file, '_ld.RDS'))
  
  block_sizes = unlist(lapply(ld_list, nrow))
  starts = cumsum(c(1, block_sizes[-length(ld_list)]))
  stops = cumsum(block_sizes) 
  
  # partial sorting
  r_list = lapply(1:length(ld_list), FUN=function(i) r[starts[i]:stops[i]])
  ix_sort = lapply(r_list, FUN=function(x) order(-abs(x)) - 1) # blockwise indices 
  
  rm(block_sizes, r_list, ld_list); invisible(gc())
  
  #####################
  ## Fit L0Learn-Sum ##
  #####################
  # parameter grid
  grid = cbind(s=rep(0, 5),
               lambda0_start = exp(seq(log(1e-5), log(1e-3), length.out=5)),
               lambda1 = rep(0, 5),
               lambda2 = exp(seq(log(1e2), log(1e-1), length.out=5)),
               alpha = seq(0.85, 0.92, length.out=5))
  n_lambda0 = 50
  
  n_par = n_lambda0 * nrow(grid)
  cat(n_par, 'parameters \n')
  
  # use tuning MAF to rescale effect sizes
  if (file.exists(tun_file)) {
    freqs = read_bim(paste0(tuning_file, '.frq'))
    rescales = with(freqs, 1 / sqrt(2*MAF*(1-MAF)))
  } else {
    rescales = 1
  }
  
  # run PRS / return beta on original scale
  beta = matrix(0, length(r), n_par)
  par_out = L0LearnSum_auto(beta, ld_list, r, ix_sort, starts-1, stops-1,
                            grid, n_lambda0, scaling=rescales, active=F, psi=F)
  colnames(par_out) = c('s', 'lambda0', 'lambda1', 'lambda2', 'M', 'L0', 'L1', 'L2', 'conv', 'totit')
  par_out = data.frame(par_out)
  
  # remove SNPs that are zero for all parameters
  nonzero = which(rowSums(beta^2) > 0)
  beta = data.table(map0[nonzero,], beta[nonzero,])
  
  fwrite(beta, paste0(out_file, '_beta.txt'), sep=' ')
  
  ############
  ## Tuning ##
  ############
  if (!is.null(tun_file)) {
    # calculate tuning PRS
    run_tun = system(paste0('plink2 --bfile ', tun_file,
                            ' --score ', out_file, '_beta.txt ',
                            'cols=+scoresums,-maybesid,-dosagesum,-scoreavgs ', 
                            '1 2 header --score-col-nums 3-', n_par+2, 
                            ' --out ', out_file, '_tun'), intern=T)
    
    # read in phenotype / covariate data
    fam = fread(paste0(tun_file, '.fam'),
                col.names=c('FID', 'IID', 'pat', 'mat', 'sex', 'pheno'))
    
    pheno = data.matrix(fread(pheno_file))[,c('FID', 'IID', pheno_name)]
    
    # evaluate tuning PRS
    prs_tun = fread(paste0(out_file, '_tun.sscore')) %>%
      select(FID, IID, contains('SCORE')) 
    
    if (!is.null(cov_file)) {
      cov = data.matrix(fread(cov_file))[,c('FID', 'IID', cov_name)]
      
      tun_data = fam %>%
        select(FID, IID) %>%
        left_join(pheno, by=c('FID', 'IID')) %>%
        left_join(prs_tun, by=c('FID', 'IID')) %>%
        left_join(cov, by=c('FID', 'IID')) %>%
        na.omit()
      
      residuals = residuals(lm(pheno ~ .))
      
      rm(fam, prs_tun, cov); invisible(gc())
      
    } else {
      tun_data = fam %>%
        select(FID, IID) %>%
        left_join(pheno, by=c('FID', 'IID')) %>%
        left_join(prs_tun, by=c('FID', 'IID'))
      
      rm(fam, prs_tun); invisible(gc())
    }
    
    tuning$r2_tun = apply(prs_tun, 2, FUN=function(x) ifelse(sd(x) > 0, cor(x, tun.pheno)^2, 0))
  }
}


print('genotype')

# SNP information
freqs = read_bim(paste0(tuning_file, '.frq'))
rescales = with(freqs, 1 / sqrt(2*MAF*(1-MAF)))

map = read_bim(paste0(tuning_file, '.bim'))

##########################
### summary statistics ###
##########################
print('summary statistics')

sumstat = fread(sumstat) 

# match summary statistics to genotype data
combined = map %>% 
  mutate(ix = row_number()) %>%
  
  # join bim with sumstats
  left_join(map, by=c('id'), suffix=c('.gen', '.sum')) %>%
  
  # remove incorrectly specified SNPs
  filter((ref.gen == ref.sum) & (alt.gen == alt.sum) | (ref.gen == alt.sum) & (alt.gen == ref.sum)) %>%
  
  # check for flipped effects
  mutate(flipped = (ref.gen == ref.sum) & (alt.gen == alt.sum))

flipped = ((combined$REF == combined$alt)&(combined$ALT == combined$ref))
combined$BETA[flipped] = -combined$BETA[flipped]
combined$T_STAT[flipped] = -combined$T_STAT[flipped]

# convert summary statistics
r = with(combined, t2cor(T_STAT, OBS_CT))

# save for plink output
map0 = combined[,c('ID', 'alt')]

rm(snplist, freqs, rescales, map,
   pheno, sumstat, combined, flipped); gc() # save memory

#################
### LD blocks ###
#################
print('LD blocks')

ld_list = readRDS(ld_file)

block_sizes = unlist(lapply(ld_list, nrow))
starts = cumsum(c(1, block_sizes[-length(ld_list)])) - 1
stops = cumsum(block_sizes) - 1

# partial sorting
r_list = lapply(1:length(ld_list), FUN=function(i) r[starts[i]:stops[i]])
ix_sort = lapply(r_list, FUN=function(x) order(-abs(x)) - 1) # blockwise indices 

rm(block_sizes, r_list); gc()

##################
### L0LearnSum ###
##################
print('L0LearnSum')

# tuning parameters
grid = cbind(s=rep(0, 5),
             lambda0_start = exp(seq(log(1e-5), log(1e-3), length.out=5)),
             lambda1 = rep(0, 5),
             lambda2 = exp(seq(log(1e2), log(1e-1), length.out=5)),
             alpha = seq(0.92, 0.96, length.out=5) )

n_lambda0 = 100

n_par = n_lambda0 * nrow(grid)

# run PRS / return beta on original scale
beta = matrix(0, length(r), n_par)
tuning = L0LearnSum_auto(beta, ld_list, r, ix_sort, starts, stops,
                         grid, n_lambda0, scaling=rescales, active=F, psi=F)

colnames(tuning) = c('s', 'lambda0', 'lambda1', 'lambda2', 'M', 'L0', 'L1', 'L2', 'conv', 'totit')
tuning = data.frame(tuning)

# save tuning results
saveRDS(tuning, paste0(output_directory, '/pars_', output_name, '.RDS'))

# format beta for plink
beta = data.table(map0, beta)
fwrite(beta, paste0(output_directory, '/beta_', output_name, '.csv'))

#####################
### calculate PRS ###
#####################
print('PRS')

# use plink
system(paste0('~/software2/plink2 --bfile ', tuning_file, ' --score ', 
              output_directory, '/beta_', output_name, '.csv',
              '1 2 header --score-col-nums 3-502 --out ',
              output_directory, '/prs_tun_', output_name))

system(paste0('~/software2/plink2 --bfile ', validation_file, ' --score ', 
              output_directory, '/beta_', output_name, '.csv',
              '1 2 header --score-col-nums 3-502 --out ',
              output_directory, '/prs_val_', output_name))

##############
### tuning ###
##############
print('tuning')

# read in phenotype
pheno = data.matrix(fread(paste0(gen,'phenotypes_hm3_rho_', rho, '_', GA, '.phen')))

tun.pheno = pheno[100001:110000, (rep+2)]; # colnames(tun.pheno) = c('FID', 'IID', 'PHENO')
val.pheno = pheno[110001:120000, (rep+2)]; # colnames(val.pheno) = c('FID', 'IID', 'PHENO')

# read in PRS
prs_tun = fread(paste0('Results/prs_tun_rho_', rho, '_size_', size, '_rep_', rep, '_GA_', GA, '.sscore'))
prs_val = fread(paste0('Results/prs_val_rho_', rho, '_size_', size, '_rep_', rep, '_GA_', GA, '.sscore'))

# calculate R2
tuning$r2_tun = apply(prs_tun, 2, FUN=function(x) ifelse(sd(x) > 0, cor(x, tun.pheno)^2, 0))
tuning$r2_val = apply(prs_val, 2, FUN=function(x) ifelse(sd(x) > 0, cor(x, val.pheno)^2, 0))

# re-save tuning results w/ R2
saveRDS(tuning, paste0('Results/results_rho_', rho, '_size_', size, '_rep_', rep, '_GA_', GA, '.RDS'))
system(paste0(paste0('rm Results/pars_rho_', rho, '_size_', size, '_rep_', rep, '_GA_', GA, '.RDS')))

##################
### ensembling ###
##################
print('ensembling')

# remove bad PRS
good_ix = which(apply(prs_tun, 2, sd) > 0) # remove any 0 PRS
prs_tun = prs_tun[,good_ix]
prs_val = prs_val[,good_ix]
tuning = tuning[good_ix,]

# sort by R2 and remove correlated PRS
cor_prs = cor(prs_tun)

prs_tun_order = order(-tuning$r2_tun)
ix_keep = prs_tun_order[1]
for (i in 2:length(prs_tun_order)) {
  if (max(abs(cor_prs[ix_keep, prs_tun_order[i]])) < 0.98) 
    ix_keep = c(ix_keep, prs_tun_order[i])
}
ix_keep = sort(ix_keep)

# run lasso on PRS
fit_glmnet = cv.glmnet(prs_tun[,ix_keep], tun.pheno, nfolds=3)

best_prs_tun = c(predict(fit_glmnet, prs_tun0[,ix_keep], s='lambda.min'))
best_prs_val = c(predict(fit_glmnet, prs_val0[,ix_keep], s='lambda.min'))
ix_prs_select = which(coef(fit_glmnet, s='lambda.min')[-1] != 0)

saveRDS(best_prs_val, paste0('Results/super_prs_val_', type, '_', pen,
                             '_rho_', rho, '_size_', size,
                             '_rep_', rep, '_GA_', GA, '.RDS'))

beta_keep = beta[,c(1, 2, (good_ix[ix_prs_select] + 2))]
nonzero = which(rowSums(beta_keep != 0) > 0)
n_select = length(nonzero)

best_results = data.frame(rho=rho, size=size, rep=rep, GA=GA,
                          n_select=n_select,
                          n_prs=length(ix_prs_select),
                          r2_tun=cor(best_prs_tun, tun.pheno)^2,
                          r2_val=cor(best_prs_val, val.pheno)^2)

saveRDS(best_results, paste0('Results/super_results_', type, '_', pen,
                             '_rho_', rho, '_size_', size,
                             '_rep_', rep, '_GA_', GA, '.RDS'))



