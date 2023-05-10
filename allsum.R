############
## Inputs ##
############
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c('-o', '--out'), type='character', 
              default='allsum', 
              help='root of output file', 
              metavar='character'),
  
  make_option(c('-s', '--sumstat'), type='character', 
              default=NULL,
              help='sumstat file [assumed header: id, chr, pos, ref, alt, stat, n]', 
              metavar='character'),
  make_option('--sumstat-name', type='character', 
              default=NULL,
              help='comma-separated column names [example: id,chr,pos,ref,alt,stat,n]', 
              metavar='character'),
  make_option('--sumstat-col', type='character', 
              default=NULL,
              help='comma-separated column numbers [example: 1,2,3,4,5,6,7]', 
              metavar='character'),
  
  make_option(c('-r', '--ref'), type='character', 
              default='Reference/1000G_EUR_hm3_mega', 
              help='LD reference (default: 1000G-EUR) [extensions: .map, _ld.RDS]', 
              metavar='character'),
  make_option('--plink2', type='character', 
              default=NULL, 
              help='path to plink2', 
              metavar='character'),
  make_option(c('-t', '--tun'), type='character', 
              default=NULL, 
              help='plink-format tuning genotype [extensions: .bed, .bim, .fam, (optional: .frq)]', 
              metavar='character'),
  make_option(c('-v', '--val'), type='character', 
              default=NULL, 
              help='plink-format validation genotype [extensions: .bed, .bim, .fam]', 
              metavar='character'),
  
  make_option(c('-p', '--pheno'), type='character', 
              default=NULL, 
              help='phenotype file [assumed header: FID IID phenotype]', 
              metavar='character'),
  make_option('--pheno-name', type='character', 
              default=NULL, 
              help='comma-separated column names [example: FID,IID,phenotype]', 
              metavar='character'),
  make_option('--pheno-col', type='character', 
              default=NULL, 
              help='comma-separated column numbers [example: 1,2,3]', 
              metavar='character'),
  
  make_option(c('-c', '--cov'), type='character', 
              default=NULL, 
              help='covariate data [assumed header: FID IID X1 -- Xd]', 
              metavar='character'),
  make_option('--cov-name', type='character', 
              default=NULL, 
              help='comma-separated column names [example: FID,IID,X1,...,Xd]', 
              metavar='character'),
  make_option('--cov-col', type='character', 
              default=NULL, 
              help='comma-separated column numbers [example: 1,2,3,...,d+2]', 
              metavar='character'),
  
  make_option('--active', type='integer', 
              default=0, 
              help='use active set for optimization [include number of repeated active sets]', 
              metavar='integer'),
  
  make_option('--match-by-pos', type='character', 
              default=NULL, 
              action='store_true',
              help='match by chr/pos instead of RSID + use RSID from map or sum', 
              metavar=NULL),
  
  make_option(c('-g', '--grid'), type='character', 
              default=NULL, 
              help='file for parameter grid [header: lambda0_start lambda1 lambda2 alpha]', 
              metavar='character'),
  make_option('--nlambda0', type='integer',
              default=50,
              help='number of lambda0 values',
              metavar='integer'),
  make_option('--seed', type='integer', 
              default=123, 
              help='random seed', 
              metavar='integer')
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser, convert_hyphens_to_underscores = T)

if (is.null(opt$match_by_pos)) opt$match_by_pos = 'map'

## To implement later ##
# --ens: other ensembling options (glmnet weighting)
# --boot: bootstrap confidence interval to prediction

##########################
### Check input errors ###
##########################
# missing required inputs
if (is.null(opt$sumstat)) {
  cat('ERROR: please provide summary statistics')
  quit('no')
}
if (is.null(opt$plink2)) {
  cat('ERROR: please provide path for plink2')
  quit('no')
}
if (is.null(opt$sumstat)) {
  cat('ERROR: please provide summary statistics')
  quit('no')
}

# incorrect header inputs
if (!is.null(opt$sumstat_name) & !is.null(opt$sumstat_col)) {
  cat('ERROR: please use only one of --sumstat-name or --sumstat-col')
  quit('no')
}
if (!is.null(opt$pheno_name) & !is.null(opt$pheno_col)) {
  cat('ERROR: please use only one of --pheno-name or --pheno-col')
  
  if (is.null(opt$pheno)) {
    cat('ERROR: please provide phenotype data')
  }
  
  quit('no')
}
if (!is.null(opt$cov_name) & !is.null(opt$cov_col)) {
  cat('ERROR: please use only one of --cov-name or --cov-col')
  
  if (is.null(opt$pheno)) {
    cat('ERROR: please provide covariate data')
  }
  quit('no')
}

{
  cat('Active arguments: \n')
  arg = parse_args(opt_parser, convert_hyphens_to_underscores = F)
  arg$help = NULL
  for (i in names(arg)) {
    cat(paste0('  --', i), arg[[i]], '\n')
  }
  
  invisible(rm(arg)); invisible(gc())
}

################
### Packages ###
################
cat('loading packages...')
suppressPackageStartupMessages({
  library(Rcpp)
  library(RcppArmadillo)
  library(data.table)
  library(dplyr)
  library(glmnet)
})

invisible(sourceCpp('L0LearnSum.cpp'))

set.seed(opt$seed)

cat('\n')
  
#####################
## Data Processing ##
#####################
# read in summary statistics
cat('GWAS: ')
sumstat = fread(opt$sumstat, showProgress=F)

if (!is.null(opt$sumstat_name)) {
  opt$sumstat_name = strsplit(opt$sumstat_name, ',')[[1]]
  
  sumstat = sumstat %>%
    select(any_of(opt$sumstat_name)) %>%
    data.frame()
  
} else if (!is.null(opt$sumstat_col)) {
  opt$sumstat_col = as.numeric(strsplit(opt$sumstat_col, ',')[[1]])
  
  sumstat = data.frame(sumstat)[,opt$sumstat_col]
}
names(sumstat) = c('id', 'chr', 'pos', 'ref', 'alt', 'stat', 'n')

# match SNPs between sumstats and LD reference
map = fread(paste0(opt$ref, '.map'),
            col.names=c('chr', 'id', 'posg', 'pos', 'alt', 'ref', 'block', 'ix'), 
            showProgress=F)

if (is.null(opt$match_by_pos)) {
  # using RSID
  combined = map %>% 
    left_join(sumstat, by=c('id'), suffix=c('.map', '.sum')) %>%
    rename(chr=chr.map) %>%
    
    # remove incorrectly specified SNPs
    filter((ref.map == ref.sum) & (alt.map == alt.sum) | 
             (ref.map == alt.sum) & (alt.map == ref.sum)) 
  
  map0 = combined %>% select(id, alt=alt.map)
} else {
  # using chr/pos
  combined = map %>% 
    left_join(sumstat, by=c('chr', 'pos'), suffix=c('.map', '.sum')) %>%
    
    # remove incorrectly specified SNPs
    filter((ref.map == ref.sum) & (alt.map == alt.sum) | 
             (ref.map == alt.sum) & (alt.map == ref.sum))
  
  map0 = combined %>% select(id=paste0('id.', opt$match_by_pos), alt=alt.map)
}
  
# flip sumstat effects as necessary (line up with LD reference)
flipped = with(combined, ref.map != ref.sum)
combined$stat[flipped] = -combined$stat[flipped]

cat(nrow(sumstat), 'variants,', 
    nrow(combined), 'matched,',
    sum(flipped), 'effects flipped \n')

# convert summary statistics
r = with(combined, stat / sqrt(n - 2 + stat^2))

invisible(rm(sumstat)); invisible(gc())

# read in LD blocks
cat('LD: ')

ld_list = readRDS(paste0(opt$ref, '_ld.RDS'))

cat('matching to sumstats, ')
map_match = map %>%
  filter(ix == 1) %>%
  select(chr, block)

ix_match_list = lapply(1:nrow(map_match), FUN=function(i) {
  combined %>% filter(chr==map_match$chr[i], 
                      block==map_match$block[i]) %>% pull(ix)
})

for (i in 1:length(ld_list)) {
  ix_match = ix_match_list[[i]]
  if (length(ix_match) == 0 | is.null(ix_match)) {
    ld_list[[i]] = integer(0)
  } else {
    ld_list[[i]] = as.matrix(ld_list[[i]][ix_match,ix_match])
  }
}

ld_list = Filter(f=function(x) length(x) > 0, ld_list)
block_sizes = unlist(lapply(ld_list, nrow))
starts = cumsum(c(1, block_sizes[-length(ld_list)]))
stops = cumsum(block_sizes) 

cat(length(ld_list), 'blocks \n')

# partial sorting (blockwise)
r_list = lapply(1:length(ld_list), FUN=function(i) r[starts[i]:stops[i]])
ix_sort = lapply(r_list, FUN=function(x) order(-abs(x)) - 1)

invisible(rm(block_sizes, r_list, combined, map, 
             ix_match, ix_match_list, map_match, flipped))
invisible(gc())

#####################
## Fit L0Learn-Sum ##
#####################

# parameter grid
if (!is.null(opt$grid)) {
  grid = data.frame(fread(opt$grid))
} else {
  grid = cbind(lambda0_start = exp(seq(log(1e-5), log(1e-3), length.out=5)),
               lambda1 = rep(0, 5),
               lambda2 = exp(seq(log(1e2), log(1e-1), length.out=5)),
               alpha = seq(0.85, 0.92, length.out=5))
}
n_lambda0 = opt$nlambda0

n_par = n_lambda0 * nrow(grid)
cat(paste0('L0Learn-Sum (', n_par, ' parameters): '))

# use tuning MAF to rescale effect sizes
if (!is.null(opt$tun) & file.exists(paste0(opt$tun, '.frq'))) {
  freqs = fread(paste0(opt$tun, '.frq'), 
                col.names=c('chr', 'id', 'alt', 'ref', 'maf', 'n'),
                showProgress=F)
  
  maf = map0 %>% 
    left_join(freqs, by='id') %>% 
    pull(maf)
  
  # convert MAF to rescaling factor
  rescales = 1 / sqrt(2*maf*(1-maf))
  
  # set any bad rescaling factors to 0
  rescales[is.na(rescales)] = 0
  rescales[is.infinite(rescales)] = 0
  rescales[is.nan(rescales)] = 0
  
  invisible(rm(freqs, maf)); invisible(gc())
} else {
  rescales = 1
}

# run PRS / return beta on original scale
beta = matrix(0, length(r), n_par)

if (opt$active > 0) {
  # active set updates
  par_out = L0LearnSum_active_auto(beta, ld_list, r, ix_sort, starts-1, stops-1,
                                   grid, n_lambda0, scaling=rescales, max_active=opt$active)
} else {
  # regular coordinate descent
  par_out = L0LearnSum_auto(beta, ld_list, r, ix_sort, starts-1, stops-1,
                            grid, n_lambda0, scaling=rescales)
}

colnames(par_out) = c('lambda0', 'lambda1', 'lambda2', 'M', 'L0', 'L1', 'L2', 'conv', 'totit')
par_out = data.frame(par_out)
beta[is.na(beta)] = 0
beta[is.nan(beta)] = 0
beta[is.infinite(beta)] = 0

fwrite(par_out, paste0(opt$out, '_pars.txt'), sep='\t')

# remove SNPs that are zero for all parameters
nonzero = which(rowSums(beta^2) > 0)
beta = data.table(map0[nonzero,], beta[nonzero,])

fwrite(beta, paste0(opt$out, '_beta.txt'), sep=' ')
cat(paste0('   Effect estimates (', length(nonzero), 'x', n_par, ') written to ', opt$out, '_beta.txt \n'))

invisible(rm(r, ld_list, ix_sort, i,
             grid, n_lambda0, 
             beta, rescales, starts, stops, nonzero))
invisible(gc())

############
## Tuning ##
############
if (!is.null(opt$tun)) {
  #################
  ## L0Learn-Sum ##
  #################
  cat('Tuning: ')
  
  # calculate tuning PRS
  cat('computing PRS, ')
  run_tun = system(paste0(opt$plink2, ' --bfile ', opt$tun,
                          ' --score ', opt$out, '_beta.txt ',
                          'cols=+scoresums,-maybesid,-dosagesum,-scoreavgs ', 
                          '1 2 header --score-col-nums 3-', n_par+2, 
                          ' --out ', opt$out, '_tun'), intern=T)
  system(paste0('rm ', opt$out, '_tun.log'))
  invisible(rm(run_tun)); invisible(gc())
  
  # only look at results that converged and are nonzero
  par_out$pred_tun = NA
  par_out$rev_beta = NA
  
  ix_conv = which(par_out$conv == 1 & par_out$L0 > 0)
  
  # read in phenotype / covariate data
  cat('reading data, ')
  fam = fread(paste0(opt$tun, '.fam'),
              col.names=c('FID', 'IID', 'pat', 'mat', 'sex', 'pheno'),
              showProgress=F) %>%
    data.frame()
  
  pheno = data.frame(fread(opt$pheno, showProgress=F))
  if (!is.null(opt$pheno_name)) {
    opt$pheno_name = strsplit(opt$pheno_name, ',')[[1]]
    
    pheno = pheno[,opt$pheno_name]
  } else if (!is.null(opt$pheno_col)) {
    opt$pheno_col = as.numeric(strsplit(opt$pheno_col, ',')[[1]])
    
    pheno = pheno[,opt$pheno_col]
  }
  names(pheno) = c('FID', 'IID', 'pheno')
  
  # check for binary outcome
  binary_outcome = (length(table(pheno$pheno)) == 2) 
  if (binary_outcome) suppressPackageStartupMessages(library(RISCA))
  
  # set up covariates
  prs_tun = fread(paste0(opt$out, '_tun.sscore'), showProgress=F) %>%
    rename(FID='#FID') %>%
    select(FID, IID, contains('SCORE')) %>%
    data.frame()
  
  # join PRS, phenotype, covariate
  if (!is.null(opt$cov)) {
    
    cov = data.frame(fread(opt$cov, showProgress=F))
    if (!is.null(opt$cov_name)) {
      opt$cov_name = strsplit(opt$cov_name, split=',')[[1]]
      
      cov = cov[,opt$cov_name]
      names(cov)[1:2] = c('FID', 'IID')
      cov_formula = paste(names(cov)[-c(1:2)], collapse='+')
      
    } else if (!is.null(opt$cov_col)) {
      opt$cov_col = as.numeric(strsplit(opt$cov_col, split=',')[[1]])
      
      cov = cov[,opt$cov_col]
      names(cov)[1:2] = c('FID', 'IID')
      cov_formula = paste(names(cov)[-c(1:2)], collapse='+')
    } else {
      names(cov)[1:2] = c('FID', 'IID')
      cov_formula = paste(names(cov)[-c(1:2)], collapse='+')
    }
    
    tun_data = fam %>%
      select(FID, IID) %>%
      left_join(pheno, by=c('FID', 'IID')) %>%
      left_join(prs_tun, by=c('FID', 'IID')) %>%
      left_join(cov, by=c('FID', 'IID')) %>%
      na.omit()
    
  } else {
    cov_formula = '1'
    
    tun_data = fam %>%
      select(FID, IID) %>%
      left_join(pheno, by=c('FID', 'IID')) %>%
      left_join(prs_tun, by=c('FID', 'IID')) %>%
      na.omit()
  }
  invisible(rm(fam, prs_tun)); invisible(gc())
  
  # evaluate fit in tuning data
  cat('evaluating, ')
  if (binary_outcome) {
    # conditional AUC
    par_out$pred_tun[ix_conv] = sapply(ix_conv, FUN=function(i)
      roc.binary(status='pheno',
                 variable=paste0('SCORE', i, '_SUM'),
                 confounders=paste0('~',cov_formula),
                 data=tun_data)$auc)
    
    # need to flip beta if AUC < 0.5
    par_out$rev_beta = NA
    par_out$rev_beta = (par_out$pred_tun < 0.5)
    
    # final AUC
    par_out$pred_tun = pmax(par_out$pred_tun, 1-par_out$pred_tun)
    
  } else {
    # adjust for covariates
    res_tun = residuals(lm(paste0('pheno ~', cov_formula), data=tun_data))
    
    # temporary correlation
    temp_cor = sapply(ix_conv, FUN=function(i) 
      cor(tun_data[,paste0('SCORE', i, '_SUM')], res_tun))
    
    # need to flip beta if temporary correlation < 0.5
    par_out$rev_beta[ix_conv] = (temp_cor < 0)
    
    # final R2 = correlation^2
    par_out$pred_tun[ix_conv] = temp_cor^2
  }
  
  # choose best beta
  best_tuning_ix = which.max(par_out$pred_tun)
  beta_col = c('id', 'alt', paste0('V', best_tuning_ix))
  
  beta = fread(paste0(opt$out, '_beta.txt'),
               select=beta_col, showProgress=F) %>% 
    data.frame()
  if (par_out[best_tuning_ix, 'rev_beta'] == T) beta[,3] = -beta[,3]
  
  best_beta = beta[which(beta[,3] != 0), ]
  cat(nrow(best_beta), 'selected SNPs, ')
  
  fwrite(best_beta, paste0(opt$out, '_best_beta.txt'), sep=' ')
  
  cat('saving results \n')
  fwrite(par_out, paste0(opt$out, '_pars.txt'), sep='\t')
  
  cat(paste0('   Tuning summary written to ', opt$out, '_pars.txt \n'))
  
  invisible(rm(beta, best_beta, beta_col)); invisible(gc())
  
  #############
  ## ALL-Sum ##
  #############
  cat('Ensembling: ')
  # run CV lasso
  if (binary_outcome) {
    fit_glmnet = cv.glmnet(data.matrix(tun_data[,paste0('SCORE', ix_conv, '_SUM')]), 
                           tun_data$pheno, nfolds=3, family='binomial')
    
  } else {
    fit_glmnet = cv.glmnet(data.matrix(tun_data[,paste0('SCORE', ix_conv, '_SUM')]), 
                           res_tun, nfolds=3, family='gaussian')
    
  }
  
  # combine beta
  best_weights = coef(fit_glmnet, s='lambda.min')[-1]
  ix_prs_select = which(best_weights != 0)
  coef_glmnet = best_weights[ix_prs_select]
  cat(length(ix_prs_select), 'selected PRS, ')
  
  if (length(ix_prs_select) == 0) {
    n_select = 0

  } else {
    beta_col = c('id', 'alt', paste0('V', ix_conv[ix_prs_select]))
    
    beta = fread(paste0(opt$out, '_beta.txt'),
                 select=beta_col, showProgress=F) %>% 
      data.frame()
    
    if (length(ix_prs_select) == 1) {
      nonzero = which(beta[,-c(1:2)]^2 > 0)
      beta_fit = data.matrix(beta[nonzero,-c(1, 2)]) * coef_glmnet
    } else {
      nonzero = which(rowSums(beta[,-c(1:2)]^2) > 0)
      beta_fit = data.matrix(beta[nonzero,-c(1, 2)]) %*% coef_glmnet
    }
    
    beta_ensemble = data.table(beta[nonzero,1:2], beta_fit)
    fwrite(beta_ensemble, paste0(opt$out, '_ensembled_beta.txt'), sep=' ')
    
    n_select = sum(beta_fit != 0)
    
    invisible(rm(beta, beta_col, beta_fit, beta_ensemble, nonzero))
  }
  
  cat(n_select, 'selected SNPs \n')
  
  invisible(rm(tun_data, best_tuning_ix, ix_conv,
               fit_glmnet, best_weights, coef_glmnet))
  invisible(gc())
}

################
## Validation ##
################
if (!is.null(opt$val)) {
  cat('Validation: ')
  # fam data
  fam = fread(paste0(opt$val, '.fam'),
              col.names=c('FID', 'IID', 'pat', 'mat', 'sex', 'pheno'),
              showProgress=F) %>%
    data.frame()
  
  # load PRS
  cat('constructing PRS, ')
  if (file.exists(paste0(opt$out, '_best_beta.txt'))) {
    run_val = system(paste0(opt$plink2, ' --bfile ', opt$val,
                            ' --score ', opt$out, '_best_beta.txt ',
                            'cols=+scoresums,-maybesid,-dosagesum,-scoreavgs ', 
                            '1 2 3 header --out ', opt$out, '_val'), intern=T)
    system(paste0('rm ', opt$out, '_val.log'))
    invisible(rm(run_val)); invisible(gc())
    
    prs_val = fread(paste0(opt$out, '_val.sscore'), showProgress=F) %>%
      rename(FID='#FID') %>%
      select(FID, IID, contains('SCORE')) %>%
      data.frame()
    names(prs_val)[3] = 'SCORE_grid'
  } else {
    prs_val = fam %>%
      select(FID, IID) %>%
      mutate(SCORE_grid = 0)
  }
  
  if (file.exists(paste0(opt$out, '_ensembled_beta.txt'))) {
    run_ensemble = system(paste0(opt$plink2, ' --bfile ', opt$val,
                                 ' --score ', opt$out, '_ensembled_beta.txt ',
                                 'cols=+scoresums,-maybesid,-dosagesum,-scoreavgs ', 
                                 '1 2 3 header --out ', opt$out, '_val_ensemble'), intern=T)
    system(paste0('rm ', opt$out, '_val_ensemble.log'))
    invisible(rm(run_ensemble)); invisible(gc())
    
    prs_val_ensemble = fread(paste0(opt$out, '_val_ensemble.sscore'), showProgress=F) %>%
      rename(FID='#FID') %>%
      select(FID, IID, contains('SCORE')) %>%
      data.frame()
    names(prs_val_ensemble)[3] = 'SCORE_ensemble'
  } else {
    prs_val_ensemble = fam %>%
      select(FID, IID) %>%
      mutate(SCORE_ensemble = 0)
  }
  
  # join PRS, phenotype, covariate
  if (!is.null(opt$cov)) {
    val_data = fam %>%
      select(FID, IID) %>%
      left_join(pheno, by=c('FID', 'IID')) %>%
      left_join(prs_val, by=c('FID', 'IID')) %>%
      left_join(prs_val_ensemble, by=c('FID', 'IID')) %>%
      left_join(cov, by=c('FID', 'IID')) %>%
      na.omit()
    
    invisible(rm(fam, prs_val, prs_val_ensemble, pheno, cov))
    invisible(gc())
    
  } else {
    val_data = fam %>%
      select(FID, IID) %>%
      left_join(pheno, by=c('FID', 'IID')) %>%
      left_join(prs_val, by=c('FID', 'IID')) %>%
      left_join(prs_val_ensemble, by=c('FID', 'IID')) %>%
      na.omit()
    
    invisible(rm(fam, prs_val, prs_val_ensemble, pheno)) 
    invisible(gc())
  }
  
  # evaluate AUC/R2
  cat('evaluating, ')
  if (binary_outcome) {
    # conditional AUC
    pred_val = roc.binary(status='pheno', variable='SCORE_grid',
                          confounders=paste0('~',cov_formula),
                          data=val_data)$auc
    
    pred_val_ensemble = roc.binary(status='pheno', variable='SCORE_ensemble',
                                   confounders=paste0('~',cov_formula),
                                   data=val_data)$auc
    
  } else {
    # adjust for covariates
    res_val = residuals(lm(paste0('pheno ~', cov_formula), data=val_data))
    
    # temporary correlation
    pred_val = cor(val_data[,'SCORE_grid'], res_val)^2
    
    pred_val_ensemble = cor(val_data[,'SCORE_ensemble'], res_val)^2
    
  }
  
  results = data.frame(
    Method = c('L0Learn-Sum', 'ALL-Sum'),
    n_snp = c(par_out %>% filter(row_number() == which.max(pred_tun)) %>% pull(L0),
              n_select),
    n_prs = c(1, length(ix_prs_select)),
    pred_val = c(pred_val, pred_val_ensemble)
  )
  
  cat('saving results \n')
  write.table(results, paste0(opt$out, '_results.txt'), 
              sep = '\t', row.names=F, col.names=T, quote=F)
  
  cat(paste0('   Prediction results written to ', opt$out, '_results.txt \n'))
}

