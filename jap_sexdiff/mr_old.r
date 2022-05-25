# i have not save the new mr code, if need, re-write accodring to that of t2d-cat
# the collect part is new
#=====================================================================================
# mr pair
#=====================================================================================
path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/'
pair = read.csv(paste0(path, 'result/collect/ldsc.csv'))
pair = pair[!grepl('ukb', pair[,1]),]

clump = read.table(paste0(path, 'para/nsnp_5e-8_clumped.txt'))
keep = c(gsub('.txt', '', clump[clump[,2]>=10, 1]), 'spracklen_nature2020_t2d_m', 'spracklen_nature2020_t2d_f')

pair = pair[as.numeric(pair$rg_p)<0.05/2, c('data1', 'data2')]
pair = rbind(pair, setNames(pair[,2:1], c('data1', 'data2')))
pair = pair[pair$data1%in%keep,]

# 87_f with <10 instruments even in 1e-5 
add = c('Case_control_87_f', 'spracklen_nature2020_t2d_f', 'Case_control_87_m', 'spracklen_nature2020_t2d_m',
        'Case_control_37_f', 'spracklen_nature2020_t2d_f', 'Case_control_37_m', 'spracklen_nature2020_t2d_m',
        'Case_control_37_f', 'Case_control_96_f', 'Case_control_37_m', 'Case_control_96_m')
add_pair = matrix(add, ncol=2, byrow=T)
add_pair = setNames(data.frame(rbind(add_pair[2,], add_pair[, c(2,1)])), c('data1', 'data2'))
pair = rbind(pair, add_pair)

pair = pair[!duplicated(pair),]
write.table(pair, paste0(path, 'para/mr_pair.txt'), row.names = F, quote = F, sep = '\t')

#=====================================================================================
# MR
#=====================================================================================
library(TwoSampleMR)
library(dplyr)
options(stringsAsFactors = F)

path_in_bbj = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/clean/'
path_in_t2d = '/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/clean/'

datas_proxy = c('Case_control_87_m', 'Case_control_87_f')

pair = read.table(paste0(path_in_bbj, '../para/mr_pair.txt'), header = 1)

range=1:nrow(pair)

for (i in range){
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_out = paste0(path_in_bbj, '../result/mr/', data1, '&', data2, '.rdata')
  if (file.exists(file_out)){
    next()
  }
  
  if (grepl('t2d', data1)){
    path_in_1 = path_in_t2d
  }  else{ path_in_1 = path_in_bbj}
  
  if (grepl('t2d', data2)){
    path_in_2 = path_in_t2d
  }  else{ path_in_2 = path_in_bbj}
  
  path_clump = paste0(path_in_1, '../clump/clump_out/')
  
  snp = read.table(paste0(path_clump, data1, '.clumped'), header=1)[,3]
  
  if (data1 %in% datas_proxy){
    thres = 1e-5
  } else{thres = 5e-8} 
  
  df1_raw = read.table(paste0(path_in_1, '../clump/for_clump/', data1, '.txt'), header=1, sep = '\t')
  df2_raw = read.table(paste0(path_in_2, data2, '.txt'), header=1, sep = '\t')
  
  df1 = df1_raw[df1_raw$SNP%in%snp,]
  df2 = df2_raw[df2_raw$SNP%in%snp,]
  
  # changed!
  df1 = df1 %>% rename('pval.exposure' = 'P', 'effect_allele.exposure' = 'A1', 'other_allele.exposure' = 'A2', 
                       'samplesize.exposure' = 'N', 
                       'beta.exposure' = 'BETA', 'se.exposure' = 'SE', 'eaf.exposure' = 'FRQ')
  df2 = df2 %>% rename('pval.outcome' = 'P', 'effect_allele.outcome' = 'A1', 'other_allele.outcome' = 'A2', 
                       'samplesize.outcome' = 'N', 
                       'beta.outcome' = 'BETA', 'se.outcome' = 'SE', 'eaf.outcome' = 'FRQ')
  
  df1[,c('id.exposure', 'exposure')] = data1
  df2[,c('id.outcome', 'outcome')] = data2
  
  dat = harmonise_data(df1, df2)
  dat = dat[dat$mr_keep==T,]
  dat = dat[order(dat$pval.exposure),]
  dat = dat[!duplicated(dat$SNP),]
  
  out = mr(dat)
  pleio = mr_pleiotropy_test(dat)[,5:7]
  colnames(pleio) = paste0('pleio_', colnames(pleio))
  
  res = c()
  for (method in out$method){
    sub = out[out$method==method, c('nsnp', 'b', 'se', 'pval')]
    colnames(sub) = sapply(colnames(sub), function(x){paste0(method, '_', x)})
    res = c(res, sub)
  }
  
  res = c(out$id.exposure[1], out$id.outcome[1], unlist(res), pleio, thres)
  
  names(res)[1:2] = c("exposure", "outcome")
  save(res, file=file_out)
}

### collect 
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
pair = read.table(paste0(path_list[['bbj']], 'para/mr_pair.txt'), header=1)
get_bbj_stratify_info()

out = c()
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1'];data2 = pair[i, 'data2']
  file_out = paste0(path_list[['bbj']], 'result/mr/', data1, '&', data2, '.rdata')
  if (!file.exists(file_out)){next}
  load(file_out)
  out = c(out, res)
}

collect = data.frame(matrix(out, ncol = length(res), byrow=T))
colnames(collect) = names(res)
collect = collect[,colnames(collect)[!grepl('Simple', colnames(collect))]]

write.csv(collect, paste0(path_list[['bbj']], 'result/collect/mr.csv') , row.names=F)
#=====================================================================================
# gsmr
# datas_proxy= c('Case_control_37_m', 'Case_control_37_f') # use proxy threshold at 1e-5
#=====================================================================================
# make_snp_list_for_r2—mat 
path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/'

path_clump_bbj = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/clump/clump_out/'
path_clump_t2d = '/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/clump/clump_out/'

pair = read.table(paste0(path, 'para/gsmr_pair.txt'), header=1)
datas = unique(unlist(pair[,1]))

snp = c()
for (data in datas){
  if (grepl('t2d', data)){
    path_clump = path_clump_t2d}
  else {
    path_clump = path_clump_bbj
  }
  snp = c(snp, read.table(paste0(path_clump, data, '.clumped'), header=1)[,3])
  snp = unique(snp)
}

write.table(data.frame(snp), paste0(path, 'gsmr/for_gsmr_snp'), row.names = F, col.names = F, quote = F)

# r2  matrix
path="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/gsmr/"
file_bfile="/home/yanglab_data/user/zhanghy/gwas/bfile/1000g/no_mhc/eas/eas_nomhc_maf_0.01"

plink --bfile $file_bfile \
--extract $path\for_gsmr_snp \
--r2 square \
--write-snplist \
--out $path\for_gsmr_r2_matrix

## gsmr.r
library(gsmr)
library(dplyr)
library(TwoSampleMR)
options(stringsAsFactors = F)
source('~/gwas/code/source.txt')
assign('gsmr', 'eas') # assign setting

path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/' 
pair = read.table(paste0(path, 'para/gsmr_pair.txt'), header=1)

# r2 matrix
mat = read.table(paste0(path, 'gsmr/for_gsmr_r2_matrix.ld'))
snp_list = read.table(paste0(path, 'gsmr/for_gsmr_r2_matrix.snplist'))[,1]
colnames(mat) = snp_list
rownames(mat) = snp_list

path_in_bbj = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/clean/'
path_in_t2d = '/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/clean/'

datas_proxy = c('Case_control_87_m', 'Case_control_87_f')

pair = read.table(paste0(path_in_bbj, '../para/gsmr_pair.txt'), header = 1)

range=1:nrow(pair)

for (i in range){
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_out = paste0(path_in_bbj, '../result/gsmr/', data1, '&', data2, '.rdata')
  if (file.exists(file_out)){
    next()
  }
  
  if (grepl('t2d', data1)){
    path_in_1 = path_in_t2d
  }  else{ path_in_1 = path_in_bbj}
  
  if (grepl('t2d', data2)){
    path_in_2 = path_in_t2d
  }  else{ path_in_2 = path_in_bbj}
  
  path_clump = paste0(path_in_1, '../clump/clump_out/')
  
  snp = read.table(paste0(path_clump, data1, '.clumped'), header=1)[,3]
  
  if (data1 %in% datas_proxy){
    thres = 1e-5
  } else{thres = 5e-8} 
  
  df1_raw = read.table(paste0(path_in_1, '../clump/for_clump/', data1, '.txt'), header=1, sep = '\t')
  df2_raw = read.table(paste0(path_in_2, data2, '.txt'), header=1, sep = '\t')
  
  df1 = df1_raw[df1_raw$SNP%in%snp,]
  df2 = df2_raw[df2_raw$SNP%in%snp,]
  
  # changed!
  df1 = df1 %>% rename('pval.exposure' = 'P', 'effect_allele.exposure' = 'A1', 'other_allele.exposure' = 'A2', 
                       'samplesize.exposure' = 'N', 
                       'beta.exposure' = 'BETA', 'se.exposure' = 'SE', 'eaf.exposure' = 'FRQ')
  df2 = df2 %>% rename('pval.outcome' = 'P', 'effect_allele.outcome' = 'A1', 'other_allele.outcome' = 'A2', 
                       'samplesize.outcome' = 'N', 
                       'beta.outcome' = 'BETA', 'se.outcome' = 'SE', 'eaf.outcome' = 'FRQ')
  
  df1[,c('id.exposure', 'exposure')] = data1
  df2[,c('id.outcome', 'outcome')] = data2
  
  dat = harmonise_data(df1, df2)
  dat = dat[dat$mr_keep==T,]
  dat = dat[order(dat$pval.exposure),]
  dat = dat[!duplicated(dat$SNP),]
  
  gsmr_data = dat%>%rename('a1'='effect_allele.exposure', 'a2'='other_allele.exposure', 
                           'bzx_pval'='pval.exposure', 'bzy_pval'='pval.outcome', 
                           'bzx'='beta.exposure', 'bzy'='beta.outcome', 
                           'bzx_se'='se.exposure', 'bzy_se'='se.outcome', 
                           'bzx_n'='samplesize.exposure', 'bzy_n'='samplesize.outcome', 
                           'a1_freq'='eaf.exposure')
  
  ldrho = mat[rownames(mat) %in% gsmr_data$SNP, colnames(mat) %in% gsmr_data$SNP]
  snp_coeff_id = rownames(ldrho)
  
  gsmr_results = try(gsmr(gsmr_data$bzx, gsmr_data$bzx_se, gsmr_data$bzx_pval, 
                          gsmr_data$bzy, gsmr_data$bzy_se, gsmr_data$bzy_pval,
                          ldrho, snp_coeff_id, n_ref, 
                          heidi_outlier_flag, gwas_thresh = thres, 
                          single_snp_heidi_thresh, multi_snp_heidi_thresh, 
                          nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, 
                          gsmr2_beta))
  
  if (grepl('non-conformable arrays', gsmr_results[[1]])){
    gsmr_results = try(gsmr(gsmr_data$bzx, gsmr_data$bzx_se, gsmr_data$bzx_pval, 
                            gsmr_data$bzy, gsmr_data$bzy_se, gsmr_data$bzy_pval,
                            ldrho, snp_coeff_id, n_ref, 
                            F, gwas_thresh = thres, 
                            single_snp_heidi_thresh, multi_snp_heidi_thresh, 
                            nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, 
                            gsmr2_beta))
    heidi = 1
  } else {
    heidi = 0
  }
  
  res = c(data1, data2, nrow(gsmr_data), length(gsmr_results$used_index),length(gsmr_results$linkage_snps),length(gsmr_results$pleio_snps),
          gsmr_results$bxy, gsmr_results$bxy_se, gsmr_results$bxy_pval, heidi, thres)
  
  save(res, file=file_out)
}

### collect 
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_bbj_stratify_info()
pair = read.table(paste0(path_list[['bbj']], 'para/mr_pair.txt'), header=1)
out = c()
for (i in 1:nrow(pair)){
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  file_in = paste0(path_list[['bbj']], 'result/gsmr/gsmr_out/', data1, '&', data2, '.rdata')
  load(file_in)
  if (length(res)==11){res = res[-c(3,5,6)]}
  out = c(out, res)
}
collect = data.frame(matrix(out, ncol = length(res), byrow=T))
colnames(collect) = c('data1', 'data2', 'gsmr_nsnp', 'gsmr_b', 'gsmr_se', 'gsmr_pval', 'heidi', 'pthres')

write.csv(collect, paste0(path_list[['bbj']], 'result/collect/gsmr.csv') , row.names=F)
#=====================================================================================
# cause
#=====================================================================================
## cause_prune.r
library(cause)
library(dplyr)

path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/'
path_ld = '/home/yanglab_data/user/zhanghy/project/mr_server/db/cause_ref/eas/'
path_in_bbj = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/clean/'
path_in_t2d = '/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/clean/'

pair = read.table(paste0(path, 'para/cause_pair.txt'), header=1)

range = 1:nrow(pair)

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_out = paste0(path, 'result/for_cause/prune/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){
    next()
  }
  
  if (grepl('t2d', data1)){
    path_in_1 = path_in_t2d
  }  else{ path_in_1 = path_in_bbj}
  
  if (grepl('t2d', data2)){
    path_in_2 = path_in_t2d
  }  else{ path_in_2 = path_in_bbj}
  
  df1_raw = read.table(paste0(path_in_1, data1, '.txt'), header=1, sep = '\t')
  df2_raw = read.table(paste0(path_in_2, data2, '.txt'), header=1, sep = '\t')
  
  df <- gwas_merge(df1_raw, df2_raw, snp_name_cols = c("SNP", "SNP"), 
                   beta_hat_cols = c("BETA", "BETA"), 
                   se_cols = c("SE", "SE"), 
                   A1_cols = c("A1", "A1"), 
                   A2_cols = c("A2", "A2"))
  
  rm(list = c('df1_raw', 'df2_raw'))
  # prune
  print('pruning ...')
  variants <- df %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
  
  pruned = c()
  
  for (chr in 1:22){
    print(paste0('chr ', chr, ' pruned snp:'))
    
    ld <- readRDS(paste0(path_ld, 'eas_', chr, '_maf0.05_0.1_maxwindow_ld.RDS'))
    snp_info <- readRDS(paste0(path_ld, 'eas_', chr, '_maf0.05_0.1_maxwindow_snpdata.RDS'))
    
    out <- ld_prune(variants = variants,
                    ld = ld, total_ld_variants = snp_info$SNP,
                    pval_cols = c("pval1"),
                    pval_thresh = c(1e-3))
    #pval_thresh = c(p_thres))
    pruned = c(pruned, out)
    print(length(pruned))
  }
  
  save(pruned, file=file_out)
}

## cause_param.r
library(cause)
library(dplyr)

path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/'
path_in_bbj = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/clean/'
path_in_t2d = '/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/clean/'

pair = read.table(paste0(path, 'para/cause_pair.txt'), header=1)

range = 1:nrow(pair)

for (i in range){ # 1:29
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_out = paste0(path, 'result/for_cause/param/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){
    next()
  }
  
  if (grepl('t2d', data1)){
    path_in_1 = path_in_t2d
  }  else{ path_in_1 = path_in_bbj}
  
  if (grepl('t2d', data2)){
    path_in_2 = path_in_t2d
  }  else{ path_in_2 = path_in_bbj}
  
  df1_raw = read.table(paste0(path_in_1, data1, '.txt'), header=1, sep = '\t')
  df2_raw = read.table(paste0(path_in_2, data2, '.txt'), header=1, sep = '\t')
  
  df <- gwas_merge(df1_raw, df2_raw, snp_name_cols = c("SNP", "SNP"), 
                   beta_hat_cols = c("BETA", "BETA"), 
                   se_cols = c("SE", "SE"), 
                   A1_cols = c("A1", "A1"), 
                   A2_cols = c("A2", "A2"))
  
  rm(list = c('df1_raw', 'df2_raw'))
  
  set.seed(100)
  varlist <- with(df, sample(snp, size=1000000, replace=FALSE))
  params <- est_cause_params(df, varlist)
  
  save(params, file=file_out)
}

## cause.r
library(cause)
library(dplyr)

path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/'
path_in_bbj = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/clean/'
path_in_t2d = '/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/clean/'

pair = read.table(paste0(path, 'para/cause_pair.txt'), header=1)

range = 1:nrow(pair)

for (i in range){ # 1:29
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_out = paste0(path, 'result/cause_out/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){
    next()
  }
  
  if (grepl('t2d', data1)){
    path_in_1 = path_in_t2d
  }  else{ path_in_1 = path_in_bbj}
  
  if (grepl('t2d', data2)){
    path_in_2 = path_in_t2d
  }  else{ path_in_2 = path_in_bbj}
  
  # load param and pruned
  print('loading params and pruned ...')
  load(paste0(path, 'result/for_cause/prune/', data1, '&', data2, '.rdata'))
  load(paste0(path, 'result/for_cause/param/', data1, '&', data2, '.rdata'))
  print('done!')
  
  df1_raw = read.table(paste0(path_in_1, data1, '.txt'), header=1, sep = '\t')
  df2_raw = read.table(paste0(path_in_2, data2, '.txt'), header=1, sep = '\t')
  
  df <- gwas_merge(df1_raw, df2_raw, snp_name_cols = c("SNP", "SNP"), 
                   beta_hat_cols = c("BETA", "BETA"), 
                   se_cols = c("SE", "SE"), 
                   A1_cols = c("A1", "A1"), 
                   A2_cols = c("A2", "A2"))
  
  rm(list = c('df1_raw', 'df2_raw'))
  
  res = try(cause(X=df, variants = pruned, param_ests = params))
  if ('try-error' %in% class(res)){
    res = cause(X=df, variants = pruned, param_ests = params, force=TRUE)
    res$force = 1
  } else{
    res$force = 0
  }
  save(res, file=file_out)
}

### collect
library(cause)
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_bbj_stratify_info()
pair = read.table(paste0(path_list[['bbj']], 'para/mr_pair.txt'), header=1)

out = c()
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  file_in = paste0(path_list[['bbj']], 'result/cause/cause_out/', data1, '&', data2, '.rdata')
  load(file_in)
  print(data1); print(data2)
  print(summary(res))
  out = c(out, data1, data2, res$stat, res$force)
}

collect = data.frame(matrix(unlist(out), ncol = length(c(data1, data2, res$stat, res$force)), byrow=T))
colnames(collect) = c('data1', 'data2', names(res$stat), 'force')

write.csv(collect, paste0(path_list[['bbj']], 'result/collect/cause.csv') , row.names=F)


#=====================================================================================
# collect eas female t2d bmi
#=====================================================================================
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_bbj_stratify_info()
path = path_list[['bbj']]
## gsmr
gsmr = read.csv(paste0(path, 'result/collect/t2d_cat/gsmr.csv'))
gsmr = gsmr%>%rename(nsnp=gsmr_nsnp, b=gsmr_b, se=gsmr_se, p=gsmr_pval)%>%mutate(method="GSMR")%>%select(method, data1, data2, nsnp, nsnp_x, b, se, p, r2_data1, r2_data2, pthres)

## mr
df = read.csv(paste0(path, 'result/collect/t2d_cat/mr.csv'))
mr = data.frame()
for (method in c('MR.Egger', 'Inverse.variance.weighted', 'Weighted.median', 'Weighted.mode')){
  sub = df[,c(names(df)[grepl(method, names(df))], 'nsnp_x', 'data1', 'data2', 'snp_r2.exposure', 'snp_r2.outcome', 'pthres')]%>%mutate(t=method)
  names(sub) = c('nsnp', 'b_raw', 'se', 'p', 'nsnp_x', 'data1', 'data2', 'r2_data1', 'r2_data2', 'pthres', 'method')
  mr = rbind(mr, sub)
}

## mrlap
mrlap = read.csv(paste0(path, 'result/collect/t2d_cat/mrlap.csv'))%>%mutate(method='MRlap')
mrlap = mrlap%>%merge(mr%>%filter(method=='Inverse.variance.weighted')%>%select(data1, data2, nsnp, nsnp_x, r2_data1, r2_data2, pthres), by=c('data1', 'data2'))

## cause
cause = read.csv(paste0(path, 'result/collect/t2d_cat/cause.csv'))
cause = cause%>%rename(nsnp=cause_nsnp, b=cause_causal_gamma_b, se=cause_causal_gamma_se, p=cause_causal_gamma_pval, r2_data1=snp_r2.exposure, r2_data2=snp_r2.outcome)%>%
  mutate(method="CAUSE", pthres=1e-3)%>%select(method, data1, data2, nsnp, nsnp_x, b, se, p, r2_data1, r2_data2, pthres)

## combine
res = data.frame(rbind(cause, gsmr, mr, mrlap))
res = res%>%mutate(b_l=b-1.96*se, b_u=b+1.96*se)%>%mutate(b_raw=b, b_l_raw=b_l, b_u_raw=b_u)
## logit to liability
res[res$data1 == 'E11_m', c('b_raw', 'b_l', 'b_u')] = res[res$data1 == 'E11_m', c('b_raw', 'b_l', 'b_u')]*get_lia_factor(prev_t2d_m_eur, prev_cataract_m_eur)
res[res$data1 == 'E11_f', c('b_raw', 'b_l', 'b_u')] = res[res$data1 == 'E11_f', c('b_raw', 'b_l', 'b_u')]*get_lia_factor(prev_t2d_f_eur, prev_cataract_f_eur)
res[res$data2 == 'E11_m', c('b_raw', 'b_l', 'b_u')] = res[res$data2 == 'E11_m', c('b_raw', 'b_l', 'b_u')]*get_lia_factor(prev_cataract_m_eur, prev_t2d_m_eur)
res[res$data2 == 'E11_f', c('b_raw', 'b_l', 'b_u')] = res[res$data2 == 'E11_f', c('b_raw', 'b_l', 'b_u')]*get_lia_factor(prev_cataract_f_eur, prev_t2d_f_eur)

res[res$data1 == 'Case_control_37_m', c('b_raw', 'b_l', 'b_u')] = res[res$data1 == 'Case_control_37_m', c('b_raw', 'b_l', 'b_u')]*get_lia_factor(prev_cataract_m_eas, prev_t2d_m_eas)
res[res$data1 == 'Case_control_37_f', c('b_raw', 'b_l', 'b_u')] =  res[res$data1 == 'Case_control_37_f', c('b_raw', 'b_l', 'b_u')]*get_lia_factor(prev_cataract_f_eas, prev_t2d_f_eas)
res[res$data2 == 'Case_control_37_m', c('b_raw', 'b_l', 'b_u')] =  res[res$data2 == 'Case_control_37_m', c('b_raw', 'b_l', 'b_u')]*get_lia_factor(prev_t2d_m_eas, prev_cataract_m_eas)
res[res$data2 == 'Case_control_37_f', c('b_raw', 'b_l', 'b_u')] =  res[res$data2 == 'Case_control_37_f', c('b_raw', 'b_l', 'b_u')]*get_lia_factor(prev_t2d_f_eas, prev_cataract_f_eas)

res[,c('or', 'or_l', 'or_u')] = exp(res[,c('b_raw', 'b_l', 'b_u')])
res[,c('or_raw', 'or_l_raw', 'or_u_raw')] = exp(res[,c('b_raw', 'b_l_raw', 'b_u_raw')])
res[res=='MR.Egger'] = 'MR-Egger'; res[res=='Inverse.variance.weighted'] = 'IVW'; res[res=='Weighted.median'] = 'Weighted median'; res[res=='Weighted.mode'] = 'Weighted mode'
write.csv(res, '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/collect/t2d_cat/collect.csv', row.names=F)

res = res%>%filter(method!='MRlap')
t2d_m_eur = res%>%filter(data1=='E11_m')%>%pull('b_raw')
t2d_f_eur = res%>%filter(data1=='E11_f')%>%pull('b_raw')
cat_m_eur = res%>%filter(data1=='cataract_m')%>%pull('b_raw')
cat_f_eur = res%>%filter(data1=='cataract_f')%>%pull('b_raw')

t2d_m_agen = res%>%filter(data1=='spracklen_nature2020_t2d_m'&data2=='Case_control_37_m')%>%pull('b_raw')
t2d_f_agen = res%>%filter(data1=='spracklen_nature2020_t2d_f'&data2=='Case_control_37_f')%>%pull('b_raw')
cat_m_agen = res%>%filter(data1=='Case_control_37_m'&data2=='spracklen_nature2020_t2d_m')%>%pull('b_raw')
cat_f_agen = res%>%filter(data1=='Case_control_37_f'&data2=='spracklen_nature2020_t2d_f')%>%pull('b_raw')

t2d_m_bbj = res%>%filter(data1=='Case_control_96_m'&data2=='Case_control_37_m')%>%pull('b_raw')
t2d_f_bbj = res%>%filter(data1=='Case_control_96_f'&data2=='Case_control_37_f')%>%pull('b_raw')
cat_m_bbj = res%>%filter(data1=='Case_control_37_m'&data2=='Case_control_96_m')%>%pull('b_raw')
cat_f_bbj = res%>%filter(data1=='Case_control_37_f'&data2=='Case_control_96_f')%>%pull('b_raw')

t.test(t2d_m_agen, t2d_f_agen, paired=T)
t.test(cat_m_agen, cat_f_agen, paired=T)

t.test(t2d_m_bbj, t2d_f_bbj, paired=T)
t.test(cat_m_bbj, cat_f_bbj, paired=T)

t.test(t2d_m_eur, t2d_f_eur, paired=T)
t.test(cat_m_eur, cat_f_eur, paired=T)