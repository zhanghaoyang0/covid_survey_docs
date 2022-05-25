# gestd
#=====================================================================================
# mr_pair
# add gest diabetes in finngen
# GEST_DIABETES	finr5_263	5687	117892	Gestational diabetes (for exclusion)
#=====================================================================================
path = '/bigdat1/user/zhanghy/gwas/summary/giant/'

expos = c('giant_bmi_f_2015', 'giant_hipadjbmi_f_2015', 'giant_hip_f_2015', 'giant_wcadjbmi_f_2015', 
          'giant_wc_f_2015', 'giant_whradjbmi_f_2015', 'giant_whr_f_2015', 'giant_height_f_2013', 'giant_weight_f_2013',
          '20015_raw_f')
outcomes = c('finr5_263', 't1d_f_bolt', 't2d_f_bolt')

pair = data.frame()
for (outcome in outcomes){
  pair = rbind(pair, cbind(expos, outcome))
}

colnames(pair) = c('data1', 'data2')

# two direct
pair = rbind(pair, setNames(pair[,2:1], c('data1', 'data2')))

write.table(pair, paste0(path, 'para/mr_pair.txt'), row.names = F, sep = '\t', quote = F)
#=====================================================================================
# MR
# files_proxy= c('Case_control_37_m', 'Case_control_37_f') # use proxy threshold at 1e-5
#=====================================================================================
## mr.r
library(TwoSampleMR)
library(dplyr)
options(stringsAsFactors = F)

# res_null : for 0 nsnp (no overlap snp)
method = c('MR Egger','Weighted median','Inverse variance weighted','Simple mode','Weighted mode')
method_col = as.vector(sapply(method, function(x){paste0(x, c('_nsnp', '_b', '_se', '_pval'))}))
pleio_col = as.vector(paste0('pleio', c('_egger_intercept', '_se', '_pval')))
res_null = list()
res_null[c('exposure', 'outcome', method_col, pleio_col, "")] = NA
res_null[names(res_null)[grepl('nsnp', names(res_null))]] = 0

# args
args = commandArgs(T)
batch = as.numeric(args[1])
size = as.numeric(args[2])
range = c((size*batch+1):(size*(batch+1)))

path_pair = '/bigdat1/user/zhanghy/gwas/summary/giant/para/mr_pair.txt' # missing pair
path_in_neale = '/bigdat1/user/zhanghy/gwas/summary/neale/clean/'
path_in_fin = '/bigdat1/user/zhanghy/gwas/summary/finngen/clean/'
path_in_giant = '/bigdat1/user/zhanghy/gwas/summary/giant/clean/'
path_in_ukb = '/bigdat1/user/zhanghy/gwas/summary/ukb/eur/clean/'
path_out = '/bigdat1/user/zhanghy/gwas/summary/giant/result/mr/'

pair = read.table(path_pair, header=1)

files_proxy= c('t1d_f_bolt', 'finr5_263') # use proxy threshold at 1e-5

range = c(1:nrow(pair))

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_out = paste0(path_out, data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){
    load(file_out)
    if (as.numeric(res[['MR Egger_nsnp']])>=3 & is.na(res[['MR Egger_b']])){
      file.remove(file_out)
    } else{
      next
    }
  }
  
  if (grepl('fin', data1)){
    path_in_1 = path_in_fin
  } else if (grepl('giant', data1)){
    path_in_1 = path_in_giant
  } else if (grepl('bolt', data1)){
    path_in_1 = path_in_ukb
  } else {
    path_in_1 = path_in_neale
  }
  
  if (grepl('fin', data2)){
    path_in_2 = path_in_fin
  } else if (grepl('giant', data2)){
    path_in_2 = path_in_giant
  } else if (grepl('bolt', data2)){
    path_in_2 = path_in_ukb
  } else {
    path_in_2 = path_in_neale
  }
  
  if (data1 %in% files_proxy){
    thres = 1e-5
    snp =  read.table(paste0(path_in_1, '../clump/clump_out/1e-5/',data1, '.clumped'), header=1, sep = '')[,3]
    df1_raw = read.table(paste0(path_in_1, '../clump/for_clump/1e-5/',data1, '.txt'), header=1, sep = '')
  } else {
    thres = 5e-8
    snp =  read.table(paste0(path_in_1, '../clump/clump_out/',data1, '.clumped'), header=1, sep = '')[,3]
    df1_raw = read.table(paste0(path_in_1, '../clump/for_clump/',data1, '.txt'), header=1, sep = '')
  }
  
  df2_raw = read.table(paste0(path_in_2, data2, '.txt'), header=1, sep = '')
  
  df1 = df1_raw[df1_raw$SNP%in%snp,]
  df2 = df2_raw[df2_raw$SNP%in%snp,]
  
  if (nrow(df2)==0){
    res_null[c('exposure', 'outcome')] = c(data1, data2)
    res = res_null
    save(res, file=file_out)
    next()
  }
  
  df1 = df1 %>% rename('pval.exposure' = 'P', 'effect_allele.exposure' = 'A1', 'other_allele.exposure' = 'A2', 
                       'samplesize.exposure' = 'N', 
                       'beta.exposure' = 'BETA', 'se.exposure' = 'SE', 'eaf.exposure' = 'FRQ')
  df2 = df2 %>% rename('pval.outcome' = 'P', 'effect_allele.outcome' = 'A1', 'other_allele.outcome' = 'A2', 
                       'samplesize.outcome' = 'N', 
                       'beta.outcome' = 'BETA', 'se.outcome' = 'SE', 'eaf.outcome' = 'FRQ')
  
  df1$id.exposure = data1
  df2$id.outcome = data2
  
  df1$exposure = data1
  df2$outcome = data2
  
  dat = harmonise_data(df1, df2)
  dat = dat[dat$mr_keep==T,]
  dat = dat[order(dat$pval.exposure),]
  dat = dat[!duplicated(dat$SNP),]
  
  if (nrow(dat)<3){
    res = c(data1, data2, rep(c(nrow(dat), rep(NA, 3)), 5), rep(NA, 3), thres)
    res = as.list(res)
    
    methods = c('MR Egger_', 'Weighted median_', 'Inverse variance weighted_', 'Simple mode_', 'Weighted mode_')
    method_names = as.vector(sapply(methods, function(x){paste0(x, c('nsnp', 'b', 'se', 'pval'))}))
    names(res) = c('exposure', 'outcome', method_names, paste0('pleio_', c('egger_intercept', 'se', 'pval')), '')
  } else{
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
  }
  
  save(res, file=file_out)
  rm(list = c('df1_raw', 'df2_raw'))
}

## mr_sbatch.sh
#!/usr/bin/bash
batch="$1"
size="$2"
Rscript mr.r $batch $size


## run
for i in {91..100}  # 1..354
do
node=6
sbatch -N1 -n1 -c2 -w compute$node mr_sbatch.sh $i 50
done


#=====================================================================================
# MR collect | gsmr_pair
#=====================================================================================
path = '/bigdat1/user/zhanghy/gwas/summary/giant/'
pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)

methods = c('MR Egger', 'Weighted median', 'Inverse variance weighted', 'Weighted mode')
item1 = as.vector(sapply(methods, function(x){paste0(x, c('_nsnp', '_b', '_se', '_pval'))}))
item2 = paste0('pleio', c('_egger_intercept', '_se', '_pval'))
col_name = c('data1', 'data2', item1, item2, 'p_thres')
write(paste(col_name, collapse = ","), paste0(path, 'result/collect/mr.csv'))

out = c()
for (i in 1:nrow(pair)){
  if (i %% 100 ==0){
    print(i)
  }
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_in = paste0(path, 'result/mr/', data1, '&', data2, '.rdata')
  
  load(file_in)
  out = unlist(res)
  out = out[!grepl('Simple', names(out))]
  if(is.na(res[1])){
    print(i)
    break()
  }
  write(paste(out, collapse=","), paste0(path, 'result/collect/mr.csv'), append=TRUE)
}
#=====================================================================================
# gsmr
#=====================================================================================
# make_snp_list_for_r2——mat 
path = '/bigdat1/user/zhanghy/gwas/summary/giant/'
file_out = '/bigdat1/user/zhanghy/gwas/summary/neale/gsmr/for_gsmr_snp'
path_in_neale = '/bigdat1/user/zhanghy/gwas/summary/neale/clean/'
path_in_fin = '/bigdat1/user/zhanghy/gwas/summary/finngen/clean/'
path_in_giant = '/bigdat1/user/zhanghy/gwas/summary/giant/clean/'
path_in_ukb = '/bigdat1/user/zhanghy/gwas/summary/ukb/eur/clean/'

pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)
datas = unique(unlist(pair[,1]))

snp = c()
for (data in datas){
  if (grepl('fin', data)){
    path_clump = paste0(path_in_fin, '../clump/clump_out/')
  } else if (grepl('giant', data)){
    path_clump = paste0(path_in_giant, '../clump/clump_out/') 
  } else if (grepl('t1d', data)){
    path_clump = paste0(path_in_ukb, '../clump/clump_out/1e-5/') # proxy
  } else if (grepl('bolt', data)){
    path_clump = paste0(path_in_ukb, '../clump/clump_out/') 
  }  else {
    path_clump = paste0(path_in_neale, '../clump/clump_out/') 
  }
  
  snp = c(snp, read.table(paste0(path_clump, data, '.clumped'), header=1)[,3])
  snp = unique(snp)
}

write.table(data.frame(snp), paste0(path, 'gsmr/for_gsmr_snp'), row.names = F, col.names = F, quote = F)

# r2  matrix
path="/bigdat1/user/zhanghy/gwas/summary/giant/gsmr/"
file_bfile="/bigdat1/user/zhanghy/gwas/bfile/1000g/no_mhc/eur/eur_nomhc_maf_0.01"

plink --bfile $file_bfile \
--extract $path\for_gsmr_snp \
--r2 square \
--write-snplist \
--out $path\for_gsmr_r2_matrix

## gsmr.r
library(gsmr)
library(TwoSampleMR)
library(dplyr)
options(stringsAsFactors = F)
source('~/gwas/code/source.txt')
assign('gsmr', 'eur') # assign setting
nsnps_thresh=3 # include more gest diabetes results

args = commandArgs(T)
batch = as.numeric(args[1])
size = as.numeric(args[2])
range = c((size*batch+1):(size*(batch+1)))

path = '/bigdat1/user/zhanghy/gwas/summary/giant/' 
path_in_neale = '/bigdat1/user/zhanghy/gwas/summary/neale/clean/'
path_in_fin = '/bigdat1/user/zhanghy/gwas/summary/finngen/clean/'
path_in_giant = '/bigdat1/user/zhanghy/gwas/summary/giant/clean/'
path_in_ukb = '/bigdat1/user/zhanghy/gwas/summary/ukb/eur/clean/'

pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)

# r2 matrix
mat = read.table(paste0(path, 'gsmr/for_gsmr_r2_matrix.ld'))
snp_list = read.table(paste0(path, 'gsmr/for_gsmr_r2_matrix.snplist'))[,1]
colnames(mat) = snp_list
rownames(mat) = snp_list

files_proxy= c('t1d_f_bolt') # use proxy threshold at 1e-5
range = 1:nrow(pair)

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_out = paste0(path, 'result/gsmr/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){
    next()
  }
  
  if (grepl('fin', data1)){
    path_in_1 = path_in_fin
  } else if (grepl('giant', data1)){
    path_in_1 = path_in_giant
  } else if (grepl('bolt', data1)){
    path_in_1 = path_in_ukb
  } else {
    path_in_1 = path_in_neale
  }
  
  if (grepl('fin', data2)){
    path_in_2 = path_in_fin
  } else if (grepl('giant', data2)){
    path_in_2 = path_in_giant
  } else if (grepl('bolt', data2)){
    path_in_2 = path_in_ukb
  } else {
    path_in_2 = path_in_neale
  }
  
  if (data1 %in% files_proxy){
    thres = 1e-5
    snp =  read.table(paste0(path_in_1, '../clump/clump_out/1e-5/',data1, '.clumped'), header=1, sep = '')[,3]
    df1_raw = read.table(paste0(path_in_1, '../clump/for_clump/1e-5/',data1, '.txt'), header=1, sep = '')
  } else {
    thres = 5e-8
    snp =  read.table(paste0(path_in_1, '../clump/clump_out/',data1, '.clumped'), header=1, sep = '')[,3]
    df1_raw = read.table(paste0(path_in_1, '../clump/for_clump/',data1, '.txt'), header=1, sep = '')
  }
  
  df2_raw = read.table(paste0(path_in_2, data2, '.txt'), header=1, sep = '')
  
  df1 = df1_raw[df1_raw$SNP%in%snp,]
  df2 = df2_raw[df2_raw$SNP%in%snp,]
  
  df1 = df1 %>% rename('pval.exposure' = 'P', 'effect_allele.exposure' = 'A1', 'other_allele.exposure' = 'A2', 
                       'samplesize.exposure' = 'N', 
                       'beta.exposure' = 'BETA', 'se.exposure' = 'SE', 'eaf.exposure' = 'FRQ')
  df2 = df2 %>% rename('pval.outcome' = 'P', 'effect_allele.outcome' = 'A1', 'other_allele.outcome' = 'A2', 
                       'samplesize.outcome' = 'N', 
                       'beta.outcome' = 'BETA', 'se.outcome' = 'SE', 'eaf.outcome' = 'FRQ')
  df1[, c('id.exposure', 'exposure')] = data1
  df2[, c('id.outcome', 'outcome')] = data2
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
          gsmr_results$bxy, gsmr_results$bxy_se, gsmr_results$bxy_pval, thres)
  save(res, file=file_out)
}

## gsmr_sbatch.sh
#!/usr/bin/bash
batch="$1"
size="$2"
Rscript gsmr.r $batch $size

## run
for i in {0..10}  # 1..354
do
node=1
sbatch -N1 -n1 -c2 -w compute$node gsmr_sbatch.sh $i 20
done

## check and save miss pair
path = '/bigdat1/user/zhanghy/gwas/summary/neale/'

col_name = c('data1', 'data2')
write(paste(col_name, collapse = " "), paste0(path, '/para/gsmr_pair1.txt'))

pair = read.table(paste0(path, '/para/gsmr_pair.txt'), header=1)

for (i in 1:nrow(pair)){
  if (i%%100==0){
    print(i)
  }
  data1 = pair[i, 1]
  data2 = pair[i, 2]
  file = paste0(path, 'result/gsmr/', data1, '&', data2, '.rdata')
  if (!file.exists(file)){
    print(file)
    write(paste(c(data1, data2), collapse = " "), paste0(path, '/para/gsmr_pair1.txt'), append=TRUE)
  }
}
#=====================================================================================
# gsmr collect | cause pair
#=====================================================================================
path = '/bigdat1/user/zhanghy/gwas/summary/giant/'

pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)

out = c()
for (i in 1:nrow(pair)){
  if (i %% 100 ==0){
    print(i)
  }
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  load(paste0(path, 'result/gsmr/', data1, '&', data2, '.rdata'))
  out = c(out, unlist(res))
}

collect = data.frame(matrix(out, ncol = length(res), byrow=T))
colnames(collect) = c('data1', 'data2', 'nsnp', 'nsnp_used', 'nsnp_linkage', 
                      'nsnp_pleio', 'bxy', 'bxy_se', 'bxy_pval', 'p_thres')

write.csv(collect, paste0(path, 'result/collect/gsmr.csv'), row.names=F)
#=====================================================================================
# cause | 1.2.0
#=====================================================================================
## cause_prune.r
library(cause)
library(dplyr)

args = commandArgs(T)
batch = as.numeric(args[1])
size = as.numeric(args[2])
range = c((size*batch+1):(size*(batch+1)))

path = '/bigdat1/user/zhanghy/gwas/summary/giant/'
path_ld = '/bigdat1/user/zhanghy/project/mr_server/db/cause_ref/eur/'
path_in_neale = '/bigdat1/user/zhanghy/gwas/summary/neale/clean/'
path_in_fin = '/bigdat1/user/zhanghy/gwas/summary/finngen/clean/'
path_in_giant = '/bigdat1/user/zhanghy/gwas/summary/giant/clean/'
path_in_ukb = '/bigdat1/user/zhanghy/gwas/summary/ukb/eur/clean/'

pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)
range=48:nrow(pair)

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_out = paste0(path, 'result/for_cause/prune/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){
    next()
  }
  
  if (grepl('fin', data1)){
    path_in_1 = path_in_fin
  } else if (grepl('giant', data1)){
    path_in_1 = path_in_giant
  } else if (grepl('bolt', data1)){
    path_in_1 = path_in_ukb
  } else {
    path_in_1 = path_in_neale
  }
  
  if (grepl('fin', data2)){
    path_in_2 = path_in_fin
  } else if (grepl('giant', data2)){
    path_in_2 = path_in_giant
  } else if (grepl('bolt', data2)){
    path_in_2 = path_in_ukb
  } else {
    path_in_2 = path_in_neale
  }
  
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
  variants = variants[variants$pval1<1e-3,] # aviod a error when only snp greater than thres in ld file.
  
  pruned = c()
  
  for (chr in 1:22){
    print(paste0('chr ', chr, ' pruned snp:'))
    
    ld <- readRDS(paste0(path_ld, 'chr', chr, '_AF0.05_0.1.RDS'))
    snp_info <- readRDS(paste0(path_ld, 'chr', chr, '_AF0.05_snpdata.RDS'))
    
    out <- ld_prune(variants = variants,
                    ld = ld, total_ld_variants = snp_info$SNP,
                    pval_cols = c("pval1"),
                    pval_thresh = c(1e-3))
    
    pruned = c(pruned, out)
    print(length(pruned))
  }
  
  save(pruned, file=file_out)
}

## cause_param.r
library(cause)
library(dplyr)

args = commandArgs(T)
batch = as.numeric(args[1])
size = as.numeric(args[2])
range = c((size*batch+1):(size*(batch+1)))

path = '/bigdat1/user/zhanghy/gwas/summary/giant/'
path_in_neale = '/bigdat1/user/zhanghy/gwas/summary/neale/clean/'
path_in_fin = '/bigdat1/user/zhanghy/gwas/summary/finngen/clean/'
path_in_giant = '/bigdat1/user/zhanghy/gwas/summary/giant/clean/'
path_in_ukb = '/bigdat1/user/zhanghy/gwas/summary/ukb/eur/clean/'

pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)
range=1:nrow(pair)

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_out = paste0(path, 'result/for_cause/param/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){
    next()
  }
  
  if (grepl('fin', data1)){
    path_in_1 = path_in_fin
  } else if (grepl('giant', data1)){
    path_in_1 = path_in_giant
  } else if (grepl('bolt', data1)){
    path_in_1 = path_in_ukb
  } else {
    path_in_1 = path_in_neale
  }
  
  if (grepl('fin', data2)){
    path_in_2 = path_in_fin
  } else if (grepl('giant', data2)){
    path_in_2 = path_in_giant
  } else if (grepl('bolt', data2)){
    path_in_2 = path_in_ukb
  } else {
    path_in_2 = path_in_neale
  }
  
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

args = commandArgs(T)
batch = as.numeric(args[1])
size = as.numeric(args[2])
range = c((size*batch+1):(size*(batch+1)))

path = '/bigdat1/user/zhanghy/gwas/summary/giant/'
path_in_neale = '/bigdat1/user/zhanghy/gwas/summary/neale/clean/'
path_in_fin = '/bigdat1/user/zhanghy/gwas/summary/finngen/clean/'
path_in_giant = '/bigdat1/user/zhanghy/gwas/summary/giant/clean/'
path_in_ukb = '/bigdat1/user/zhanghy/gwas/summary/ukb/eur/clean/'

pair = read.table(paste0(path, '/para/mr_pair.txt'), header=1)
range = 1:nrow(pair)

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_out = paste0(path, 'result/cause_out/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){
    next()
  }
  
  if (grepl('fin', data1)){
    path_in_1 = path_in_fin
  } else if (grepl('giant', data1)){
    path_in_1 = path_in_giant
  } else if (grepl('bolt', data1)){
    path_in_1 = path_in_ukb
  } else {
    path_in_1 = path_in_neale
  }
  
  if (grepl('fin', data2)){
    path_in_2 = path_in_fin
  } else if (grepl('giant', data2)){
    path_in_2 = path_in_giant
  } else if (grepl('bolt', data2)){
    path_in_2 = path_in_ukb
  } else {
    path_in_2 = path_in_neale
  }
  
  # load param and pruned
  print('loading params and pruned ...')
  load(paste0(path, '/result/for_cause/prune/', data1, '&', data2, '.rdata'))
  load(paste0(path, '/result/for_cause/param/', data1, '&', data2, '.rdata'))
  print('done!')
  
  df1_raw = read.table(paste0(path_in_1, data1, '.txt'), header=1, sep = '\t')
  df2_raw = read.table(paste0(path_in_2, data2, '.txt'), header=1, sep = '\t')
  
  df <- gwas_merge(df1_raw, df2_raw, snp_name_cols = c("SNP", "SNP"), 
                   beta_hat_cols = c("BETA", "BETA"), 
                   se_cols = c("SE", "SE"), 
                   A1_cols = c("A1", "A1"), 
                   A2_cols = c("A2", "A2"))
  
  
  res = try(cause(X=df, variants = pruned, param_ests = params))
  if ('try-error' %in% class(res)){
    res = cause(X=df, variants = pruned, param_ests = params, force=TRUE)
    res$force = 1
  } else{
    res$force = 0
  }
  save(res, file=file_out)
}

## cause_prune_sbatch.sh
#!/usr/bin/bash
batch="$1"
size="$2"
Rscript cause_prune.r $batch $size

## cause_param_sbatch.sh
#!/usr/bin/bash
batch="$1"
size="$2"
Rscript cause_param.r $batch $size

## cause_sbatch.sh
#!/usr/bin/bash
batch="$1"
size="$2"
Rscript cause.r $batch $size

## run
for i in {11..15}  # 1..354
do
node=6
# sbatch -N1 -n1 -c2 -w compute$node cause_prune_sbatch.sh $i 1
# sbatch -N1 -n1 -c2 -w compute$node cause_param_sbatch.sh $i 5
sbatch -N1 -n1 -c2 -w compute$node cause_sbatch.sh $i 5
done

## check and save miss pair
path = '/bigdat1/user/zhanghy/gwas/summary/neale/'

col_name = c('data1', 'data2')
write(paste(col_name, collapse = " "), paste0(path, 'para/cause_pair1.txt'))
pair = read.table(paste0(path, 'para/cause_pair.txt'), header=1)

for (i in 1:nrow(pair)){
  if (i%%100==0){
    print(i)
  }
  data1 = pair[i, 1]
  data2 = pair[i, 2]
  file = paste0(path, 'result/cause_out/', data1, '&',data2, '.rdata')
  if (!file.exists(file)){
    print(file)
    write(paste(c(data1, data2), collapse = " "), paste0(path, 'para/cause_pair1.txt'), append=TRUE)
  }
}
#=====================================================================================
# cause collect
#=====================================================================================
library(cause)
path = '/bigdat1/user/zhanghy/gwas/summary/giant/'

pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)
skip_pairs = c('20015_raw_f&100002_raw_f', '1558_m&100022_m', '1528_f&100150_f')

out = c()

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  if (paste0(data1, '&', data2) %in% skip_pairs){
    next()
  }  
  
  load(paste0(path, 'result/cause_out/', data1, '&', data2, '.rdata'))
  
  nsnp = length(res$data$snp)
  coef = summary(res)[4]$tab
  stat = c(coef[1,3:4], coef[2, 2:4])
  out = c(out, data1, data2, nsnp, stat, summary(res)$p, res$force)
}

collect = data.frame(matrix(out, ncol = 10, byrow=T))
colnames(collect) = c('data1', 'data2', 'nsnp', 'eta_sharing', 'q_sharing', 'gamma_causal', 'eta_causal', 'q_causal', 'p', 'force')
collect$p = as.numeric(collect$p)

coef = matrix(as.numeric(unlist(strsplit(gsub('\\(|,|\\)', '', collect$gamma_causal), ' '))), ncol=3, byrow=T)

z = coef[,1]/((coef[,3]-coef[,2])/(2*1.96))
p = pnorm(abs(z), lower.tail=FALSE)*2
collect$p_gamma_causal = p
collect[is.na(collect$p_gamma_causal), 'p_gamma_causal'] = 1 # beta and z is both 0

write.csv(collect, paste0(path, 'result/collect/cause.csv'), row.names=F)
