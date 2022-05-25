#=====================================================================================
# mr_pair
#=====================================================================================
require(data.table)
library(readxl)

path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/'
file_info_gwascata = paste0(path, 'info/gwascatalog_info.xlsx') 

info_gwascata = data.frame(read_excel(file_info_gwascata, sheet=1))
info_gwascata = info_gwascata[is.na(info_gwascata$drop),]
h2_gwascata = read.table(paste0(path, 'para/h2.txt'), sep = '\t', header = 1)
h2 = rbind(h2_gwascata)
keep1 = gsub('.txt', '', h2[h2$h2>0&h2$p<0.05, 'data'])

nsnp_clump_gwascata = read.table(paste0(path, 'para/nsnp_5e-8_clumped.txt'))
nsnp_clump = rbind(nsnp_clump_gwascata)
keep2 = nsnp_clump[nsnp_clump[,2]>=10, 1]
keep2 = sapply(keep2, function(x){gsub('.txt', '', x)})

pair = expand.grid(info_gwascata[info_gwascata$group=='marker', 'label'], info_gwascata[info_gwascata$group=='disease', 'label'])
pair = rbind(pair, setNames(pair[,c(2:1)], c('Var1', 'Var2')))

pair = pair[pair[,1]%in%keep1 & pair[,2]%in%keep1,] # both h2 p < 0.05
pair = pair[pair[,1]%in%keep2,] # exposure with >=10 instruments
names(pair) = c('data1', 'data2')

write.table(pair, paste0(path, 'para/mr_pair.txt'), sep='\t', quote=F, row.names=F)
#=====================================================================================
# MR
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

# args
args = commandArgs(T)
batch = as.numeric(args[1])
size = as.numeric(args[2])
range = c((size*batch+1):(size*(batch+1)))

path_pair = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/para/mr_pair.txt'
path_in_gwascata = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/clean/'
path_out = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/result/mr/'

pair = read.table(path_pair, header=1)

files_proxy= c('Case_control_37_m', 'Case_control_37_f') # use proxy threshold at 1e-5

range = 1:nrow(pair)

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_out = paste0(path_out, data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){
    load(file_out)
    if (length(res)==0){
      file.remove(file_out)
    } else if (!is.na(res[[1]])){
      next()
    }
    rm(res)
  }
  
  if (grepl('GCST', data1)){
    path_in_1 = path_in_gwascata
  } else {
    path_in_1 = path_in_gwascata
  }
  
  if (grepl('GCST', data2)){
    path_in_2 = path_in_gwascata
  } else {
    path_in_2 = path_in_gwascata
  }
  
  path_clump = paste0(path_in_1, '../clump/clump_out/5e-8/')
  snp = read.table(paste0(path_clump, data1, '.clumped'), header=1)[,3]
  
  if (data1 %in% files_proxy){
    thres = 1e-5
  } else {thres = 5e-8}
  
  df1_raw = read.table(paste0(path_in_1, '../clump/for_clump/5e-8/', data1, '.txt'), header=1, sep = '\t')
  df2_raw = read.table(paste0(path_in_2, data2, '.txt'), header=1, sep = '\t')
  
  df1 = df1_raw[df1_raw$SNP%in%snp,]
  df2 = df2_raw[df2_raw$SNP%in%snp,]
  
  if (nrow(df2)<=3){
    res_null[c('exposure', 'outcome')] = c(data1, data2)
    res = res_null
    res[names(res)[grepl('nsnp', names(res))]] = nrow(df2)
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
  dat = dat[dat$se.outcome&dat$se.exposure!=0,] # se=0 make egger error
  
  rm(list = c('df1_raw', 'df2_raw'))
  
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
  
  names(res)[1:2] = c("data1", "data2")
  save(res, file=file_out)
}

## mr_sbatch.sh
#!/usr/bin/bash
batch="$1"
size="$2"
Rscript mr.r $batch $size


## run
for i in {81..100}  # 1..354
do
node=10
sbatch -N1 -n1 -c2 -w compute$node mr_sbatch.sh $i 100
done


## check and save miss pair
path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/'

col_name = c('data1', 'data2')
write(paste(col_name, collapse = " "), paste0(path, 'para/mr_pair1.txt'))
pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)

for (i in 1:nrow(pair)){
  if (i%%100==0){
    print(i)
  }
  file1 = pair[i, 1]
  file2 = pair[i, 2]
  file = paste0(path, 'result/mr/', file1, '&',file2, '.rdata')
  if (!file.exists(file)){
    print(file)
    write(paste(c(file1, file2), collapse = " "), paste0(path, 'para/mr_pair1.txt'), append=TRUE)
  }
}

#=====================================================================================
# MR collect | gsmr_pair
#=====================================================================================
path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/'
pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)

methods = c('MR Egger', 'Weighted median', 'Inverse variance weighted', 'Weighted mode')
item1 = as.vector(sapply(methods, function(x){paste0(x, c('_nsnp', '_b', '_se', '_pval'))}))
item2 = paste0('pleio', c('_egger_intercept', '_se', '_pval'))
col_name = c('exposure', 'outcome', item1, item2, 'p_thres')
write(paste(col_name, collapse = ","), paste0(path, 'result/collect/mr.csv'))

for (i in 1:nrow(pair)){
  if (i %% 100 ==0){
    print(i)
  }
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_in = paste0(path, 'result/mr/', data1, '&', data2, '.rdata')
  
  # if (!file.exists(file_in)){
  #   next() # remove next!!!
  # }
  
  load(file_in)
  out = unlist(res)
  out = out[!grepl('Simple', names(out))]
  if(is.na(res[1])){
    print('error !!!')
    break()
  }
  write(paste(out, collapse=","), paste0(path, 'result/collect/mr.csv'), append=TRUE)
}

# ## gsmr_pair
# collect = read.csv(paste0(path, 'result/collect/mr.csv'))
# method = c('MR.Egger_pval', 'Inverse.variance.weighted_pval', 'Weighted.mode_pval', 'Weighted.median_pval')
# keep = (!is.na(collect['pleio_pval'])) & (rowSums(collect[, method]<0.05/12)==4) & (collect['Weighted.mode_nsnp']>=10) 
# pair = collect[keep, c('exposure', 'outcome')]
# colnames(pair) = c('data1', 'data2')
# 
# write.table(pair, paste0(path, 'para/gsmr_pair.txt'), row.names = F, sep = '\t', quote = F)
#=====================================================================================
# gsmr
#=====================================================================================
# make_snp_list_for_r2—mat 
path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/'

path_clump_gwascata = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/clump/clump_out/5e-8/'

pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)
datas = unique(unlist(pair[,1]))

snp = c()
for (data in datas){
  if (grepl('GCST', data)){
    path_clump = path_clump_gwascata
  } else {
    path_clump = path_clump_gwascata
  }
  snp = c(snp, read.table(paste0(path_clump, data, '.clumped'), header=1)[,3])
  snp = unique(snp)
}

write.table(data.frame(snp), paste0(path, 'gsmr/for_gsmr_snp'), row.names = F, col.names = F, quote = F)

# r2  matrix
path="/bigdat1/user/zhanghy/gwas/summary/gwascatalog/gsmr/"
file_bfile="/bigdat1/user/zhanghy/gwas/bfile/1000g/no_mhc/eur/eur_nomhc"

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
assign('gsmr', 'eur') # assign setting
nsnps_thresh = 1

args = commandArgs(T)
batch = as.numeric(args[1])
size = as.numeric(args[2])
range = c((size*batch+1):(size*(batch+1)))

path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/' 
path_in_gwascata = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/clean/'

pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)

# r2 matrix
mat = read.table(paste0(path, 'gsmr/for_gsmr_r2_matrix.ld'))
snp_list = read.table(paste0(path, 'gsmr/for_gsmr_r2_matrix.snplist'))[,1]
colnames(mat) = snp_list
rownames(mat) = snp_list

files_proxy= c('Case_control_37_m', 'Case_control_37_f') # use proxy threshold at 1e-5

range = 1:nrow(pair)
for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_out = paste0(path, 'result/gsmr/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){
    load(file_out)
    if (!is.na(res[[1]])){
      next()
    }
    rm(res)
  }
  
  if (grepl('GCST', data1)){
    path_in_1 = path_in_gwascata
  } else {
    path_in_1 = path_in_gwascata
  }
  
  if (grepl('GCST', data2)){
    path_in_2 = path_in_gwascata
  } else {
    path_in_2 = path_in_gwascata
  }
  
  path_clump = paste0(path_in_1, '../clump/clump_out/5e-8/')
  snp = read.table(paste0(path_clump, data1, '.clumped'), header=1)[,3]
  
  if (data1 %in% files_proxy){
    thres = 1e-5
  } else {thres = 5e-8}
  
  df1_raw = read.table(paste0(path_in_1, '../clump/for_clump/5e-8/', data1, '.txt'), header=1, sep = '\t')
  df2_raw = read.table(paste0(path_in_2, data2, '.txt'), header=1, sep = '\t')
  
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
for i in {171..180}  # 1..354
do
node=10
sbatch -N1 -n1 -c2 -w compute$node gsmr_sbatch.sh $i 20
done

## check and save miss pair
path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/'

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
path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/'

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

# ## cause pair
# keep = !is.na(collect$bxy_pval)&as.numeric(collect$bxy_pval)<0.05/12
# pair = collect[keep, c('data1', 'data2')]
# write.table(pair, paste0(path, 'para/cause_pair.txt'), row.names = F, sep = '\t', quote = F)

#=====================================================================================
# cause
#=====================================================================================
## cause_prune.r
library(cause)
library(dplyr)

args = commandArgs(T)
batch = as.numeric(args[1])
size = as.numeric(args[2])
range = c((size*batch+1):(size*(batch+1)))

path = '/home/yanglab_data/user/zhanghy/gwas/summary/gwascatalog/'
path_ld = '/home/yanglab_data/user/zhanghy/project/mr_server/db/cause_ref/eur/'
path_in_gwascata = '/home/yanglab_data/user/zhanghy/gwas/summary/gwascatalog/clean/'

pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)

range = 1:nrow(pair)

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_out = paste0(path, 'result/cause/for_cause/prune/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){
    res = try(load(file_out))
    if (!'try-error' %in% class(res)){
      next
    }
  }
  
  if (grepl('GCST', data1)){
    path_in_1 = path_in_gwascata
  } else {
    path_in_1 = path_in_gwascata
  }
  
  if (grepl('GCST', data2)){
    path_in_2 = path_in_gwascata
  } else {
    path_in_2 = path_in_gwascata
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
  
  pruned = c()
  
  for (chr in 1:22){
    print(paste0('chr ', chr, ' pruned snp:'))
    
    ld <- readRDS(paste0(path_ld, 'chr', chr, '_AF0.05_0.1.RDS'))
    snp_info <- readRDS(paste0(path_ld, 'chr', chr, '_AF0.05_snpdata.RDS'))
    
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

# args = commandArgs(T)
# batch = as.numeric(args[1])
# size = as.numeric(args[2])
# range = c((size*batch+1):(size*(batch+1)))

path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/'
path_in_gwascata = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/clean/'

pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)

range = 1:nrow(pair)

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  file_out = paste0(path, 'result/cause/for_cause/param/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){
    res = try(load(file_out))
    if (!'try-error' %in% class(res)){
      next
    }
  }
  
  if (grepl('GCST', data1)){
    path_in_1 = path_in_gwascata
  } else {
    path_in_1 = path_in_gwascata
  }
  
  if (grepl('GCST', data2)){
    path_in_2 = path_in_gwascata
  } else {
    path_in_2 = path_in_gwascata
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

# args = commandArgs(T)
# batch = as.numeric(args[1])
# size = as.numeric(args[2])
# range = c((size*batch+1):(size*(batch+1)))

path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/'
path_ld = '/bigdat1/user/zhanghy/gwas/plink_file/cause_ref/eur/'
path_in_gwascata = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/clean/'

pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)
range = 1:nrow(pair)

skip_pair = c('1697_b&100001_raw_b', '1697_b&100002_raw_b', '20015_raw_b&100002_raw_b')

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  if (paste0(data1, '&', data2)%in%skip_pair){
    next
  }
  
  file_out = paste0(path, 'result/cause/cause_out/', data1, '&', data2, '.rdata')
  if (file.exists(file_out)){
    next()
  }
  
  if (grepl('GCST', data1)){
    path_in_1 = path_in_gwascata
  } else {
    path_in_1 = path_in_gwascata
  }
  
  if (grepl('GCST', data2)){
    path_in_2 = path_in_gwascata
  } else {
    path_in_2 = path_in_gwascata
  }
  
  # load param and pruned
  print('loading params and pruned ...')
  load(paste0(path, 'result/cause/for_cause/prune/', data1, '&', data2, '.rdata'))
  load(paste0(path, 'result/cause/for_cause/param/', data1, '&', data2, '.rdata'))
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
for i in {0..15}  
do
node=$((i%10+1))

sbatch -N1 -n1 -c2 -w compute$node cause_prune_sbatch.sh $i 1
# sbatch -N1 -n1 -c2 -w compute$node cause_param_sbatch.sh $i 2
# sbatch -N1 -n1 -c2 -w compute$node cause_sbatch.sh $i 2
done

## check and save miss pair
path = '/home/yanglab_data/user/zhanghy/gwas/summary/gwascatalog/'

col_name = c('data1', 'data2')
file_out = paste0(path, 'para/mr_pair1.txt')

write(paste(col_name, collapse = " "), file_out)
pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)

for (i in 1:nrow(pair)){
  if (i%%100==0){
    print(i)
  }
  data1 = pair[i, 1]
  data2 = pair[i, 2]
  file = paste0(path, 'result/cause/for_cause/prune/', data1, '&',data2, '.rdata')
  # file = paste0(path, 'result/cause_out/', data1, '&',data2, '.rdata')
  if (!file.exists(file)){
    print(file)
    write(paste(c(data1, data2), collapse = " "), file_out, append=TRUE)
  }
}

#=====================================================================================
# cause collect
#=====================================================================================
library(cause)
path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/'

pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)

skip_pair = c('1697_b&100001_raw_b', '1697_b&100002_raw_b', '20015_raw_b&100002_raw_b')

out = c()
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file = paste0(path, 'result/cause/cause_out/', data1, '&', data2, '.rdata')
  
  if (paste0(data1, '&', data2)%in%skip_pair){
    next
  }
  
  if (! file.exists(file)){
    next()
  }
  
  load(file)
  
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

write.csv(collect, paste0(path, 'result/collect/cause.csv'), row.names=F)





#=====================================================================================
# collect all | scatter plot
#=====================================================================================
library(dplyr)
library(readxl)
library(stringr)

path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/'
file_info = paste0(path, 'info/gwascatalog_info.xlsx') 

info = data.frame(read_excel(file_info, sheet=1))

pair = read.table(paste0(path, 'para/mr_pair.txt'), sep='\t', header=1)
mr = read.csv(paste0(path, 'result/collect/mr.csv'))
colnames(mr)[1:2] = c('data1', 'data2')

gsmr = read.csv(paste0(path, 'result/collect/gsmr.csv'))
colnames(gsmr)[3:9] = paste0('gsmr_', colnames(gsmr)[3:9])
gsmr$p_thres=mr$p_thres=NULL
df = mr%>%merge(gsmr, by=c('data1', 'data2'))

cause = read.csv(paste0(path, 'result/collect/cause.csv'))
cause$force = NULL
cause = cause%>%rename(cause=gamma_causal, cause_pval=p_gamma_causal, cause_p_compare=p, cause_nsnp=nsnp)
cols = c('eta_sharing', 'q_sharing', 'cause', 'eta_causal', 'q_causal')
for (col in cols){
  sub = str_split_fixed(gsub(',', '', gsub(')', '', gsub('(', '', cause[,col], fixed=T))), ' ', 3)
  colnames(sub) = paste0(col, c('_b', '_b_l', '_b_u'))
  cause = cbind(cause, sub)
  cause[col] = NULL
}
df = df%>%merge(cause, by=c('data1', 'data2'))

df$data1_mean = df$data1
df$data2_mean = df$data2

for (i in c(unique(c(df$data1, df$data2)))){
  df[,c('data1_mean', 'data2_mean')][df[,c('data1_mean', 'data2_mean')]==i] = info[info$label==i, 'trait']
}


write.csv(df, paste0(path, 'result/collect/collect.csv'), row.names=F)


df = read.csv(paste0(path, 'result/collect/ldsc.csv'))
df$data1_mean = df$data1
df$data2_mean = df$data2

for (i in c(unique(c(df$data1, df$data2)))){
  df[,c('data1_mean', 'data2_mean')][df[,c('data1_mean', 'data2_mean')]==i] = info[info$label==i, 'trait']
}
write.csv(df, paste0(path, 'result/collect/ldsc_mean.csv'), row.names=F)