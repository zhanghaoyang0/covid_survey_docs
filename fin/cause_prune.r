library(cause)
library(dplyr)

path_ld = '/home/yanglab_data/user/zhanghy/project/mr_server/db/cause_ref/eur/'
path_list = list('fin'='/home/yanglab_data/user/zhanghy/gwas/summary/finngen/clean/', 
  '_b'='/home/yanglab_data/user/zhanghy/gwas/summary/neale/clean/', 
  'xue'='/home/yanglab_data/user/zhanghy/gwas/summary/other/eur/clean/')
pair = read.table(paste0(path_list[['fin']], '../para/gsmr_pair.txt'), header=1)%>%arrange(order(data1))
datas_proxy= c('Case_control_37_m', 'Case_control_37_f') # use proxy threshold at 1e-5
data1_prior = data2_prior = 'none'

if (T){
  args = commandArgs(T)
  batch = as.numeric(args[1])-1
  size = as.numeric(args[2])
  range = c((size*batch+1):(size*(batch+1)))
} else {range=1:nrow(pair)}


for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})]
  path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})]
  file_out = paste0(path_list[['fin']], '../result/for_cause/prune/', data1, '&', data2, '.rdata')

  if (file.exists(file_out)){
    rm(res)
    load(file_out)
    if (exists('res')){next}}

  if (data1 != data1_prior){df1_raw = read.table(paste0(path_in_1, data1, '.txt.gz'), header=1, sep = '\t')}
  if (data2 != data2_prior){df2_raw = read.table(paste0(path_in_2, data2, '.txt.gz'), header=1, sep = '\t')}
  
  df <- gwas_merge(df1_raw, df2_raw, snp_name_cols = c("SNP", "SNP"), beta_hat_cols = c("BETA", "BETA"), 
    se_cols = c("SE", "SE"), A1_cols = c("A1", "A1"), A2_cols = c("A2", "A2"))
  variants <- df %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
  pruned = c()

  for (chr in 1:22){
    print(paste0('chr ', chr))
    ld <- readRDS(paste0(path_ld, 'chr', chr, '_AF0.05_0.1.RDS'))
    snp_info <- readRDS(paste0(path_ld, 'chr', chr, '_AF0.05_snpdata.RDS'))
    out <- ld_prune(variants = variants, ld = ld, total_ld_variants = snp_info$SNP, pval_cols = c("pval1"),pval_thresh = c(1e-3))
    pruned = c(pruned, out)
  }

  res = list()
  res$df = df
  res$pruned = pruned
  save(res, file=file_out)
  data1_prior = data1; data2_prior = data2
}
