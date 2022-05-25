library(TwoSampleMR)
library(dplyr)
options(stringsAsFactors = F)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source.r')

path_bbj = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/clean/'
path_list = list('QTL|Case'='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/clean/',
  't2d'='/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/clean/')

pair = read.table(paste0(path_bbj, '../para/mr_pair.txt'), header = 1)%>%arrange(order(data1))

datas_proxy = c('QTL_43', 'QTL_89', 'QTL_117', 'Case_control_36')
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
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')

  path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})][[1]]
  path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})][[1]]
  file_out1 = paste0(path_bbj, '../result/clump/for_clump/', pthres, '/', data1, '&', data2, '.rdata')
  file_out2 = paste0(path_bbj, '../result/clump/for_clump/', pthres, '/', data1, '&', data2, '.txt')

  if (file.exists(file_out1)&(file.exists(file_out2)|file.exists(paste0(file_out2, '.gz')))){next}
  if (data1 != data1_prior){df1_raw = read.table(paste0(path_in_1, data1, '.txt.gz'), header=1, sep = '\t')}
  if (data2 != data2_prior){df2_raw = read.table(paste0(path_in_2, data2, '.txt.gz'), header=1, sep = '\t')}

  out = get_clump_res(df1_raw, df2_raw, pthres)
  rdata = out[[1]]

  save(rdata, file=file_out1)
  write.table(out[[2]], file_out2, row.names=F, sep='\t', quote=F)
  data1_prior = data1; data2_prior = data2
}
