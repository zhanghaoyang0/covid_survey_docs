# finngen
#=====================================================================================
# mr_pair
#=====================================================================================
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/gwas/summary/finngen/'

pair = read.csv(paste0(path, 'result/collect/ldsc.csv'))%>%filter(rg_p<0.001)
pair = pair[,1:2]

# two direct
pair = rbind(pair, setNames(pair[,2:1], c('data1', 'data2')))

# drop expo without 10 instrument
nsnp_clump_fin = read.table(paste0(path, 'para/nsnp_5e-8_clumped.txt'))
nsnp_clump_neale = read.table(paste0(path, '../neale/para/nsnp_5e-8_clumped.txt'))
nsnp_clump = rbind(nsnp_clump_fin, nsnp_clump_neale)

keep = nsnp_clump[nsnp_clump[,2]>=10, 1]
keep = sapply(keep, function(x){gsub('.txt', '', x)})
keep = c(keep, 'xueal_nc2018_t2d')

pair = pair[pair$data1%in%keep&!grepl('covid|ecg', pair$data2), ]

write.table(pair, paste0(path, 'para/mr_pair.txt'), row.names = F, sep = '\t', quote = F)

#=====================================================================================
# clump
#=====================================================================================
## clump.r
library(TwoSampleMR)
library(dplyr)
options(stringsAsFactors = F)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source.r')

path_pair = '/home/yanglab_data/user/zhanghy/gwas/summary/finngen/para/mr_pair.txt'
path_list = list('fin'='/home/yanglab_data/user/zhanghy/gwas/summary/finngen/clean/',
  'covid'='/home/yanglab_data/user/zhanghy/gwas/summary/covid/clean/',
  'ecg'='/home/yanglab_data/user/zhanghy/gwas/summary/ecg/verweij/clean/',
  '_b'='/home/yanglab_data/user/zhanghy/gwas/summary/neale/clean/',
  'xue'='/home/yanglab_data/user/zhanghy/gwas/summary/other/eur/clean/')
pair = read.table(path_pair, header=1)%>%arrange(order(data1))

datas_proxy= c('Case_control_37_m', 'Case_control_37_f') # use proxy threshold at 1e-5
data1_prior = data2_prior = 'none'

# args
if (F){
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
  file_out1 = paste0(path_list[['fin']], '../result/clump/for_clump/', pthres, '/', data1, '&', data2, '.rdata')
  file_out2 = paste0(path_list[['fin']], '../result/clump/for_clump/', pthres, '/', data1, '&', data2, '.txt')

  if (file.exists(file_out1)&(file.exists(file_out2)|file.exists(paste0(file_out2, '.gz')))){next}
  if (data1 != data1_prior){df1_raw = read.table(paste0(path_in_1, data1, '.txt.gz'), header=1, sep = '\t')}
  if (data2 != data2_prior){df2_raw = read.table(paste0(path_in_2, data2, '.txt.gz'), header=1, sep = '\t')}

  out = get_clump_res(df1_raw, df2_raw, pthres)
  df = out[[1]]

  save(df, file=file_out1)
  write.table(out[[2]], file_out2, row.names=F, sep='\t', quote=F)
  data1_prior = data1; data2_prior = data2
}

## clump.sh
#!/usr/bin/bash
batch="$1"
size="$2"
Rscript clump.r $batch $size

## run | 34
for i in {331..345}  
do
node=6
sbatch -N1 -n1 -c1 -w compute$node clump.sh $i 100
done

### clump
file_bfile="/home/yanglab_data/user/zhanghy/gwas/bfile/1000g/no_mhc/eur/eur_nomhc_maf_0.01"
file_pair="/home/yanglab_data/user/zhanghy/gwas/summary/finngen/para/mr_pair.txt"
path_clump="/home/yanglab_data/user/zhanghy/gwas/summary/finngen/result/clump/"
datas_proxy=""
nrow=`wc -l $file_pair | awk '{print $1}'`

for i in $(seq 2 $nrow)
do
echo $i
data1=`awk -v i=$i 'NR==i {print $1}' $file_pair`
data2=`awk -v i=$i 'NR==i {print $2}' $file_pair`

if [ `echo $datas_proxy | grep -w $data1 | wc -l` == 1 ]
then
pthres="1e-5"
else
pthres="5e-8"
fi

file_in=${path_clump}for_clump/${pthres}/${data1}\&${data2}.txt.gz
file_out=${path_clump}clump_out/${pthres}/${data1}\&${data2}

if [ -f $file_out.clumped ]
then
continue
fi

plink --allow-no-sex \
--clump $file_in \
--bfile $file_bfile \
--clump-kb 1000 --clump-p1 1.0 --clump-p2 1.0 --clump-r2 0.05 \
--out $file_out
done
#=====================================================================================
# MR
#=====================================================================================
## mr.r
library(TwoSampleMR)
library(dplyr)
options(stringsAsFactors = F)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source.r')

if (F){
  args = commandArgs(T)
  batch = as.numeric(args[1])-1
  size = as.numeric(args[2])
  range = c((size*batch+1):(size*(batch+1)))
} else {range=1:nrow(pair)}

path_list = list('fin'='/home/yanglab_data/user/zhanghy/gwas/summary/finngen/clean/',
    'covid'='/home/yanglab_data/user/zhanghy/gwas/summary/covid/clean/',
    'ecg'='/home/yanglab_data/user/zhanghy/gwas/summary/ecg/verweij/clean/',
    '_b'='/home/yanglab_data/user/zhanghy/gwas/summary/neale/clean/',
    'xue'='/home/yanglab_data/user/zhanghy/gwas/summary/other/eur/clean/')
pair = read.table(paste0(path_list[['fin']], '../para/mr_pair.txt'), header=1)%>%arrange(order(data1))

datas_proxy= c('Case_control_37_m', 'Case_control_37_f') # use proxy threshold at 1e-5
data1_prior = data2_prior = 'none'
range = 1:nrow(pair)

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')

  path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})]
  path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})]
  file_in = paste0(path_list[['fin']], '../result/clump/for_clump/', pthres, '/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[['fin']], '../result/mr/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){next}
  
  file_clump = paste0(path_list[['fin']], '../result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')
  snp = read.table(file_clump, header=1)[,3]
  load(file_in)
  res = get_mr_res(df, snp)
  
  save(res, file=file_out)
  rm(df); rm(snp)
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
path = '/home/yanglab_data/user/zhanghy/gwas/summary/finngen/'

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
path = '/home/yanglab_data/user/zhanghy/gwas/summary/finngen/'
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

## gsmr_pair | 
collect = read.csv(paste0(path, 'result/collect/mr.csv'))
method = c('MR.Egger_pval', 'Inverse.variance.weighted_pval', 'Weighted.mode_pval', 'Weighted.median_pval')
keep = (!is.na(collect['pleio_pval'])) & (rowSums(collect[, method]<0.05/12)>0) & (collect['Weighted.mode_nsnp']>=10)  # modify sig method==4 to >0 for new method project.
pair = collect[keep, c('exposure', 'outcome')]
colnames(pair) = c('data1', 'data2')

write.table(pair, paste0(path, 'para/gsmr_pair.txt'), row.names = F, sep = '\t', quote = F) 
#=====================================================================================
# gsmr
# files_proxy= c('Case_control_37_m', 'Case_control_37_f') # use proxy threshold at 1e-5
#=====================================================================================
# make_snp_list_for_r2—mat 
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/gwas/summary/finngen/'
pair = read.table(paste0(path, 'para/gsmr_pair.txt'), header=1)%>%arrange(order(data1))
datas = unique(unlist(pair[,1]))
datas_proxy= c('Case_control_37_m', 'Case_control_37_f') # use proxy threshold at 1e-5

out = c()
for (i in 1:nrow(pair)){
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')
  file_clump = paste0(path, 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')
  snp = read.table(file_clump, header=1)[,3]
  out = c(out, snp)
  out = unique(out)
}

write.table(data.frame(out), paste0(path, 'gsmr/for_gsmr_snp'), row.names = F, col.names = F, quote = F)

# r2  matrix
path="/home/yanglab_data/user/zhanghy/gwas/summary/finngen/gsmr/"
file_bfile="/home/yanglab_data/user/zhanghy/gwas/bfile/1000g/no_mhc/eur/eur_nomhc"

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
source('/home/yanglab_data/user/zhanghy/project/temp/code/source.r')
assign('gsmr', 'eur') # assign setting

if (F){
  args = commandArgs(T)
  batch = as.numeric(args[1])-1
  size = as.numeric(args[2])
  range = c((size*batch+1):(size*(batch+1)))
} else {range=1:nrow(pair)}

path = '/home/yanglab_data/user/zhanghy/gwas/summary/finngen/' 
pair = read.table(paste0(path, 'para/gsmr_pair.txt'), header=1)%>%arrange(order(data1))
mat = read.table(paste0(path, 'gsmr/for_gsmr_r2_matrix.ld'))
snp_list = read.table(paste0(path, 'gsmr/for_gsmr_r2_matrix.snplist'))[,1]
colnames(mat) = rownames(mat) = snp_list

datas_proxy= c('Case_control_37_m', 'Case_control_37_f') # use proxy threshold at 1e-5

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')

  file_in = paste0(path, 'result/clump/for_clump/', pthres, '/', data1, '&', data2, '.rdata')
  file_out = paste0(path, 'result/gsmr/', data1, '&', data2, '.rdata')
  file_clump = paste0(path, 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')
  
  if (file.exists(file_out)){next}

  load(file_in)
  snp = read.table(file_clump, header=1)[,3]
  res = get_gsmr_res(df, snp, mat)
  
  save(res, file=file_out)
}

## gsmr_sbatch.sh
#!/usr/bin/bash
batch="$1"
size="$2"
Rscript gsmr.r $batch $size

## run
for i in {271..290}  
do
node=8
sbatch -N1 -n1 -c2 -w compute$node gsmr_sbatch.sh $i 50
done

## check and save miss pair
path = '/home/yanglab_data/user/zhanghy/gwas/summary/finngen/'

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
path = '/home/yanglab_data/user/zhanghy/gwas/summary/finngen/'

pair = read.table(paste0(path, 'para/gsmr_pair.txt'), header=1)

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

## cause pair
keep = !is.na(collect$bxy_pval)&as.numeric(collect$bxy_pval)<0.05/12
pair = collect[keep, c('data1', 'data2')]
write.table(pair, paste0(path, 'para/cause_pair.txt'), row.names = F, sep = '\t', quote = F)

#=====================================================================================
# cause
#=====================================================================================
## cause_prune.r
library(cause)
library(dplyr)

path_ld = '/home/yanglab_data/user/zhanghy/project/mr_server/db/cause_ref/eur/'
path_list = list('fin'='/home/yanglab_data/user/zhanghy/gwas/summary/finngen/clean/', 
  '_b'='/home/yanglab_data/user/zhanghy/gwas/summary/neale/clean/', 
  'xue'='/home/yanglab_data/user/zhanghy/gwas/summary/other/eur/clean/')
pair = read.table(paste0(path_list[['fin']], '../para/gsmr_pair.txt'), header=1)%>%arrange(order(data1))
datas_proxy= c('Case_control_37_m', 'Case_control_37_f') # use proxy threshold at 1e-5
data1_prior = data2_prior = 'none'

if (F){
  args = commandArgs(T); batch = as.numeric(args[1])-1; size = as.numeric(args[2]); range = c((size*batch+1):(size*(batch+1)))
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

## cause_param.r
library(cause)
library(dplyr)

path = '/bigdat1/user/zhanghy/gwas/summary/finngen/'  
pair = read.table(paste0(path, 'para/gsmr_pair.txt'), header=1)

if (F){
  args = commandArgs(T); batch = as.numeric(args[1])-1; size = as.numeric(args[2]); range = c((size*batch+1):(size*(batch+1)))
} else {range=1:nrow(pair)}

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  file_in = paste0(path, 'result/for_cause/prune/', data1, '&', data2, '.rdata')
  file_out = paste0(path, 'result/for_cause/param/', data1, '&', data2, '.rdata')
  if (file.exists(file_out)){next()}
  load(file_in)
  if (!exists('res')){
    next
  }
  set.seed(100)
  varlist <- with(res$df, sample(snp, size=1000000, replace=FALSE))
  params <- est_cause_params(res$df, varlist)
  save(params, file=file_out)
  rm(res)
}

## cause.r
library(cause)
library(dplyr)

if (F){
  args = commandArgs(T)
  batch = as.numeric(args[1])-1
  size = as.numeric(args[2])
  range = c((size*batch+1):(size*(batch+1)))
} else {range=1:nrow(pair)}

path_ld = '/home/yanglab_data/user/zhanghy/gwas/plink_file/cause_ref/eur/'
path_list = list('fin'='/home/yanglab_data/user/zhanghy/gwas/summary/finngen/clean/',
  '_b'='/home/yanglab_data/user/zhanghy/gwas/summary/neale/clean/',
  'xue'='/home/yanglab_data/user/zhanghy/gwas/summary/other/eur/clean/')

pair = read.table(paste0(path, 'para/gsmr_pair.txt'), header=1)
range = 1:nrow(pair)

skip_pair = c('1697_b&100001_raw_b', '1697_b&100002_raw_b', '20015_raw_b&100002_raw_b')

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  if (paste0(data1, '&', data2)%in%skip_pair){
    next
  }
  
  file_out = paste0(path, 'result/cause_out/', data1, '&', data2, '.rdata')
  if (file.exists(file_out)){
    next()
  }
  
  path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})]
  path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})]
  
  # load param and pruned
  print('loading params and pruned ...')
  load(paste0(path, 'result/for_cause/prune/', data1, '&', data2, '.rdata'))
  load(paste0(path, 'result/for_cause/param/', data1, '&', data2, '.rdata'))
  print('done!')
  
  df1_raw = read.table(paste0(path_in_1, data1, '.txt.gz'), header=1, sep = '\t')
  df2_raw = read.table(paste0(path_in_2, data2, '.txt.gz'), header=1, sep = '\t')
  
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

## cause_prune.sh
#!/usr/bin/bash
batch="$1"
size="$2"
Rscript cause_prune.r $batch $size

## cause_param.sh
#!/usr/bin/bash
batch="$1"
size="$2"
Rscript cause_param.r $batch $size

## cause.sh
#!/usr/bin/bash
batch="$1"
size="$2"
Rscript cause.r $batch $size

## run
for i in {191..202}  
do
# node=$((i%10+1))
node=3
sbatch -N1 -n1 -c2 -w compute$node cause_prune.sh $i 50
# sbatch -N1 -n1 -c2 -w compute$node cause_param.sh $i 50
# sbatch -N1 -n1 -c2 -w compute$node cause.sh $i 5
done


## check and save miss pair
path = '/home/yanglab_data/user/zhanghy/gwas/summary/finngen/'

col_name = c('data1', 'data2')
write(paste(col_name, collapse = " "), paste0(path, 'para/cause_pair3.txt'))
pair = read.table(paste0(path, 'para/cause_pair.txt'), header=1)

for (i in 1:nrow(pair)){
  if (i%%100==0){
    print(i)
  }
  data1 = pair[i, 1]
  data2 = pair[i, 2]
  #file = paste0(path, 'result/for_cause/param/', data1, '&',data2, '.rdata')
  file = paste0(path, 'result/cause_out/', data1, '&',data2, '.rdata')
  if (!file.exists(file)){
    print(file)
    write(paste(c(data1, data2), collapse = " "), paste0(path, 'para/cause_pair3.txt'), append=TRUE)
  }
}

#=====================================================================================
# cause collect
#=====================================================================================
library(cause)
path = '/home/yanglab_data/user/zhanghy/gwas/summary/finngen/'

pair = read.table(paste0(path, 'para/cause_pair.txt'), header=1)

skip_pair = c('1697_b&100001_raw_b', '1697_b&100002_raw_b', '20015_raw_b&100002_raw_b')

out = c()
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file = paste0(path, 'result/cause_out/', data1, '&', data2, '.rdata')
  
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