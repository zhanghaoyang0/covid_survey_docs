#=====================================================================================
# mr_pair | from ldsc
#=====================================================================================
source('/home/yanglab_data/user/zhanghy/project/temp/code/source_mr.r')
path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/'

pair = read.csv(paste0(path, 'result/collect/ldsc_t2d.csv'))
pair = pair[pair$rg_p<0.025,]
pair = pair[,1:2]

# add blood marker
marker = get_bloodmarker('marker')
add = rbind(cbind('spracklen_nature2020_t2d', marker), 
            cbind('spracklen_nature2020_t2d_bmisubset', marker))
colnames(add) = c('data1', 'data2')
pair = rbind(pair, add)

# two direct
pair = rbind(pair, setNames(pair[,2:1], c('data1', 'data2')))

# drop expo without 10 instrument
nsnp_clump = read.table(paste0(path, 'para/nsnp_5e-8_clumped.txt'))

keep = nsnp_clump[nsnp_clump[,2]>=10, 1]
keep = sapply(keep, function(x){gsub('.txt', '', x)})
keep = c(keep, marker, 'spracklen_nature2020_t2d', 'spracklen_nature2020_t2d_bmisubset')

pair = pair[pair$data1%in%keep, ]
pair = pair[!duplicated(pair),]

write.table(pair, paste0(path, 'para/mr_pair.txt'), row.names = F, sep = '\t', quote = F)
#=====================================================================================
# clump
#=====================================================================================
## clump.r
library(TwoSampleMR)
options(stringsAsFactors = F)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source_mr.r')
get_bbj_info()

pair = read.table(paste0(path_list[['bbj']], 'para/mr_pair.txt'), header = 1)%>%arrange(order(data1))
data1_prior = data2_prior = 'none'

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
  file_out1 = paste0(path_list[['bbj']], 'result/clump/for_clump/', pthres, '/', data1, '&', data2, '.rdata')
  file_out2 = paste0(path_list[['bbj']], 'result/clump/for_clump/', pthres, '/', data1, '&', data2, '.txt')
  # file_clump = paste0(path_list[['bbj']], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')
  # t1 = read.table(file_out2, header=1)
  # t2 = read.table(file_clump, header=1)
  # if (nrow(t2)<10){
  #   print(i);print(data1)
  #   # unlink(file_out1); unlink(file_out2); unlink(file_clump)
  # }

  if (file.exists(file_out1)&(file.exists(file_out2)|file.exists(paste0(file_out2, '.gz')))){next}
  if (data1 != data1_prior){df1_raw = read.table(paste0(path_in_1, data1, '.txt.gz'), header=1, sep = '\t')}
  if (data2 != data2_prior){df2_raw = read.table(paste0(path_in_2, data2, '.txt.gz'), header=1, sep = '\t')}

  out = get_harmo_res(df1_raw, df2_raw, pthres)
  df = out[[1]]

  save(df, file=file_out1)
  write.table(out[[2]], file_out2, row.names=F, sep='\t', quote=F)
  data1_prior = data1; data2_prior = data2
}

### clump.sh
#!/usr/bin/bash
batch="$1"
size="$2"
Rscript clump.r $batch $size

## run | 30
for i in {1..30}  
do
node=2
sbatch -N1 -n1 -c1 -w compute$node clump.sh $i 50
done

### clump
file_bfile="/home/yanglab_data/user/zhanghy/gwas/bfile/1000g/no_mhc/eas/eas_nomhc_maf_0.01"
file_pair="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/para/mr_pair.txt"
path_clump="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/result/clump/"
datas_proxy="QTL_43 QTL_89 QTL_117 Case_control_32 Case_control_36 Case_control_38 Case_control_45 Case_control_51 
 Case_control_61 Case_control_71 Case_control_75 Case_control_79 Case_control_88"
nrow=`wc -l $file_pair | awk '{print $1}'`

for i in $(seq 1 $nrow)
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
### mr.r
library(TwoSampleMR)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source_mr.r')
markers = get_bloodmarker('marker')
get_bbj_info()

# args
if (F){
  args = commandArgs(T)
  batch = as.numeric(args[1])-1
  size = as.numeric(args[2])
  range = c((size*batch+1):(size*(batch+1)))
} else {range =1:nrow(pair)}

ld = read.table(paste0(path_list[['bbj']], 'clump/instrument_r2_0.6.ld'), header=1) # remove snp overlap or ld in bi-direct test
pair = read.table(paste0(path_list[['bbj']], 'para/mr_pair.txt'), header = 1)%>%arrange(order(data1))

data1_prior = data2_prior = 'none'
range=1:nrow(pair)

for (i in range){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')

  path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})][[1]]
  path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})][[1]]
  file_in = paste0(path_list[['bbj']], 'result/clump/for_clump/', pthres, '/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[['bbj']], 'result/mr/', data1, '&', data2, '.rdata')
  file_clump = paste0(path_list[['bbj']], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')

  if (file.exists(file_out)){next}
  if (F){ # old version, clump at data1. now we clump data1 after merging data2
    path_clump = paste0(path_in_1, '../clump/clump_out/')
    snp = read.table(paste0(path_clump, data1, '.clumped'), header=1)[,3]
  }
  if (F) { # comment from t2d-blood aje paper, remove pleio snp with old version instruments
    if ((data2=='spracklen_nature2020_t2d'&data1%in%markers)|(data1=='spracklen_nature2020_t2d'&data2%in%markers)){
    snp_outcome = read.table(paste0(path_in_2, '../clump/clump_out/', data2, '.clumped'), header=1)[,3]
    snp_pleio = unique(c(snp_outcome, unlist(ld%>%filter(SNP_A%in%snp_outcome|SNP_B%in%snp_outcome)%>%select(SNP_A, SNP_B))))
    snp = snp[!snp%in%snp_pleio]
    }
  }

  load(file_in)
  snp = read.table(file_clump, header=1)[,3]
  res = get_mr_res(df, snp)

  save(res, file=file_out)
  data1_prior = data1; data2_prior = data2
  rm(df); rm(snp)
}

### mr.sh
#!/usr/bin/bash
batch="$1"
size="$2"
Rscript mr.r $batch $size

## run
for i in {0..2}  
do
node=3
sbatch -N1 -n1 -c1 -w compute$node mr.sh $i 100
done

#=====================================================================================
# MR collect | gsmr_pair
#=====================================================================================
path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/'
pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)

out = c()
for (i in 1:nrow(pair)){
  data1 = pair[i, 'data1']
  data2 = pair[i, 'data2']
  
  file_in = paste0(path, 'result/mr/', data1, '&', data2, '.rdata')
  
  load(file_in)
  
  out = c(out, unlist(res))
  if(is.na(res[1])){
    print(i)
    break()
  }
}

collect = data.frame(matrix(out, ncol = length(res), byrow=T))
colnames(collect) = names(res)
colnames(collect)[ncol(collect)] = 'p_thres'
collect[,3:26] = sapply(collect[,3:26], as.numeric)
collect = collect[,!grepl('Simple', colnames(collect))]
write.csv(collect, paste0(path, 'result/collect/mr_t2d.csv'), row.names=F)

## gsmr_pair
method = c('MR Egger_pval', 'Inverse variance weighted_pval', 'Weighted mode_pval', 'Weighted median_pval')
keep = (!is.na(collect['pleio_pval'])) & (rowSums(collect[, method]<0.05/12)>0) & (collect['Weighted mode_nsnp']>=10) # modifiy sig method==4 to > 0 for new method project
pair = collect[keep, c('data1', 'data2')]

# add blood marker
if (F){
  source('~/gwas/code/source_mr.r')
  marker = get_bloodmarker('marker')
  add = rbind(cbind('spracklen_nature2020_t2d', marker), 
              cbind('spracklen_nature2020_t2d_bmisubset', marker))
  colnames(add) = c('data1', 'data2') # two direct
  add = rbind(add, setNames(add[,2:1], c('data1', 'data2')))
  pair = rbind(pair, add)
}

pair = pair[!duplicated(pair),]

write.table(pair, paste0(path, 'para/gsmr_pair.txt'), row.names = F, sep = '\t', quote = F)

#=====================================================================================
# gsmr
# datas_proxy= c('Case_control_37_m', 'Case_control_37_f') # use proxy threshold at 1e-5
#=====================================================================================
# make_snp_list_for_r2—mat 
source('/home/yanglab_data/user/zhanghy/project/temp/code/source_mr.r')
get_bbj_info()
pair = read.table(paste0(path_list[['bbj']], 'para/mr_pair.txt'), header=1)
snp = c()
for (i in 1:nrow(pair)){
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')
  file_clump = paste0(path_list[['bbj']], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')
  snp = c(snp, read.table(file_clump, header=1)[,3])
  snp = unique(snp)
}
write.table(data.frame(snp), paste0(path_list[['bbj']], 'result/gsmr/for_gsmr/for_gsmr_snp'), row.names = F, col.names = F, quote = F)


# r2  matrix
path="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/result/gsmr/for_gsmr/"
file_bfile="/home/yanglab_data/user/zhanghy/gwas/bfile/1000g/no_mhc/eas/eas_nomhc"

plink --bfile $file_bfile \
--extract $path\for_gsmr_snp \
--r2 square \
--write-snplist \
--out $path\for_gsmr_r2_matrix

## gsmr.r
library(gsmr)
library(TwoSampleMR)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source_mr.r')
get_bbj_info()
get_gsmr_para('gsmr', 'eas') # assign setting
pair = read.table(paste0(path_list[['bbj']], 'para/mr_pair.txt'), header=1)

markers = get_bloodmarker('marker')
nsnps_thresh = 1

# ld = read.table(paste0(path_list[['bbj']], 'result/clump/instrument_r2_0.6.ld'), header=1) # remove snp overlap or ld in bi-direct test
mat = read.table(paste0(path_list[['bbj']], 'result/gsmr/for_gsmr/for_gsmr_r2_matrix.ld'))
snp_list = read.table(paste0(path_list[['bbj']], 'result/gsmr/for_gsmr/for_gsmr_r2_matrix.snplist'))[,1]
colnames(mat) = rownames(mat) = snp_list

range = 1:nrow(pair)
for (i in range){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')
  
  path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})][[1]]
  path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})][[1]]
  file_in = paste0(path_list[['bbj']], 'result/clump/for_clump/', pthres, '/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[['bbj']], 'result/gsmr/gsmr_out/', data1, '&', data2, '.rdata')
  file_clump = paste0(path_list[['bbj']], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')

  if (file.exists(file_out)){next()}
  if (F){ # old version, clump at data1. now we clump data1 after merging data2
    path_clump = paste0(path_in_1, '../clump/clump_out/')
    snp = read.table(paste0(path_clump, data1, '.clumped'), header=1)[,3]
  }
  if (F) { # comment from t2d-blood aje paper, remove pleio snp with old version instruments
    if ((data2=='spracklen_nature2020_t2d'&data1%in%markers)|(data1=='spracklen_nature2020_t2d'&data2%in%markers)){
    snp_outcome = read.table(paste0(path_in_2, '../clump/clump_out/', data2, '.clumped'), header=1)[,3]
    snp_pleio = unique(c(snp_outcome, unlist(ld%>%filter(SNP_A%in%snp_outcome|SNP_B%in%snp_outcome)%>%select(SNP_A, SNP_B))))
    snp = snp[!snp%in%snp_pleio]
    }
  }

  load(file_in)
  snp = read.table(file_clump, header=1)[,3]
  res = get_gsmr_res(df, snp, mat)
  
  save(res, file=file_out)
}

## gsmr.sh
#!/usr/bin/bash
i="$1"
Rscript gsmr.r $i

for i in {281..287}  # 1..287
do
node=$((i%10+1))
node=3
sbatch -N1 -n1 -c2 -w compute$node gsmr.sh $i
done

#=====================================================================================
# gsmr collect | cause pair
#=====================================================================================
path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/' 
pair = read.table(paste0(path, 'para/mr_pair.txt'), header=1)

out = c()
for (i in 1:nrow(pair)){
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  load(paste0(path, 'result/gsmr/gsmr_out/', data1, '&', data2, '.rdata'))
  if (length(res)==11){
    res = res[c(1,2,4,7,8,9,10, 11)]
  }
  out = c(out, unlist(res))
}

collect = data.frame(matrix(out, ncol = length(res), byrow=T))
colnames(collect) = c('data1', 'data2', 'nsnp', 'bxy', 'bxy_se', 'bxy_pval', 'heidi_flag', 'p_thres')
write.csv(collect, paste0(path, 'result/collect/gsmr_t2d.csv'), row.names=F)


## cause pair
keep = !is.na(collect$bxy_pval)&as.numeric(collect$bxy_pval)<0.05/12
pair = collect[keep, c('data1', 'data2')]

## add
source('~/gwas/code/source_mr.r')
marker = get_bloodmarker('marker')
add = rbind(cbind('spracklen_nature2020_t2d', marker), cbind('spracklen_nature2020_t2d_bmisubset', marker))
colnames(add) = c('data1', 'data2')
add = rbind(add, setNames(add[,2:1], c('data1', 'data2')))
pair = rbind(pair, add)
pair = pair[!duplicated(pair),]

write.table(pair, paste0(path, 'para/cause_pair.txt'), row.names = F, sep = '\t', quote = F)
#=====================================================================================
# cause
#=====================================================================================
## cause_prune.r
library(cause)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source_mr.r')
get_bbj_info()

pair = read.table(paste0(path_list[['bbj']], 'para/mr_pair.txt'), header=1)%>%arrange(order(data1))
data1_prior = data2_prior = 'none'

if (F){
  args = commandArgs(T)
  batch = as.numeric(args[1])-1
  size = as.numeric(args[2])
  range = c((size*batch+1):(size*(batch+1)))
} else {range=17:nrow(pair)}

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})][[1]]
  path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})][[1]]
  file_out = paste0(path_list[['bbj']], 'result/cause/for_cause/prune/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){
    rm(res)
    load(file_out)
    if (exists('res')){next}
  }
  
  if (data1 != data1_prior){df1_raw = read.table(paste0(path_in_1, data1, '.txt.gz'), header=1, sep = '\t')}
  if (data2 != data2_prior){df2_raw = read.table(paste0(path_in_2, data2, '.txt.gz'), header=1, sep = '\t')}
  
  df <- gwas_merge(df1_raw, df2_raw, snp_name_cols = c("SNP", "SNP"), beta_hat_cols = c("BETA", "BETA"), 
    se_cols = c("SE", "SE"), A1_cols = c("A1", "A1"), A2_cols = c("A2", "A2"))
  variants <- df %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
  pruned = c()

  for (chr in 1:22){
    print(paste0('chr ', chr))
    ld <- readRDS(paste0(path_list[['cause_ref']], 'eas_', chr, '_maf0.05_0.1_maxwindow_ld.RDS'))
    snp_info <- readRDS(paste0(path_list[['cause_ref']], 'eas_', chr, '_maf0.05_0.1_maxwindow_snpdata.RDS'))
    out <- ld_prune(variants = variants, ld = ld, total_ld_variants = snp_info$SNP, pval_cols = c("pval1"),pval_thresh = c(1e-3))
    pruned = c(pruned, out)
    print(length(pruned))
  }
  
  res = list()
  res$df = df
  res$pruned = pruned
  save(res, file=file_out)
  data1_prior = data1; data2_prior = data2
}

## cause_param.r
library(cause)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source_mr.r')
get_bbj_info()

pair = read.table(paste0(path_list[['bbj']], 'para/mr_pair.txt'), header=1)%>%arrange(order(data1))

if (F){
  args = commandArgs(T)
  batch = as.numeric(args[1])-1
  size = as.numeric(args[2])
  range = c((size*batch+1):(size*(batch+1)))
} else {range=1:nrow(pair)}

for (i in range){ 
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  file_in = paste0(path_list[['bbj']], 'result/cause/for_cause/prune/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[['bbj']], 'result/cause/for_cause/param/', data1, '&', data2, '.rdata')
  if (file.exists(file_out)){next()}
  load(file_in)
  set.seed(100)
  varlist <- with(res$df, sample(snp, size=1000000, replace=FALSE))
  params <- est_cause_params(res$df, varlist)
  save(params, file=file_out)
  rm(res)
}

## cause.r
library(cause)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source_mr.r')
get_bbj_info()

markers = get_bloodmarker('marker')
# remove snp overlap or ld in bi-direct test
# ld = read.table(paste0(path, 'clump/instrument_r2_0.6.ld'), header=1)
pair = read.table(paste0(path_list[['bbj']], 'para/mr_pair.txt'), header=1)%>%arrange(order(data1))

if (F){
  args = commandArgs(T)
  batch = as.numeric(args[1])-1
  size = as.numeric(args[2])
  range = c((size*batch+1):(size*(batch+1)))
} else {range=1:nrow(pair)}

for (i in range){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  file_in1 = paste0(path_list[['bbj']], 'result/cause/for_cause/prune/', data1, '&', data2, '.rdata')
  file_in2 = paste0(path_list[['bbj']], 'result/cause/for_cause/param/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[['bbj']], 'result/cause/cause_out/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){next()}
  load(file_in1); load(file_in2)
  df = res$df; pruned = res$pruned

  if (F){
    if ((data2=='spracklen_nature2020_t2d'&data1%in%markers)|
      (data1=='spracklen_nature2020_t2d'&data2%in%markers)){
    load(paste0(path, 'result/cause/for_cause/prune/', data2, '&', data1, '.rdata'))
    snp_outcome = pruned
    load(paste0(path, 'result/cause/for_cause/prune/', data1, '&', data2, '.rdata'))
    snp_pleio = unique(c(snp_outcome, unlist(ld%>%filter(SNP_A%in%snp_outcome|SNP_B%in%snp_outcome)%>%select(SNP_A, SNP_B))))
    pruned = pruned[!pruned%in%snp_pleio]
   }
  }
  res = try(cause(X=df, variants = pruned, param_ests = params))
  if ('try-error' %in% class(res)){
    res = cause(X=df, variants = pruned, param_ests = params, force=TRUE)
    res$force = 1
  } else{res$force = 0}
  save(res, file=file_out)
  rm(df); rm(pruned); rm(params)
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

for i in {51..57}  # 57
do
# node=$((i%10+1))
node=1
sbatch -N1 -n1 -c1 -w compute$node cause_prune.sh $i 10
done

## check missing pair
path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/'
col_name = c('data1', 'data2')
file_out = paste0(path, 'para/cause_pair1.txt')
write(paste(col_name, collapse = " "), file_out)
pair = read.table(paste0(path, 'para/gsmr_pair.txt'), header=1)

for (i in 1:nrow(pair)){
  if (i%%100==0){print(i)}
  data1 = pair[i, 1]; data2 = pair[i, 2]
  file = paste0(path, 'result/cause/for_cause/prune/', data1, '&',data2, '.rdata')
  if (file.exists(file)){
    rm(res); load(file)
    if (exists('res')){next}
  }
  #file = paste0(path, 'result/cause/for_cause/param/', data1, '&',data2, '.rdata')
  #file = paste0(path, 'result/cause/cause_out/', data1, '&',data2, '.rdata')
  if (!file.exists(file)){
    print(file)
    write(paste(c(data1, data2), collapse = " "), file_out, append=TRUE)
  }
}

#=====================================================================================
# cause collect
#=====================================================================================
library(cause)
path_pair = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/para/mr_pair.txt'
path_in = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/result/cause/cause_out/'
file_out = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/result/collect/cause_t2d.csv'
pair = read.table(path_pair, header=1)

out = c()

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  
  load(paste0(path_in, data1, '&', data2, '.rdata'))
  
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

write.csv(collect, file_out, row.names=F)
