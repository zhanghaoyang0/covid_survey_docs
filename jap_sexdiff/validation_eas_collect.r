# eas bbj, spracklen
#=====================================================================================
# clean agen t2d
#=====================================================================================
library(dplyr)
library(stringr)

path = '/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/'
file_bim = '/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/chn.bim'
bim = read.table(file_bim)[, c(1, 4, 2)]%>%rename(SNP=V2, CHR=V1, POS=V4)%>%filter(!(CHR==6&POS<=33448354&POS>=28477797))

files = list.files(paste0(path, 'raw/'))
files = files[grepl('_m|_f', files)]

for (file in files){
  print(file)
  df = read.table(gzfile(paste0(path, 'raw/', file)), header=1, fill=T)
  df = df%>%rename(CHR=Chr, POS=Pos, A1=EA, A2=NEA, FRQ=EAF, N=Neff, BETA=Beta)%>%merge(bim, by = c('CHR', 'POS'))%>%select(CHR, POS, SNP, FRQ, A1, A2, BETA, SE, P, N)
  df = df%>%mutate(A1=toupper(A1), A2=toupper(A2), P=as.numeric(P))%>%na.omit()%>%filter(FRQ>0.01&FRQ<0.99)
  df = df%>%mutate(id=1:nrow(df))%>%arrange(P)%>%filter(!duplicated(SNP))%>%arrange(id)%>%select(-id)
  write.table(df, paste0(path, '/clean/', gsub('.gz', '', file)), row.names=F, quote=F, sep='\t')
}
#=====================================================================================
# ldsc
#=====================================================================================
prev_t2d_m_pop=0.086; prev_t2d_f_pop=0.068
prev_cataract_m_pop=0.406; prev_cataract_f_pop=0.4233

prev_t2d_m_sam_sprack=0.2388549; prev_t2d_f_sam_sprack=0.1685085
prev_t2d_m_sam_bbj=0.2369583; prev_t2d_f_sam_bbj=0.1420604
prev_cataract_m_sam=0.1065; prev_cataract_f_sam=0.1259

trait1=Case_control_96_f
trait2=Case_control_37_f

if [[ $trait1 = *"_m_"* ]]; then sex="m"; else sex="f"; fi
if [[ $trait1 = *"sprack"* ]]
then file1="/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/munge/${trait1}.sumstats.gz"
else file1="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/munge/${trait1}.sumstats.gz"
fi
file2="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/munge/${trait2}.sumstats.gz"

path_ld="/home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/${sex}/"
file_out="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ldsc/t2d_cat/${trait1}&${trait2}"

/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/ldsc.py \
--rg  $file1,$file2 \
--ref-ld-chr $path_ld \
--w-ld-chr $path_ld \
--samp-prev $prev_t2d_f_sam_bbj,$prev_cataract_f_sam \
--pop-prev $prev_t2d_f_pop,$prev_cataract_f_pop \
--out ${file_out}_lia

/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/ldsc.py \
--rg  $file1,$file2 \
--ref-ld-chr $path_ld \
--w-ld-chr $path_ld \
--samp-prev $prev_t2d_f_sam_bbj,$prev_cataract_f_sam \
--pop-prev $prev_t2d_f_pop,$prev_cataract_f_pop \
--intercept-h2 1,1 \
--out ${file_out}_constrain_lia


## compare rg with z test, in ng. use this test, risk of eur female>male, risk of eas male > female
rg1 = 0.18
rg2 = 0.26
se1 = 0.06
se2 = 0.06

z = (rg1-rg2)/(se1^2+se2^2)
p = 2*pnorm(abs(z), lower.tail=FALSE)

z
p
#=====================================================================================
# hetero
#=====================================================================================
for sex in m f 
do for pop in eas eur 
do
trait1=spracklen_nature2020_t2d_${sex}
trait2=E11_${sex}

# trait1=Case_control_37_${sex}
# trait2=cataract_${sex}

if [[ $pop == "eas" ]]; then path_ld="/home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/${sex}/"
else path_ld="/home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr/${sex}/"; fi

if [[ $trait1 = *"sprack"* ]]; then file1="/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/munge/${trait1}.sumstats.gz"
else file1="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/munge/${trait1}.sumstats.gz"; fi
file2="/home/yanglab_data/user/zhanghy/gwas/summary/ukb/eur/munge/${trait2}.sumstats.gz"
file_out="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ldsc/t2d_cat/${trait1}&${trait2}_ref_${pop}"

# /home/yanglab_data/user/zhanghy/soft_slurm/ldsc/ldsc.py \
# --rg  $file1,$file2 \
# --ref-ld-chr $path_ld \
# --w-ld-chr $path_ld \
# --out ${file_out}

done;done

#=====================================================================================
# clump
#=====================================================================================
## clump.r
library(TwoSampleMR)
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_t2dcat_diff_info('eas')
data1_prior = data2_prior = 'none'

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')

  path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})][[1]]
  path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})][[1]]

  file_out1 = paste0(path_list[['bbj']], 'result/clump/for_clump/', pthres, '/', data1, '&', data2, '.rdata')
  file_out2 = paste0(path_list[['bbj']], 'result/clump/for_clump/', pthres, '/', data1, '&', data2, '.txt')

  if (file.exists(file_out1)&(file.exists(file_out2)|file.exists(paste0(file_out2, '.gz')))){next}
  if (data1 != data1_prior){df1_raw = read.table(paste0(path_in_1, 'clean/', data1, '.txt.gz'), header=1, sep = '\t')}
  if (data2 != data2_prior){df2_raw = read.table(paste0(path_in_2, 'clean/', data2, '.txt.gz'), header=1, sep = '\t')}

  out = get_harmo_res(df1_raw, df2_raw, pthres)
  df = out[[1]]

  save(df, file=file_out1)
  write.table(out[[2]], file_out2, row.names=F, sep='\t', quote=F)
  data1_prior = data1; data2_prior = data2
}

### clump
path_clump="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/clump/"
datas_proxy="Case_control_37_f Case_control_37_m"

for sex in m f
do
file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/chn_${sex}"

data1="spracklen_nature2020_bmiadjt2d_${sex}"
#data2="Case_control_96_${sex}"
data2="Case_control_37_${sex}"

if [ `echo $datas_proxy | grep -w $data1 | wc -l` == 1 ]
then pthres="1e-5"
else pthres="5e-8"
fi

file_in=${path_clump}for_clump/${pthres}/${data1}\&${data2}.txt.gz
file_out=${path_clump}clump_out/${pthres}/${data1}\&${data2}

if [ -f $file_out.clumped ]
then continue
fi

plink --allow-no-sex \
--clump $file_in \
--bfile $file_bfile \
--clump-kb 1000 --clump-p1 1.0 --clump-p2 1.0 --clump-r2 0.05 \
--out $file_out

done
#=====================================================================================
# ld of instruments, remove for bi-dir test
#=====================================================================================
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_t2dcat_diff_info('eas', F)

out = c()
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')
  snp1 = read.table(paste0(path_list[[1]], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped'), header=1)[,3]
  load(paste0(path_list[['bbj']], 'result/cause/for_cause/prune/t2d_cat/', data1, '&', data2, '.rdata'))
  snp2 = res$pruned
  out = c(out, snp1, snp2)
}
out = unique(out)

write.table(data.frame(out), paste0(path_list[[1]], 'result/ld/for_ld/t2d-cat_eas'), row.names = F, col.names = F, quote = F)

## r2 
for sex in m f
do
file_snp="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ld/for_ld/t2d-cat_eas"
file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/chn_${sex}"
file_out="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ld/ld_out/t2d-cat_${sex}_eas_0.6"

plink --bfile $file_bfile \
--extract $file_snp \
--r2  \
--ld-window-r2 0.6 \
--out $file_out
done
#=====================================================================================
# MR
#=====================================================================================
### mr.r
library(TwoSampleMR)
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_t2dcat_diff_info('eas')

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1'];data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')
  sex = ifelse(grepl('_m', data1), 'm', 'f')

  file_in = paste0(path_list[[1]], 'result/clump/for_clump/', pthres, '/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[[1]], 'result/mr/t2d_cat/', data1, '&', data2, '.rdata')
  file_clump = paste0(path_list[[1]], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')

  if (file.exists(file_out)){next()}
  load(file_in)
  snp_pleio = get_snp_pleio(data2, data1, path_list[[1]], 'mr', ld[[sex]]) # remove for bi-test
  snp = read.table(file_clump, header=1)[,3]
  sex = ifelse(grepl('_m', data1), 'm', 'f')
  snp = snp[!snp%in%snp_pleio]

  res = get_mr_res(df, snp)
  
  df$r.exposure = get_r_from_lor(df$beta.exposure, df$eaf.exposure, get_r2_para(data1)$ncase, get_r2_para(data1)$nctrl, get_r2_para(data1)$prev) 
  df$r.outcome = get_r_from_lor(df$beta.outcome, df$eaf.outcome, get_r2_para(data2)$ncase, get_r2_para(data2)$nctrl, get_r2_para(data2)$prev) 
  r2 = directionality_test(df)[,5:8]
  
  res[['nsnp_x']] = sum(snp%in%bim[,2])
  res = c(res, r2); res = unlist(res)
  save(res, file=file_out)
}

# collect.r
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_t2dcat_diff_info('all')

out = c()
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1'];data2 = pair[i, 'data2']
  file_out = paste0(path_list[['bbj']], 'result/mr/t2d_cat/', data1, '&', data2, '.rdata')
  load(file_out)
  out = c(out, res)
}

collect = data.frame(matrix(out, ncol = length(res), byrow=T))
colnames(collect) = names(res)
collect = collect[,colnames(collect)[!grepl('Simple', colnames(collect))]]

write.csv(collect, paste0(path_list[['bbj']], 'result/collect/t2d_cat/mr.csv') , row.names=F)
#=====================================================================================
# gsmr
#=====================================================================================
## make_snp_list_for_extract.r 
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_t2dcat_diff_info('eas', F)

out = c()
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1'];data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')

  file_clump = paste0(path_list[[1]], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')
  snp = read.table(file_clump, header=1)[,3]
  out = c(out, snp)
}
out = unique(out)
write.table(data.frame(out), paste0(path_list[[1]], 'result/gsmr/for_gsmr/r2_mat_t2d_cat_snp'), row.names = F, col.names = F, quote = F)

## r2  matrix
for sex in m f
do
file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/chn_${sex}"

plink --bfile $file_bfile \
--extract /home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/gsmr/for_gsmr/r2_mat_t2d_cat_snp \
--r2 square \
--write-snplist \
--out /home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/gsmr/for_gsmr/r2_mat_t2d_cat_${sex}

done

### gsmr.r
library(gsmr)
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_gsmr_para('gsmr', 'eas') # get_gsmr_para setting
get_t2dcat_diff_info('eas')

# r2 matrix
mat = list()
for (sex in c('m', 'f')){
  mat[[sex]] = read.table(paste0(path_list[[1]], 'result/gsmr/for_gsmr/', 'r2_mat_t2d_cat_', sex, '.ld'))
  colnames(mat[[sex]]) = rownames(mat[[sex]]) = read.table(paste0(path_list[[1]], 'result/gsmr/for_gsmr/', 'r2_mat_t2d_cat_', sex, '.snplist'))[,1]
}

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')
  gwas_thresh = as.numeric(pthres)
  sex = ifelse(grepl('_m', data1), 'm', 'f')

  file_in = paste0(path_list[[1]], 'result/clump/for_clump/', pthres, '/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[[1]], 'result/gsmr/gsmr_out/t2d_cat/', data1, '&', data2, '.rdata')
  file_clump = paste0(path_list[[1]], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')
  load(file_in)

  snp_pleio = get_snp_pleio(data2, data1, path_list[[1]], 'mr', ld[[sex]]) # remove for bi-test
  snp = read.table(file_clump, header=1)[,3]
  snp = snp[!snp%in%snp_pleio]
  res = get_gsmr_res(df, snp, mat[[sex]])

  df$r.exposure = get_r_from_lor(df$beta.exposure, df$eaf.exposure, get_r2_para(data1)$ncase, get_r2_para(data1)$nctrl, get_r2_para(data1)$prev) 
  df$r.outcome = get_r_from_lor(df$beta.outcome, df$eaf.outcome, get_r2_para(data2)$ncase, get_r2_para(data2)$nctrl, get_r2_para(data2)$prev) 
  r2 = directionality_test(df)[,5:8]
  res[['nsnp_x']] = sum(snp%in%bim[,2])
  res = c(res, r2); res = unlist(res)
  save(res, file=file_out)
}

# collect.r
get_t2dcat_diff_info('all')
out = c()
for (i in 1:nrow(pair)){
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  file_in = paste0(path_list[['bbj']], 'result/gsmr/gsmr_out/t2d_cat/', data1, '&', data2, '.rdata')
  load(file_in)
  out = c(out, res)
}
collect = data.frame(matrix(out, ncol = length(res), byrow=T))
colnames(collect) = names(res)

write.csv(collect, paste0(path_list[['bbj']], 'result/collect/t2d_cat/gsmr.csv') , row.names=F)
#=====================================================================================
# cause
# move prune and param from original analysis
#=====================================================================================
# cause_prune.r 
library(cause)
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_t2dcat_diff_info('eas', F)
data1_prior = data2_prior = 'none'
range = 1:nrow(pair)
for (i in 10){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  sex = ifelse(grepl('_m', data1), 'm', 'f')
  path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})][[1]]
  path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})][[1]]
  file_out = paste0(path_list[['bbj']], 'result/cause/for_cause/prune/t2d_cat/', data1, '&', data2, '.rdata')

  if (file.exists(file_out)){next}

  if (data1 != data1_prior){df1_raw = read.table(paste0(path_in_1, 'clean/', data1, '.txt.gz'), header=1, sep = '\t')}
  if (data2 != data2_prior){df2_raw = read.table(paste0(path_in_2, 'clean/', data2, '.txt.gz'), header=1, sep = '\t')}

  df <- gwas_merge(df1_raw, df2_raw, snp_name_cols = c("SNP", "SNP"), beta_hat_cols = c("BETA", "BETA"), 
    se_cols = c("SE", "SE"), A1_cols = c("A1", "A1"), A2_cols = c("A2", "A2"))
  variants <- df %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
  pruned = c()

  for (chr in 1:23){
    print(paste0('chr ', chr, ' pruned snp:'))
    ld <- readRDS(paste0(path_list[['ref']], 'chn_', sex, '_', chr, '_ld.RDS'))
    snp_info <- readRDS(paste0(path_list[['ref']], 'chn_', sex, '_', chr, '_snp.RDS'))
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

# cause_param.r
library(cause)
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_t2dcat_diff_info('eas', F)

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  file_in = paste0(path_list[[1]], 'result/cause/for_cause/prune/t2d_cat/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[[1]], 'result/cause/for_cause/param/t2d_cat/', data1, '&', data2, '.rdata')

  if (file.exists(file_out)){next()}
  load(file_in)
  set.seed(100)
  varlist <- with(res$df, sample(snp, size=1000000, replace=FALSE))
  params <- est_cause_params(res$df, varlist)
  save(params, file=file_out)
  rm(res)
}

### cause.r
library(TwoSampleMR)
library(cause)
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_t2dcat_diff_info('eas')

data1_prior = data2_prior = 'none'
for (i in 10){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  file_in1 = paste0(path_list[[1]], 'result/cause/for_cause/prune/t2d_cat/', data1, '&', data2, '.rdata')
  file_in2 = paste0(path_list[[1]], 'result/cause/for_cause/param/t2d_cat/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[[1]], 'result/cause/cause_out/t2d_cat/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){next()}
  load(file_in1); load(file_in2)

  snp_pleio = get_snp_pleio(data2, data1, path_list[[1]], 'cause', ld[[sex]]) # remove for bi-test
  df = res$df; pruned = res$pruned; pruned = pruned[!pruned%in%snp_pleio]

  res = get_cause_res(df, pruned, params)

  # get r2
  path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})][[1]]
  path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})][[1]]
  if (data1 != data1_prior){df1_raw = read.table(paste0(path_in_1, 'clean/', data1, '.txt.gz'), header=1, sep = '\t')}
  if (data2 != data2_prior){df2_raw = read.table(paste0(path_in_2, 'clean/', data2, '.txt.gz'), header=1, sep = '\t')}
  out = get_harmo_res(df1_raw%>%filter(SNP%in%pruned), df2_raw, '1e-3')
  df = out[[1]]
  df$r.exposure = get_r_from_lor(df$beta.exposure, df$eaf.exposure, get_r2_para(data1)$ncase, get_r2_para(data1)$nctrl, get_r2_para(data1)$prev) 
  df$r.outcome = get_r_from_lor(df$beta.outcome, df$eaf.outcome, get_r2_para(data2)$ncase, get_r2_para(data2)$nctrl, get_r2_para(data2)$prev) 
  r2 = directionality_test(df)[,5:8]
  stat = c(res$stat, r2)
  stat[['nsnp_x']] = sum(pruned%in%bim[,2])
  res$stat = unlist(stat)

  save(res, file=file_out)
  rm(df); rm(pruned); rm(params)
}


# collect.r
library(cause)
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_t2dcat_diff_info('all')

out = c()
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  file_in = paste0(path, 'result/cause/cause_out/t2d_cat/', data1, '&', data2, '.rdata')
  load(file_in)
  print(data1); print(data2)
  print(summary(res))
  out = c(out, data1, data2, res$stat, res$force)
}

collect = data.frame(matrix(unlist(out), ncol = length(c(data1, data2, res$stat, res$r2, res$force)), byrow=T))
colnames(collect) = c('data1', 'data2', names(res$stat), 'force')

write.csv(collect, paste0(path_list[['bbj']], 'result/collect/t2d_cat/cause.csv') , row.names=F)

#=====================================================================================
# collect mr statistics
#=====================================================================================
## logit to liability factor
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_t2dcat_diff_info('all')

## gsmr
gsmr = read.csv(paste0(path, 'result/collect/t2d_cat/gsmr.csv'))
gsmr = gsmr%>%rename(nsnp=gsmr_nsnp, b=gsmr_b, se=gsmr_se, p=gsmr_pval, r2_data1=snp_r2.exposure, r2_data2=snp_r2.outcome)%>%mutate(method="GSMR")%>%select(method, data1, data2, nsnp, nsnp_x, b, se, p, r2_data1, r2_data2, pthres)

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

# test stat is updated in supp tab, other is not updated 


## note
## i drop para of sharing and epld test in supptab, bc some epld from cat-t2d is sig. 