# t2d_cat in-house ukb 
#=====================================================================================
# clean bolt gwas
#=====================================================================================
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/gwas/summary/ukb/eur/'

bim = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt.bim')[,c(2,1,4)]
bim = bim%>%rename(SNP=V2, CHR=V1, POS=V4)%>%filter(!(CHR==6&POS<=33448354&POS>=28477797))

nsample = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/data/brt_pheno_for_bolt.txt', header=1)
nsample1 = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/data/brt_m_0.8_pheno_for_bolt.txt', header=1)
nsample2 = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/data/brt_f_0.8_pheno_for_bolt.txt', header=1)
nsample_0.8 = rbind(nsample1, nsample2)

datas = list.files(paste0(path, 'raw/'))
datas = datas[!grepl('log|old', datas)]
datas = datas[grepl('0.8', datas)]

for (data in datas[3:4]){
  print(data)
  file_out = paste0(path, 'clean/', gsub('.gz', '', data))
  if (file.exists(file_out)){next}

  df = read.table(paste0(path, 'raw/', data), header=1, fill=T)
  df = df[, c('CHR', 'BP', 'ALLELE1', 'ALLELE0', 'A1FREQ', 'BETA', 'SE', 'P_BOLT_LMM_INF')]

  if (grepl('sitting', data)){
    if (grepl('_m|_f', data)){
      trait = strsplit(data, ifelse(grepl('_m', data), '_m', '_f'))[[1]][1]
      N = nrow(nsample%>%filter(sex==ifelse(grepl('_m', data), 1, 0))) # problematic, some samples with missing sitting height
    }
    else {
      trait = gsub('.txt.gz', '', data)
      N = nrow(nsample)
    }
  } else if (grepl('0.8', data)){
    trait = ifelse(grepl('cat', data), 'cataract', 'E11')
    frq = table(nsample_0.8%>%filter(sex==ifelse(grepl('_m', data), 1, 0))%>%pull(trait))
  }
  else {
    if (grepl('_m|_f', data)){
      trait = ifelse(grepl('cat', data), 'cataract', 'E11')
      frq = table(nsample%>%filter(sex==ifelse(grepl('_m', data), 1, 0))%>%pull(trait))
    } else {
      trait = gsub('.txt.gz', '', data)
      frq = table(nsample%>%pull(trait))
    }
  }
  ncase = frq[['1']]; nctrl = frq[['0']]
  N = ncase+nctrl
  u = ncase/(ncase+nctrl)

  df = df%>%rename(POS=BP, A1=ALLELE1, A2=ALLELE0, FRQ=A1FREQ, P=P_BOLT_LMM_INF)%>%filter(FRQ>0.01&FRQ<0.99)%>%merge(bim, by = c('CHR', 'POS'))%>%mutate(N=N)%>%select(SNP, CHR, POS, FRQ, A1, A2, BETA, SE, P, N)%>%na.omit()
  df = df%>%mutate(id=1:nrow(df))%>%arrange(P)%>%filter(!duplicated(SNP))%>%arrange(id)%>%select(-id)

  # bolt need to convert beta to logit scale bc it use linear regression
  # log(or) = beta/(u(1-u)), u is case fraction
  if (!grepl('sitting', data)){
    df = df%>%mutate(N=N, BETA=BETA/u, SE=SE/u)
  }
  write.table(df, file_out, row.names = F, sep = '\t', quote = F)
  rm(list = c('ncase', 'nctrl'))
}

t1 = table(nsample%>%filter(sex==1)%>%pull(E11))
t2 = table(nsample%>%filter(sex==0)%>%pull(E11))
t3 = table(nsample%>%filter(sex==1)%>%pull(cataract))
t4 = table(nsample%>%filter(sex==0)%>%pull(cataract))

t = data.frame(matrix(c(t1, t2, t3, t4), ncol=2, byrow=T))
t = t%>%rename(nctrl=X1, ncase=X2)%>%mutate(prop_sam=ncase/(ncase+nctrl), trait=c('t2d', 't2d', 'cat', 'cat'), sex=c('m', 'f', 'm', 'f'))
#    nctrl ncase   prop_sam trait sex
# 1 181145 14059 0.07202209   t2d   m
# 2 221440  8754 0.03802879   t2d   f
# 3 183067 12137 0.06217598   cat   m
# 4 214481 15713 0.06825982   cat   f
#=====================================================================================
# ldsc
#=====================================================================================
## munge
path="/home/yanglab_data/user/zhanghy/gwas/summary/ukb/eur/"
file_ref="/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/w_hm3.snplist_withx"

files=`ls ${path}/clean/ | grep txt`
for file in $files
do
echo $file
out=`echo $file | sed -e "s/.txt.gz//g"`
if [ ! -f ${path}/munge/$out.sumstats.gz ]
then
/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/munge_sumstats.py \
--sumstats ${path}clean/$file \
--out ${path}/munge/$out \
--merge-alleles $file_ref
fi
done

## ldsc
path_in="/home/yanglab_data/user/zhanghy/gwas/summary/ukb/eur/munge/"

prev_t2d_m_pop=0.038; prev_t2d_f_pop=0.031
prev_cataract_m_pop=0.302; prev_cataract_f_pop=0.4069
prev_t2d_m_sam=0.07202209; prev_t2d_f_sam=0.03802879
prev_cataract_m_sam=0.06217598; prev_cataract_f_sam=0.06825982
prev_t2d_08_m_sam=0.07186722; prev_t2d_08_f_sam=0.03803318
prev_cataract_08_m_sam=0.06210178; prev_cataract_08_f_sam=0.06848036

sex='f'
data1=E11_${sex}
data2=cataract_${sex}
data1=E11_0.8_${sex}
data2=cataract-adjE11_0.8_${sex}

sam_prev1=`eval echo '$'"prev_t2d_${sex}_sam"`; sam_prev2=`eval echo '$'"prev_cataract_${sex}_sam"`
pop_prev1=`eval echo '$'"prev_t2d_${sex}_pop"`; pop_prev2=`eval echo '$'"prev_cataract_${sex}_pop"`
sam_prev1=`eval echo '$'"prev_t2d_08_${sex}_sam"`;sam_prev2=`eval echo '$'"prev_cataract_08_${sex}_sam"`

path_ld="/home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr/${sex}/"
file_out="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ldsc/t2d_cat/${data1}&${data2}"

/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/ldsc.py \
--rg  $path_in$data1.sumstats.gz,$path_in$data2.sumstats.gz \
--ref-ld-chr $path_ld --w-ld-chr $path_ld \
--samp-prev $sam_prev1,$sam_prev2 \
--pop-prev $pop_prev1,$pop_prev2 \
--out ${file_out}_lia

/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/ldsc.py \
--rg  $path_in$data1.sumstats.gz,$path_in$data2.sumstats.gz \
--ref-ld-chr $path_ld --w-ld-chr $path_ld \
--samp-prev $sam_prev1,$sam_prev2 \
--pop-prev $pop_prev1,$pop_prev2 \
--intercept-h2 1,1 \
--out ${file_out}_constrain_lia

#=====================================================================================
# clump
#=====================================================================================
## clump.r
library(TwoSampleMR)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('eur')
data1_prior = data2_prior = 'none'

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')

  file_out1 = paste0(path_list[['ukb']], 'result/clump/for_clump/', pthres, '/', data1, '&', data2, '.rdata')
  file_out2 = paste0(path_list[['ukb']], 'result/clump/for_clump/', pthres, '/', data1, '&', data2, '.txt')

  if (file.exists(file_out1)&(file.exists(file_out2)|file.exists(paste0(file_out2, '.gz')))){next}
  if (data1 != data1_prior){df1_raw = read.table(paste0(path_list[['ukb']], 'clean/', data1, '.txt.gz'), header=1, sep = '\t')}
  if (data2 != data2_prior){df2_raw = read.table(paste0(path_list[['ukb']], 'clean/', data2, '.txt.gz'), header=1, sep = '\t')}

  out = get_harmo_res(df1_raw, df2_raw, pthres)
  df = out[[1]]

  save(df, file=file_out1)
  write.table(out[[2]], file_out2, row.names=F, sep='\t', quote=F)
  data1_prior = data1; data2_prior = data2
}

### clump
path_clump="/home/yanglab_data/user/zhanghy/gwas/summary/ukb/eur/result/clump/"
datas_proxy="Case_control_37_f Case_control_37_m"

for sex in m f
do
file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/brt_${sex}"

data1="E11_0.8_${sex}"
data2="cataract-adjE11_0.8_${sex}"

if [ `echo $datas_proxy | grep -w $data1 | wc -l` == 1 ]
then pthres="1e-5"
else pthres="5e-8"
fi

file_in=${path_clump}for_clump/${pthres}/${data1}\&${data2}.txt.gz
file_out=${path_clump}clump_out/${pthres}/${data1}\&${data2}

if [ -f $file_out.clumped ]; then continue; fi

plink --allow-no-sex \
--clump $file_in --bfile $file_bfile --out $file_out \
--clump-kb 1000 --clump-p1 1.0 --clump-p2 1.0 --clump-r2 0.05 \

done
#=====================================================================================
# ld of instruments, remove for bi-dir test
#=====================================================================================
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('eur', F)

out = c()
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')
  snp1 = read.table(paste0(path_list[['ukb']], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped'), header=1)[,3]
  load(paste0(path_list[['bbj']], 'result/cause/for_cause/prune/t2d_cat/', data1, '&', data2, '.rdata'))
  snp2 = res$pruned
  out = c(out, snp1, snp2)
}
out = unique(out)

write.table(data.frame(out), paste0(path_list[['bbj']], 'result/ld/for_ld/t2d-cat_eur'), row.names = F, col.names = F, quote = F)

## r2 
for sex in m f
do
file_snp="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ld/for_ld/t2d-cat_eur"
file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/brt_${sex}"
file_out="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ld/ld_out/t2d-cat_${sex}_eur_0.6"

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
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('eur')

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1'];data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')
  sex = ifelse(grepl('_m', data1), 'm', 'f')

  file_in = paste0(path_list[['ukb']], 'result/clump/for_clump/', pthres, '/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[['bbj']], 'result/mr/t2d_cat/', data1, '&', data2, '.rdata')
  file_clump = paste0(path_list[['ukb']], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')

  if (file.exists(file_out)){next()}
  load(file_in)
  snp_pleio = get_snp_pleio(data2, data1, path_list[['ukb']], 'mr', ld[[sex]]) # remove for bi-test
  snp = read.table(file_clump, header=1)[,3]
  snp = snp[!snp%in%snp_pleio]
    
  res = get_mr_res(df, snp)
  
  df$r.exposure = get_r_from_lor(df$beta.exposure, df$eaf.exposure, get_r2_para(data1)$ncase, get_r2_para(data1)$nctrl, get_r2_para(data1)$prev) 
  df$r.outcome = get_r_from_lor(df$beta.outcome, df$eaf.outcome, get_r2_para(data2)$ncase, get_r2_para(data2)$nctrl, get_r2_para(data2)$prev) 
  r2 = directionality_test(df)[,5:8]
  
  res[['nsnp_x']] = sum(snp%in%bim[,2])
  res = c(res, r2); res = unlist(res)
  save(res, file=file_out)
}
#=====================================================================================
# gsmr
#=====================================================================================
## make_snp_list_for_extract.r 
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('eur', F)

out = c()
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1'];data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')
  file_clump = paste0(path_list[['ukb']], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')
  snp = read.table(file_clump, header=1)[,3]
  out = c(out, snp)
}
out = unique(out)
write.table(data.frame(out), paste0(path_list[['ukb']], 'result/gsmr/for_gsmr/r2_mat_t2d_cat_snp'), row.names = F, col.names = F, quote = F)

# r2  matrix
for sex in m f
do
file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/brt_${sex}"

plink --bfile $file_bfile \
--extract /home/yanglab_data/user/zhanghy/gwas/summary/ukb/eur/result/gsmr/for_gsmr/r2_mat_t2d_cat_snp \
--r2 square \
--write-snplist \
--out /home/yanglab_data/user/zhanghy/gwas/summary/ukb/eur/result/gsmr/for_gsmr/r2_mat_t2d_cat_${sex}

done

### gsmr.r
library(TwoSampleMR)
library(gsmr)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_gsmr_para('gsmr', 'eur') # get_gsmr_para setting
get_t2dcat_diff_info('eur')

# r2 matrix
mat = list()
for (sex in c('m', 'f')){
  mat[[sex]] = read.table(paste0(path_list[['ukb']], 'result/gsmr/for_gsmr/', 'r2_mat_t2d_cat_', sex, '.ld'))
  colnames(mat[[sex]]) = rownames(mat[[sex]]) = read.table(paste0(path_list[['ukb']], 'result/gsmr/for_gsmr/', 'r2_mat_t2d_cat_', sex, '.snplist'))[,1]
}

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')
  gwas_thresh = as.numeric(pthres)
  sex = ifelse(grepl('_m', data1), 'm', 'f')

  file_in = paste0(path_list[['ukb']], 'result/clump/for_clump/', pthres, '/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[['bbj']], 'result/gsmr/gsmr_out/t2d_cat/', data1, '&', data2, '.rdata')
  file_clump = paste0(path_list[['ukb']], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')
  if (file.exists(file_out)){next}
  load(file_in)
  snp_pleio = get_snp_pleio(data2, data1, path_list[['ukb']], 'mr', ld[[sex]]) # remove for bi-test
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

#=====================================================================================
# cause
#=====================================================================================
# cause_prune.r 
library(cause)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('eur', F)
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  sex = ifelse(grepl('_m', data1), 'm', 'f')
  file_out = paste0(path_list[['bbj']], '/result/cause/for_cause/prune/t2d_cat/', data1, '&', data2, '.rdata')
  if (file.exists(file_out)){next}

  df1_raw = read.table(paste0(path_list[['ukb']], 'clean/', data1, '.txt.gz'), header=1, sep = '\t')
  df2_raw = read.table(paste0(path_list[['ukb']], 'clean/', data2, '.txt.gz'), header=1, sep = '\t')

  df <- gwas_merge(df1_raw, df2_raw, snp_name_cols = c("SNP", "SNP"), beta_hat_cols = c("BETA", "BETA"), 
    se_cols = c("SE", "SE"), A1_cols = c("A1", "A1"), A2_cols = c("A2", "A2"))
  variants <- df %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
  pruned = c()

  for (chr in 1:23){
    print(paste0('chr ', chr, ' pruned snp:'))
    ld <- readRDS(paste0(path_list[['ref']], 'brt_', sex, '_', chr, '_ld.RDS'))
    snp_info <- readRDS(paste0(path_list[['ref']], 'brt_', sex, '_', chr, '_snp.RDS'))
    out <- ld_prune(variants = variants, ld = ld, total_ld_variants = snp_info$SNP, pval_cols = c("pval1"),pval_thresh = c(1e-3))
    pruned = c(pruned, out)
    print(length(pruned))
  }

  res = list()
  res$df = df
  res$pruned = pruned
  save(res, file=file_out)
}

# cause_param.r
library(cause)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('eur', F)

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  file_in = paste0(path_list[['bbj']], 'result/cause/for_cause/prune/t2d_cat/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[['bbj']], 'result/cause/for_cause/param/t2d_cat/', data1, '&', data2, '.rdata')

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
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('eur')

data1_prior = data2_prior = 'none'
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  file_in1 = paste0(path_list[['bbj']], 'result/cause/for_cause/prune/t2d_cat/', data1, '&', data2, '.rdata')
  file_in2 = paste0(path_list[['bbj']], 'result/cause/for_cause/param/t2d_cat/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[['bbj']], 'result/cause/cause_out/t2d_cat/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){next()}
  load(file_in1); load(file_in2)

  snp_pleio = get_snp_pleio(data2, data1, path_list[['bbj']], 'cause', ld[[sex]]) # remove for bi-test
  df = res$df; pruned = res$pruned; pruned = pruned[!pruned%in%snp_pleio]
  res = get_cause_res(df, pruned, params)

  # get r2
  if (data1 != data1_prior){df1_raw = read.table(paste0(path_list[['ukb']], 'clean/', data1, '.txt.gz'), header=1, sep = '\t')}
  if (data2 != data2_prior){df2_raw = read.table(paste0(path_list[['ukb']], 'clean/', data2, '.txt.gz'), header=1, sep = '\t')}
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