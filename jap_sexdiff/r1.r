#=====================================================================================
# R1 cause check 
#=====================================================================================
### cause
library(cause)
path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/cause/cause_out/t2d_cat/'
trait1 = 'spracklen_nature2020_t2d'
trait2 = 'Case_control_37'
sex = 'm'
load(paste0(path, trait1, '_', sex, '&', trait2, '_', sex, '.rdata'))
summary(res)

#=====================================================================================
# R1 check sampleoverlap  
#=====================================================================================
library(MRlap)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/mrlap.r') # modified a function to aviod error
get_t2dcat_diff_info('all')
hm3_file = normalizePath('/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/w_hm3.snplist')

for (i in 1:nrow(pair)){
    print(i)
    data1 = pair[i, 'data1'];data2 = pair[i, 'data2']
    pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')
    sex = ifelse(grepl('_m', data1), 'm', 'f')
    ld_file = ifelse(grepl('Case|sprack', data1), '/home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores', 
      '/home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr')
    file_in = paste0(path_list[[1]], 'result/clump/for_clump/', pthres, '/', data1, '&', data2, '.rdata')
    file_out = paste0(path_list[[1]], 'result/mrlap/t2d_cat/', data1, '&', data2, '.rdata')
    file_clump = paste0(path_list[[1]], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')
    if (file.exists(file_out)){next()}
    path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})][[1]]
    path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})][[1]]
    df1_raw = read.table(paste0(path_in_1, 'clean/', data1, '.txt.gz'), header=1, sep = '\t')
    df2_raw = read.table(paste0(path_in_2, 'clean/', data2, '.txt.gz'), header=1, sep = '\t')
    ## ldsc
    exposure_name = data1; outcome_name = data2
    path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})][[1]]
    path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})][[1]]
    exposure = df1_raw%>%rename(rs=SNP, a1=A1, a2=A2, beta=BETA, se=SE); attributes(exposure)$GName =  deparse(substitute(exposure))
    outcome = df2_raw%>%rename(rs=SNP, a1=A1, a2=A2, beta=BETA, se=SE); attributes(outcome)$GName =  deparse(substitute(outcome))
    exposure_data = MRlap:::tidy_inputGWAS(exposure, T); outcome_data = MRlap:::tidy_inputGWAS(outcome, T)
    LDSC_results = MRlap:::run_LDSC(exposure_data, exposure_name, outcome_data, outcome_name, ld_file, hm3_file, F, T)
    
    ## correction
    snp_pleio = get_snp_pleio(data2, path_list[[1]], 'mr', ld[[ifelse(grepl('sprack|Case', data1), 'eas', 'eur')]][[sex]]) # remove for bi-test
    snp = read.table(file_clump, header=1)[,3]; snp = snp[!snp%in%snp_pleio]
    load(file_in); df = df%>%filter(SNP%in%snp)

    res_MR_TwoSampleMR =  TwoSampleMR::mr_ivw(df$beta.exposure, df$beta.outcome, df$se.exposure, df$se.outcome)
    MR_results = list("alpha_obs" = res_MR_TwoSampleMR$b,"alpha_obs_se" = res_MR_TwoSampleMR$se,
                "n_exp" = mean(df$samplesize.exposure), "n_out" = mean(df$samplesize.outcome),
                "IVs" = df%>%select(beta.exposure, se.exposure)%>%rename(std_beta.exp=beta.exposure, std_SE.exp=se.exposure), "IVs_rs" = df$SNP)
    res = get_correction(
        MR_results$IVs, LDSC_results$lambda, LDSC_results$lambda_se, LDSC_results$h2_LDSC, LDSC_results$h2_LDSC_se, 
        MR_results$alpha_obs, MR_results$alpha_obs_se, MR_results$n_exp, MR_results$n_out, as.numeric(pthres), T) # the function modified by zhanghy
    save(res, file=file_out)
}

## exchange sex iv | drop 
library(TwoSampleMR)
library(MRlap)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/mrlap.r') # modified a function to aviod error
get_t2dcat_diff_info('all')
hm3_file = normalizePath('/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/w_hm3.snplist')

for (i in 8){
    print(i)
    data1 = pair[i, 'data1'];data2 = pair[i, 'data2']
    pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')
    sex = ifelse(grepl('_m', data1), 'f', 'm') # exchange
    ld_file = ifelse(grepl('Case|sprack', data1), '/home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores', 
      '/home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr')
    data1_ex = ifelse(sex=='m', gsub('_f', '_m', data1), gsub('_m', '_f', data1))
    data2_ex = ifelse(sex=='m', gsub('_f', '_m', data2), gsub('_m', '_f', data2))

    file_out = paste0(path_list[[1]], 'result/mrlap/t2d_cat/same_iv/', data1, '&', data2, '.rdata')
    file_clump = paste0(path_list[[1]], 'result/clump/clump_out/', pthres, '/', data1_ex, '&', data2_ex, '.clumped')
    if (file.exists(file_out)){next()}
    path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})][[1]]
    path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})][[1]]
    df1_raw = read.table(paste0(path_in_1, 'clean/', data1, '.txt.gz'), header=1, sep = '\t')
    df2_raw = read.table(paste0(path_in_2, 'clean/', data2, '.txt.gz'), header=1, sep = '\t')
    ## harmo df
    out = get_harmo_res(df1_raw%>%filter(SNP%in%snp), df2_raw%>%filter(SNP%in%snp), '1')
    df = out[[1]]

    ## ldsc
    exposure_name = data1; outcome_name = data2
    path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})][[1]]
    path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})][[1]]
    exposure = df1_raw%>%rename(rs=SNP, a1=A1, a2=A2, beta=BETA, se=SE); attributes(exposure)$GName =  deparse(substitute(exposure))
    outcome = df2_raw%>%rename(rs=SNP, a1=A1, a2=A2, beta=BETA, se=SE); attributes(outcome)$GName =  deparse(substitute(outcome))
    exposure_data = MRlap:::tidy_inputGWAS(exposure, T); outcome_data = MRlap:::tidy_inputGWAS(outcome, T)
    LDSC_results = MRlap:::run_LDSC(exposure_data, exposure_name, outcome_data, outcome_name, ld_file, hm3_file, F, T)

    ## correction
    snp_pleio = get_snp_pleio(data2, path_list[[1]], 'mr', ld[[ifelse(grepl('sprack|Case', data1), 'eas', 'eur')]][[sex]]) # remove for bi-test
    snp = read.table(file_clump, header=1)[,3]; snp = snp[!snp%in%snp_pleio]


    res_MR_TwoSampleMR =  TwoSampleMR::mr_ivw(df$beta.exposure, df$beta.outcome, df$se.exposure, df$se.outcome)
    MR_results = list("alpha_obs" = res_MR_TwoSampleMR$b,"alpha_obs_se" = res_MR_TwoSampleMR$se,
                "n_exp" = mean(df$samplesize.exposure), "n_out" = mean(df$samplesize.outcome),
                "IVs" = df%>%select(beta.exposure, se.exposure)%>%rename(std_beta.exp=beta.exposure, std_SE.exp=se.exposure), "IVs_rs" = df$SNP)
    res = get_correction(
        MR_results$IVs, LDSC_results$lambda, LDSC_results$lambda_se, LDSC_results$h2_LDSC, LDSC_results$h2_LDSC_se, 
        MR_results$alpha_obs, MR_results$alpha_obs_se, MR_results$n_exp, MR_results$n_out, as.numeric(pthres), T) # the function modified by zhanghy
    save(res, file=file_out)
}


## collect
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('all')
out = c()
for (i in 1:nrow(pair)){
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  file_in = paste0(path_list[['bbj']], 'result/mrlap/t2d_cat/', data1, '&', data2, '.rdata')
  load(file_in)
  out = c(out, data1, data2, res$alpha_corrected, res$alpha_corrected_se)
}
collect = data.frame(matrix(out, ncol = 4, byrow=T))%>%rename(data1=X1, data2=X2, b=X3, se=X4)%>%
  mutate(p=2*stats::pnorm(-abs(as.numeric(b)/as.numeric(se))))

write.csv(collect, paste0(path_list[['bbj']], 'result/collect/t2d_cat/mrlap.csv') , row.names=F)

#=====================================================================================
# MR | same iv
#=====================================================================================
### eas
library(TwoSampleMR)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('eas', F)

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1'];data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')
  sex = ifelse(grepl('_m', data1), 'm', 'f'); sex_ex = ifelse(grepl('_m', data1), 'f', 'm')
  data1_ex = ifelse(sex=='m', gsub('_m', '_f', data1), gsub('_f', '_m', data1))
  data2_ex = ifelse(sex=='m', gsub('_m', '_f', data2), gsub('_f', '_m', data2))

  file_out = paste0(path_list[['bbj']], 'result/mr/t2d_cat/same_iv/', data1, '&', data2, '.rdata')
  file_clump1 = paste0(path_list[['bbj']], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')
  file_clump2 = paste0(path_list[['bbj']], 'result/clump/clump_out/', pthres, '/', data1_ex, '&', data2_ex, '.clumped')
  if (file.exists(file_out)){next()}

  snp_pleio1 = get_snp_pleio(data2, data1, path_list[['bbj']], 'mr', ld[[sex]]) # remove for bi-test
  snp1 = read.table(file_clump1, header=1)[,3]; snp1 = snp1[!snp1%in%snp_pleio1]
  snp_pleio2 = get_snp_pleio(data2_ex, data1_ex, path_list[['bbj']], 'mr', ld[[sex_ex]]) # remove for bi-test
  snp2 = read.table(file_clump2, header=1)[,3]; snp2 = snp2[!snp2%in%snp_pleio2]
  snp = unique(c(snp1, snp2))

  path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})][[1]]
  path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})][[1]]
  df1_raw = read.table(paste0(path_in_1, 'clean/', data1, '.txt.gz'), header=1, sep = '\t')%>%filter(SNP%in%snp)
  df2_raw = read.table(paste0(path_in_2, 'clean/', data2, '.txt.gz'), header=1, sep = '\t')%>%filter(SNP%in%snp)
  df = get_harmo_res(df1_raw, df2_raw, '1')[[1]]
  res = get_mr_res(df, snp); res = unlist(res)
  save(res, file=file_out)
}

### eur
library(TwoSampleMR)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('eur', F)

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1'];data2 = pair[i, 'data2']
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8')
  sex = ifelse(grepl('_m', data1), 'm', 'f'); sex_ex = ifelse(grepl('_m', data1), 'f', 'm')
  data1_ex = ifelse(sex=='m', gsub('_m', '_f', data1), gsub('_f', '_m', data1))
  data2_ex = ifelse(sex=='m', gsub('_m', '_f', data2), gsub('_f', '_m', data2))

  file_out = paste0(path_list[['bbj']], 'result/mr/t2d_cat/same_iv/', data1, '&', data2, '.rdata')
  file_clump1 = paste0(path_list[['ukb']], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')
  file_clump2 = paste0(path_list[['ukb']], 'result/clump/clump_out/', pthres, '/', data1_ex, '&', data2_ex, '.clumped')
  if (file.exists(file_out)){next()}

  snp_pleio1 = get_snp_pleio(data2, data1, path_list[['bbj']], 'mr', ld[[sex]]) # remove for bi-test
  snp1 = read.table(file_clump1, header=1)[,3]; snp1 = snp1[!snp1%in%snp_pleio1]
  snp_pleio2 = get_snp_pleio(data2_ex, data1_ex, path_list[['bbj']], 'mr', ld[[sex_ex]]) # remove for bi-test
  snp2 = read.table(file_clump2, header=1)[,3]; snp2 = snp2[!snp2%in%snp_pleio2]
  snp = unique(c(snp1, snp2))
    
  df1_raw = read.table(paste0(path_list[['ukb']], 'clean/', data1, '.txt.gz'), header=1, sep = '\t')%>%filter(SNP%in%snp)
  df2_raw = read.table(paste0(path_list[['ukb']], 'clean/', data2, '.txt.gz'), header=1, sep = '\t')%>%filter(SNP%in%snp)
  df = get_harmo_res(df1_raw, df2_raw, '1')[[1]]
  res = get_mr_res(df, snp); res = unlist(res)
  save(res, file=file_out)
}


### collect.r
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('all')

out = c()
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1'];data2 = pair[i, 'data2']
  file_out = paste0(path_list[['bbj']], 'result/mr/t2d_cat/same_iv/', data1, '&', data2, '.rdata')
  load(file_out)
  out = c(out, res)
}

collect = data.frame(matrix(out, ncol = length(res), byrow=T))
colnames(collect) = names(res)
collect = collect[,colnames(collect)[!grepl('Simple', colnames(collect))]]

write.csv(collect, paste0(path_list[['bbj']], 'result/collect/t2d_cat/mr_sameiv.csv') , row.names=F)
#=====================================================================================
# gsmr | same iv
#=====================================================================================
### eas
library(TwoSampleMR)
library(gsmr)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_gsmr_para('gsmr', 'eas') # get_gsmr_para setting
get_t2dcat_diff_info('eas', F)
gwas_thresh = 1
# r2 matrix
mat = list()
for (sex in c('m', 'f')){
  mat[[sex]] = read.table(paste0(path_list[[1]], 'result/gsmr/for_gsmr/', 'r2_mat_t2d_cat_', sex, '.ld'))
  colnames(mat[[sex]]) = rownames(mat[[sex]]) = read.table(paste0(path_list[[1]], 'result/gsmr/for_gsmr/', 'r2_mat_t2d_cat_', sex, '.snplist'))[,1]
}

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  sex = ifelse(grepl('_m', data1), 'm', 'f'); sex_ex = ifelse(grepl('_m', data1), 'f', 'm')
  data1_ex = ifelse(sex=='m', gsub('_m', '_f', data1), gsub('_f', '_m', data1))
  data2_ex = ifelse(sex=='m', gsub('_m', '_f', data2), gsub('_f', '_m', data2))
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8'); pthres2 = ifelse(data1_ex %in% datas_proxy, '1e-5', '5e-8')

  file_out = paste0(path_list[[1]], 'result/gsmr/gsmr_out/t2d_cat/same_iv/', data1, '&', data2, '.rdata')
  file_clump1 = paste0(path_list[['bbj']], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')
  file_clump2 = paste0(path_list[['bbj']], 'result/clump/clump_out/', pthres2, '/', data1_ex, '&', data2_ex, '.clumped')
  if (file.exists(file_out)){next()}

  snp_pleio1 = get_snp_pleio(data2, data1, path_list[['bbj']], 'mr', ld[[sex]]) # remove for bi-test
  snp1 = read.table(file_clump1, header=1)[,3]; snp1 = snp1[!snp1%in%snp_pleio1]
  snp_pleio2 = get_snp_pleio(data2_ex, data1_ex, path_list[['bbj']], 'mr', ld[[sex_ex]]) # remove for bi-test
  snp2 = read.table(file_clump2, header=1)[,3]; snp2 = snp2[!snp2%in%snp_pleio2]
  snp = unique(c(snp1, snp2))

  path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})][[1]]
  path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})][[1]]
  df1_raw = read.table(paste0(path_in_1, 'clean/', data1, '.txt.gz'), header=1, sep = '\t')%>%filter(SNP%in%snp)
  df2_raw = read.table(paste0(path_in_2, 'clean/', data2, '.txt.gz'), header=1, sep = '\t')%>%filter(SNP%in%snp)
  df = get_harmo_res(df1_raw, df2_raw, '1')[[1]]

  res = get_gsmr_res(df, snp, mat[[sex]]); res = unlist(res)
  save(res, file=file_out)
}

### eur
library(TwoSampleMR)
library(gsmr)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_gsmr_para('gsmr', 'eur') # get_gsmr_para setting
get_t2dcat_diff_info('eur', F)
gwas_thresh = 1
# r2 matrix
mat = list()
for (sex in c('m', 'f')){
  mat[[sex]] = read.table(paste0(path_list[['ukb']], 'result/gsmr/for_gsmr/', 'r2_mat_t2d_cat_', sex, '.ld'))
  colnames(mat[[sex]]) = rownames(mat[[sex]]) = read.table(paste0(path_list[['ukb']], 'result/gsmr/for_gsmr/', 'r2_mat_t2d_cat_', sex, '.snplist'))[,1]
}
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  data1_ex = ifelse(sex=='m', gsub('_m', '_f', data1), gsub('_f', '_m', data1))
  data2_ex = ifelse(sex=='m', gsub('_m', '_f', data2), gsub('_f', '_m', data2))
  pthres = ifelse(data1 %in% datas_proxy, '1e-5', '5e-8'); pthres2 = ifelse(data1_ex %in% datas_proxy, '1e-5', '5e-8')
  sex = ifelse(grepl('_m', data1), 'm', 'f'); sex_ex = ifelse(grepl('_m', data1), 'f', 'm')

  file_out = paste0(path_list[['bbj']], 'result/gsmr/gsmr_out/t2d_cat/same_iv/', data1, '&', data2, '.rdata')
  file_clump1 = paste0(path_list[['ukb']], 'result/clump/clump_out/', pthres, '/', data1, '&', data2, '.clumped')
  file_clump2 = paste0(path_list[['ukb']], 'result/clump/clump_out/', pthres2, '/', data1_ex, '&', data2_ex, '.clumped')
  if (file.exists(file_out)){next}

  snp_pleio1 = get_snp_pleio(data2, data1, path_list[['ukb']], 'mr', ld[[sex]]) # remove for bi-test
  snp1 = read.table(file_clump1, header=1)[,3]; snp1 = snp1[!snp1%in%snp_pleio1]
  snp_pleio2 = get_snp_pleio(data2_ex, data1_ex, path_list[['ukb']], 'mr', ld[[sex_ex]]) # remove for bi-test
  snp2 = read.table(file_clump2, header=1)[,3]; snp2 = snp2[!snp2%in%snp_pleio2]
  snp = unique(c(snp1, snp2))

  df1_raw = read.table(paste0(path_list[['ukb']], 'clean/', data1, '.txt.gz'), header=1, sep = '\t')%>%filter(SNP%in%snp)
  df2_raw = read.table(paste0(path_list[['ukb']], 'clean/', data2, '.txt.gz'), header=1, sep = '\t')%>%filter(SNP%in%snp)
  df = get_harmo_res(df1_raw, df2_raw, '1')[[1]]
  res = get_gsmr_res(df, snp, mat[[sex]]); res = unlist(res)
  save(res, file=file_out)
}

# collect.r
get_t2dcat_diff_info('all')
out = c()
for (i in 1:nrow(pair)){
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  file_in = paste0(path_list[['bbj']], 'result/gsmr/gsmr_out/t2d_cat/same_iv/', data1, '&', data2, '.rdata')
  load(file_in)
  out = c(out, res)
}
collect = data.frame(matrix(out, ncol = length(res), byrow=T))
colnames(collect) = names(res)

write.csv(collect, paste0(path_list[['bbj']], 'result/collect/t2d_cat/gsmr_sameiv.csv') , row.names=F)
#=====================================================================================
# cause | same iv
#=====================================================================================
### eas
library(TwoSampleMR)
library(cause)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('eas', F)

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  sex = ifelse(grepl('_m', data1), 'm', 'f'); sex_ex = ifelse(grepl('_m', data1), 'f', 'm')
  data1_ex = ifelse(sex=='m', gsub('_m', '_f', data1), gsub('_f', '_m', data1))
  data2_ex = ifelse(sex=='m', gsub('_m', '_f', data2), gsub('_f', '_m', data2))

  file_in1 = paste0(path_list[[1]], 'result/cause/for_cause/prune/t2d_cat/', data1, '&', data2, '.rdata')
  file_in2 = paste0(path_list[[1]], 'result/cause/for_cause/prune/t2d_cat/', data1_ex, '&', data2_ex, '.rdata')
  file_in3 = paste0(path_list[[1]], 'result/cause/for_cause/param/t2d_cat/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[[1]], 'result/cause/cause_out/t2d_cat/same_iv/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){next()}
  load(file_in1) 
  snp_pleio1 = get_snp_pleio(data2, data1, path_list[[1]], 'cause', ld[[sex]]) # remove for bi-test
  pruned1 = res$pruned; pruned1 = pruned1[!pruned1%in%snp_pleio1]
  load(file_in2)
  snp_pleio2 = get_snp_pleio(data2_ex, data1_ex, path_list[[1]], 'cause', ld[[sex_ex]]) # remove for bi-test
  pruned2 = res$pruned; pruned2 = pruned2[!pruned2%in%snp_pleio2]
  pruned = unique(c(pruned1, pruned2))

  load(file_in3)

  path_in_1 = path_list[sapply(names(path_list), function(x){grepl(x, data1)})][[1]]
  path_in_2 = path_list[sapply(names(path_list), function(x){grepl(x, data2)})][[1]]
  df1_raw = read.table(paste0(path_in_1, 'clean/', data1, '.txt.gz'), header=1, sep = '\t')%>%filter(SNP%in%pruned)
  df2_raw = read.table(paste0(path_in_2, 'clean/', data2, '.txt.gz'), header=1, sep = '\t')%>%filter(SNP%in%pruned)
  df <- gwas_merge(df1_raw, df2_raw, snp_name_cols = c("SNP", "SNP"), beta_hat_cols = c("BETA", "BETA"), 
    se_cols = c("SE", "SE"), A1_cols = c("A1", "A1"), A2_cols = c("A2", "A2"))

  res = get_cause_res(df, df$snp, params)

  save(res, file=file_out)
  rm(df); rm(pruned); rm(params)
}

## eur
library(TwoSampleMR)
library(cause)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('eur', F)

for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  sex = ifelse(grepl('_m', data1), 'm', 'f'); sex_ex = ifelse(grepl('_m', data1), 'f', 'm')
  data1_ex = ifelse(sex=='m', gsub('_m', '_f', data1), gsub('_f', '_m', data1))
  data2_ex = ifelse(sex=='m', gsub('_m', '_f', data2), gsub('_f', '_m', data2))

  file_in1 = paste0(path_list[['bbj']], 'result/cause/for_cause/prune/t2d_cat/', data1, '&', data2, '.rdata')
  file_in2 = paste0(path_list[['bbj']], 'result/cause/for_cause/prune/t2d_cat/', data1_ex, '&', data2_ex, '.rdata')
  file_in3 = paste0(path_list[['bbj']], 'result/cause/for_cause/param/t2d_cat/', data1, '&', data2, '.rdata')
  file_out = paste0(path_list[['bbj']], 'result/cause/cause_out/t2d_cat/same_iv/', data1, '&', data2, '.rdata')
  
  if (file.exists(file_out)){next()}
  load(file_in1) 
  snp_pleio1 = get_snp_pleio(data2, data1, path_list[[1]], 'cause', ld[[sex]]) # remove for bi-test
  pruned1 = res$pruned; pruned1 = pruned1[!pruned1%in%snp_pleio1]
  load(file_in2)
  snp_pleio2 = get_snp_pleio(data2_ex, data1_ex, path_list[[1]], 'cause', ld[[sex_ex]]) # remove for bi-test
  pruned2 = res$pruned; pruned2 = pruned2[!pruned2%in%snp_pleio2]
  pruned = unique(c(pruned1, pruned2))

  load(file_in3)

  df1_raw = read.table(paste0(path_list[['ukb']], 'clean/', data1, '.txt.gz'), header=1, sep = '\t')%>%filter(SNP%in%pruned)
  df2_raw = read.table(paste0(path_list[['ukb']], 'clean/', data2, '.txt.gz'), header=1, sep = '\t')%>%filter(SNP%in%pruned)
  df <- gwas_merge(df1_raw, df2_raw, snp_name_cols = c("SNP", "SNP"), beta_hat_cols = c("BETA", "BETA"), 
    se_cols = c("SE", "SE"), A1_cols = c("A1", "A1"), A2_cols = c("A2", "A2"))

  res = get_cause_res(df, pruned, params)

  save(res, file=file_out)
  rm(df); rm(pruned); rm(params)
}


# collect.r
library(cause)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('all')

out = c()
for (i in 1:nrow(pair)){
  print(i)
  data1 = pair[i, 'data1']; data2 = pair[i, 'data2']
  file_in = paste0(path, 'result/cause/cause_out/t2d_cat/same_iv/', data1, '&', data2, '.rdata')
  load(file_in)
  print(data1); print(data2)
  print(summary(res))
  out = c(out, data1, data2, res$stat, res$force)
}

collect = data.frame(matrix(unlist(out), ncol = length(c(data1, data2, res$stat, res$r2, res$force)), byrow=T))
colnames(collect) = c('data1', 'data2', names(res$stat), 'force')

write.csv(collect, paste0(path_list[['bbj']], 'result/collect/t2d_cat/cause_sameiv.csv') , row.names=F)

#=====================================================================================
# collect mr statistics
#=====================================================================================
## logit to liability factor
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('all')

## gsmr
gsmr = read.csv(paste0(path, 'result/collect/t2d_cat/gsmr_sameiv.csv'))
gsmr = gsmr%>%rename(nsnp=gsmr_nsnp, b=gsmr_b, se=gsmr_se, p=gsmr_pval)%>%mutate(method="GSMR")%>%select(method, data1, data2, nsnp, b, se, p, pthres)

## mr
df = read.csv(paste0(path, 'result/collect/t2d_cat/mr_sameiv.csv'))
mr = data.frame()
for (method in c('MR.Egger', 'Inverse.variance.weighted', 'Weighted.median', 'Weighted.mode')){
  sub = df[,c(names(df)[grepl(method, names(df))], 'data1', 'data2', 'pthres')]%>%mutate(t=method)
  names(sub) = c('nsnp', 'b', 'se', 'p','data1', 'data2', 'pthres', 'method')
  mr = rbind(mr, sub)
}

## cause
cause = read.csv(paste0(path, 'result/collect/t2d_cat/cause_sameiv.csv'))
cause = cause%>%rename(nsnp=cause_nsnp, b=cause_causal_gamma_b, se=cause_causal_gamma_se, p=cause_causal_gamma_pval, r2_data1=snp_r2.exposure, r2_data2=snp_r2.outcome)%>%
  mutate(method="CAUSE", pthres=1e-3)%>%select(method, data1, data2, nsnp, nsnp_x, b, se, p, r2_data1, r2_data2, pthres)

## combine
res = data.frame(rbind(cause, gsmr, mr))
res = res%>%mutate(b_l=b-1.96*se, b_u=b+1.96*se)
## logit to liability
res[res$data1 == 'E11_m', c('b', 'b_l', 'b_u')] = res[res$data1 == 'E11_m', c('b', 'b_l', 'b_u')]*get_lia_factor(prev_t2d_m_eur, prev_cataract_m_eur)
res[res$data1 == 'E11_f', c('b', 'b_l', 'b_u')] = res[res$data1 == 'E11_f', c('b', 'b_l', 'b_u')]*get_lia_factor(prev_t2d_f_eur, prev_cataract_f_eur)
res[res$data2 == 'E11_m', c('b', 'b_l', 'b_u')] = res[res$data2 == 'E11_m', c('b', 'b_l', 'b_u')]*get_lia_factor(prev_cataract_m_eur, prev_t2d_m_eur)
res[res$data2 == 'E11_f', c('b', 'b_l', 'b_u')] = res[res$data2 == 'E11_f', c('b', 'b_l', 'b_u')]*get_lia_factor(prev_cataract_f_eur, prev_t2d_f_eur)

res[res$data1 == 'Case_control_37_m', c('b', 'b_l', 'b_u')] = res[res$data1 == 'Case_control_37_m', c('b', 'b_l', 'b_u')]*get_lia_factor(prev_cataract_m_eas, prev_t2d_m_eas)
res[res$data1 == 'Case_control_37_f', c('b', 'b_l', 'b_u')] =  res[res$data1 == 'Case_control_37_f', c('b', 'b_l', 'b_u')]*get_lia_factor(prev_cataract_f_eas, prev_t2d_f_eas)
res[res$data2 == 'Case_control_37_m', c('b', 'b_l', 'b_u')] =  res[res$data2 == 'Case_control_37_m', c('b', 'b_l', 'b_u')]*get_lia_factor(prev_t2d_m_eas, prev_cataract_m_eas)
res[res$data2 == 'Case_control_37_f', c('b', 'b_l', 'b_u')] =  res[res$data2 == 'Case_control_37_f', c('b', 'b_l', 'b_u')]*get_lia_factor(prev_t2d_f_eas, prev_cataract_f_eas)

res[,c('or', 'or_l', 'or_u')] = exp(res[,c('b', 'b_l', 'b_u')])
res[res=='MR.Egger'] = 'MR-Egger'; res[res=='Inverse.variance.weighted'] = 'IVW'; res[res=='Weighted.median'] = 'Weighted median'; res[res=='Weighted.mode'] = 'Weighted mode'
write.csv(res, '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/collect/t2d_cat/collect.csv', row.names=F)


t2d_m_eur = res%>%filter(data1=='E11_m')%>%pull('b')
t2d_f_eur = res%>%filter(data1=='E11_f')%>%pull('b')
cat_m_eur = res%>%filter(data1=='cataract_m')%>%pull('b')
cat_f_eur = res%>%filter(data1=='cataract_f')%>%pull('b')

t2d_m_agen = res%>%filter(data1=='spracklen_nature2020_t2d_m'&data2=='Case_control_37_m')%>%pull('b')
t2d_f_agen = res%>%filter(data1=='spracklen_nature2020_t2d_f'&data2=='Case_control_37_f')%>%pull('b')
cat_m_agen = res%>%filter(data1=='Case_control_37_m'&data2=='spracklen_nature2020_t2d_m')%>%pull('b')
cat_f_agen = res%>%filter(data1=='Case_control_37_f'&data2=='spracklen_nature2020_t2d_f')%>%pull('b')

t2d_m_bbj = res%>%filter(data1=='Case_control_96_m'&data2=='Case_control_37_m')%>%pull('b')
t2d_f_bbj = res%>%filter(data1=='Case_control_96_f'&data2=='Case_control_37_f')%>%pull('b')
cat_m_bbj = res%>%filter(data1=='Case_control_37_m'&data2=='Case_control_96_m')%>%pull('b')
cat_f_bbj = res%>%filter(data1=='Case_control_37_f'&data2=='Case_control_96_f')%>%pull('b')

t.test(t2d_m_agen, t2d_f_agen, paired=T)
t.test(cat_m_agen, cat_f_agen, paired=T)

t.test(t2d_m_bbj, t2d_f_bbj, paired=T)
t.test(cat_m_bbj, cat_f_bbj, paired=T)

t.test(t2d_m_eur, t2d_f_eur, paired=T)
t.test(cat_m_eur, cat_f_eur, paired=T)

# test stat is updated in supp tab, other is not updated 


## note
## i drop para of sharing and epld test in supptab, bc some epld from cat-t2d is sig. 





