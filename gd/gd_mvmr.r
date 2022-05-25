library(dplyr)
h2 = read.table('/home/yanglab_data/user/zhanghy/gwas/summary/neale/para/h2.txt', header=1)
keep2 = h2%>%filter(grepl('_f', data)&p<0.05)%>%pull(data)
keep2 = gsub('_f.txt', '', keep2)
nsnp = read.table('/home/yanglab_data/user/zhanghy/gwas/summary/neale/para/nsnp_5e-8.txt', header=0)
keep1 = nsnp%>%filter(grepl('_f', V1)&nsnp[,2]>=10)%>%pull(V1)
keep1 = gsub('_f.txt', '', keep1)
keep = intersect(keep1, keep2)

df = readxl::read_xlsx("/home/yanglab_data/user/zhanghy/gwas/summary/neale/info/ukbb_v3_20180731.xlsx", sheet = 4, col_names=F) 
df = setNames(data.frame(df), c('trait', 'label'))
df = df%>%filter(trait%in%keep)%>%distinct()
trait = unique(df$trait)


files = list.files('/home/yanglab_data/user/zhanghy/gwas/summary/neale/clean/')
files[files%in%paste0(trait, '_f.txt.gz')]



## mvmr
# covar
# bmi: QTL_2_m, QTL_4_f
# smoke: cigaratter per day: QTL_131_m, QTL_133_f
# pad: Case_control_87_m, Case_control_87_f
# cad: Case_control_104_m, Case_control_104_f
# osteoporosis: Case_control_80_m, Case_control_80_f
# Rheumatoid arthritis Case_control_94_m, Case_control_94_m
#=====================================================================================
# clump
#=====================================================================================
## for_clump.r
## combine exposure for clump
source('/home/yanglab_data/user/zhanghy/project/temp/code/source.r')
library(TwoSampleMR)
options(stringsAsFactors = F)

path_list = list()
path_list[['bbj']] = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/'
path_list[['agen']] = '/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/'
covar_list = list()
covar_list[['m']] = c('QTL_2_m', 'QTL_131_m', 'Case_control_87_m', 'Case_control_104_m', 'Case_control_80_m', 'Case_control_94_m')
covar_list[['f']] = c('QTL_4_f', 'QTL_133_f', 'Case_control_87_f', 'Case_control_104_f', 'Case_control_80_f', 'Case_control_94_f')

for (sex in c('m', 'f')){
  data1 = paste0('spracklen_nature2020_t2d_', sex); data3 = paste0('Case_control_37_', sex)
  df1_raw = read.table(paste0(path_list[['agen']], 'clean/', data1, '.txt.gz'), sep='\t', header=1)
  df3_raw = read.table(paste0(path_list[['bbj']], 'clean/', data3, '.txt.gz'), sep='\t', header=1)
  for (data2 in covar_list[[sex]]){ # here data2 is covar
    df2_raw = read.table(paste0(path_list[['bbj']], 'clean/', data2, '.txt.gz'), sep='\t', header=1)
    for (pthres in c('1e-5', '5e-8')){
      out = get_mvharmo_res(df1_raw, df2_raw, df3_raw, data1, data2, data3, pthres)
      file_out1 = paste0(path_list[['bbj']], 'result/clump/for_clump/mvmr/', pthres, '/', data1, '&', data3, '_co_', data2, '.txt')
      file_out2 = paste0(path_list[['bbj']], 'result/clump/for_clump/mvmr/', pthres, '/', data1, '&', data3, '_co_', data2, '.rdata')
      df_comb = rbind(df1_raw%>%filter(SNP%in%out[[2]]), df2_raw%>%filter(SNP%in%out[[2]]))%>%arrange(P)%>%filter (!duplicated(SNP))
      write.table(df_comb, file_out1, row.names=F, sep='\t', quote=F)
      df = out[[1]]
      save(df, file=file_out2)
    }
  }
}

## clump
path_clump="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/clump/"

for sex in m f
do file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/chn_${sex}"
for pthres in 1e-5 5e-8
do for data in `ls ${path_clump}/for_clump/mvmr/${pthres} | grep txt`
do 
file_in=${path_clump}/for_clump/mvmr/${pthres}/$data
file_out=${path_clump}/clump_out/mvmr/${pthres}/${data/.txt/}

plink --allow-no-sex \
--clump $file_in --bfile $file_bfile --out $file_out \
--clump-kb 1000 --clump-p1 1.0 --clump-p2 1.0 --clump-r2 0.05 \

done; done; done

#=====================================================================================
# mvmr
#=====================================================================================
source('/home/yanglab_data/user/zhanghy/project/temp/code/source.r')
library(TwoSampleMR)
options(stringsAsFactors = F)

path_list = list()
path_list[['bbj']] = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/'
path_list[['agen']] = '/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/'
covar_list = list()
covar_list[['m']] = c('QTL_2_m', 'QTL_131_m', 'Case_control_87_m', 'Case_control_104_m', 'Case_control_80_m', 'Case_control_94_m')
covar_list[['f']] = c('QTL_4_f', 'QTL_133_f', 'Case_control_87_f', 'Case_control_104_f', 'Case_control_80_f', 'Case_control_94_f')

res = data.frame()
for (sex in c('m', 'f')){
  data1 = paste0('spracklen_nature2020_t2d_', sex); data3 = paste0('Case_control_37_', sex)
  for (data2 in covar_list[[sex]]){ # here data2 is covar
    for (pthres in c('1e-5', '5e-8')){
      file_in = paste0(path_list[['bbj']], 'result/clump/for_clump/mvmr/', pthres, '/', data1, '&', data3, '_co_', data2, '.rdata')
      file_clump = paste0(path_list[['bbj']], 'result/clump/clump_out/mvmr/', pthres, '/', data1, '&', data3, '_co_', data2, '.clumped')
      load(file_in) # df
      snp = read.table(file_clump, header=1)[,3]
      keep = rownames(df[[1]])%in%snp
      df[['exposure_beta']] = df[['exposure_beta']][keep,]; df[['exposure_pval']] = df[['exposure_pval']][keep,]; df[['exposure_se']] = df[['exposure_se']][keep,]
      df[['outcome_beta']] = df[['outcome_beta']][keep]; df[['outcome_pval']] = df[['outcome_pval']][keep]; df[['outcome_se']] = df[['outcome_se']][keep]

      out = mv_multiple(df, pval_threshold = as.numeric(pthres))[[1]]%>%select(exposure, outcome, nsnp, b, se, pval)%>%mutate(sex=sex, covar=data2, pthres=as.numeric(pthres))
      res = rbind(res, out)
    }
  }
}

write.csv(res, paste0(path_list[['bbj']], 'result/collect/t2d_cat/mvmr.csv'), row.names=F)

