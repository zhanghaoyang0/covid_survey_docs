############################################################################
# MR
############################################################################
#### 第一步，选出暴露p<5e-8的跟结局合并在一起。
library(TwoSampleMR)
library(dplyr)

df1_raw = read.csv('/home/yanglab_data/user/xiuxh/gwas/ecg/batch3/bolt/out/lead10_114.tab', sep='\t')
df2_raw = read.csv('/home/yanglab_data/user/zhanghy/gwas/summary/finngen/clean/finr5_348.txt.gz', sep='\t')

df1 = df1_raw%>%filter(P_BOLT_LMM_INF<5e-8)%>%
    rename(pval.exposure=P_BOLT_LMM_INF, effect_allele.exposure=ALLELE1, other_allele.exposure=ALLELE0, beta.exposure=BETA, se.exposure=SE, eaf.exposure=A1FREQ)%>%
    mutate(id.exposure='I_214', exposure='I_214', samplesize.exposure=34519)
df2 = df2_raw%>%
    rename(pval.outcome=P, effect_allele.outcome=A1, other_allele.outcome=A2, samplesize.outcome=N, beta.outcome=BETA, se.outcome=SE, eaf.outcome=FRQ)%>%
    mutate(id.outcome='cahd', outcome='cahd')

dat = harmonise_data(df1, df2)

df_for_clump = dat%>%rename(P = pval.exposure)%>%select(SNP, P)
write.table(df_for_clump, paste0('/home/zhanghy/wangxf/df_for_clump.txt'), row.names=F, sep='\t', quote=F)

### 第二步，选独立的SNP
file_bfile="/home/yanglab_data/user/zhanghy/gwas/bfile/1000g/no_mhc/eur/eur_nomhc"
file_in="/home/zhanghy/wangxf/df_for_clump.txt"
file_out="/home/zhanghy/wangxf/df_for_clump"

/home/yanglab_data/user/zhanghy/soft_slurm/plink_1.9/plink \
--allow-no-sex \
--clump $file_in --bfile $file_bfile --out $file_out \
--clump-kb 1000 --clump-p1 1.0 --clump-p2 1.0 --clump-r2 0.05 

### 第三步，用clump后的SNP做MR
snp = read.table('/home/zhanghy/wangxf/df_for_clump.clumped', header=1)[,3]
df_mr = dat%>%filter(SNP%in%snp)
mr(df_mr)


############################################################################
# GSMR
############################################################################
#### 第一步，选出暴露p<5e-8的跟结局合并在一起, R环境
library(TwoSampleMR)
library(dplyr)

df1_raw = read.csv('/home/yanglab_data/user/xiuxh/gwas/ecg/batch3/bolt/out/I_214.tab', sep='\t')
df2_raw = read.csv('/home/yanglab_data/user/zhanghy/gwas/summary/finngen/clean/finr5_348.txt.gz', sep='\t')

df1 = df1_raw%>%filter(P_BOLT_LMM_INF<5e-8)%>%
    rename(pval.exposure=P_BOLT_LMM_INF, effect_allele.exposure=ALLELE1, other_allele.exposure=ALLELE0, beta.exposure=BETA, se.exposure=SE, eaf.exposure=A1FREQ)%>%
    mutate(id.exposure='I_214', exposure='I_214', samplesize.exposure=34519)
df2 = df2_raw%>%
    rename(pval.outcome=P, effect_allele.outcome=A1, other_allele.outcome=A2, samplesize.outcome=N, beta.outcome=BETA, se.outcome=SE, eaf.outcome=FRQ)%>%
    mutate(id.outcome='cahd', outcome='cahd')

dat = harmonise_data(df1, df2) # 你可以把这个harmoniise的数据保存下载，gsmr跟mr都要用。

df_for_clump = dat%>%rename(P = pval.exposure)%>%select(SNP, P)
write.table(df_for_clump, paste0('/home/zhanghy/wangxf/df_for_clump.txt'), row.names=F, sep='\t', quote=F)

### 第二步，选独立的SNP,生成这些snp的LD矩阵, bash 环境
## clump
file_bfile="/home/yanglab_data/user/zhanghy/gwas/bfile/1000g/no_mhc/eur/eur_nomhc"
file_in="/home/zhanghy/wangxf/df_for_clump.txt"
file_out="/home/zhanghy/wangxf/df_for_clump"

/home/yanglab_data/user/zhanghy/soft_slurm/plink_1.9/plink \
--allow-no-sex \
--clump $file_in --bfile $file_bfile --out $file_out \
--clump-kb 1000 --clump-p1 1.0 --clump-p2 1.0 --clump-r2 0.05 

## ld mat
file_snp="/home/zhanghy/wangxf/df_for_clump.txt"
file_mat="/home/zhanghy/wangxf/ld_mat"
awk '{print $3}' ${file_out}.clumped | grep -v SNP > $file_snp

plink --bfile $file_bfile \
--extract $file_snp \
--r2 square \
--write-snplist \
--out $file_mat

### 第三步，用clump后的SNP做GSMR, R环境
library(gsmr)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source_mr.r')
get_gsmr_para('gsmr', 'eur') # 提取参数设置，你可以打开这个函数看看具体设置
nsnps_thresh = 1
pthres = 5e-8

file_snp = '/home/zhanghy/wangxf/df_for_clump.clumped'
file_ldmat = '/home/zhanghy/wangxf/ld_mat.ld'
file_ldmat_snp = '/home/zhanghy/wangxf/ld_mat.snplist'

snp = read.table(file_snp, header=1)[,3]
mat = read.table(paste0(file_ldmat))
snp_list = read.table(paste0(file_ldmat_snp))[,1]
colnames(mat) = rownames(mat) = snp_list

df_mr = dat%>%filter(SNP%in%snp)

gsmr_data = df_mr%>%rename('a1'='effect_allele.exposure', 'a2'='other_allele.exposure', 'bzx_pval'='pval.exposure', 'bzy_pval'='pval.outcome', 'bzx'='beta.exposure',
                    'bzy'='beta.outcome', 'bzx_se'='se.exposure', 'bzy_se'='se.outcome', 'bzx_n'='samplesize.exposure', 'bzy_n'='samplesize.outcome', 'a1_freq'='eaf.exposure')

ldrho = mat[rownames(mat) %in% gsmr_data$SNP, colnames(mat)%in%gsmr_data$SNP]
snp_coeff_id = rownames(ldrho)

gsmr_results = try(gsmr(gsmr_data$bzx, gsmr_data$bzx_se, gsmr_data$bzx_pval, gsmr_data$bzy, gsmr_data$bzy_se, gsmr_data$bzy_pval,
                        ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh = as.numeric(pthres), single_snp_heidi_thresh, 
                        multi_snp_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta))

############################################################################
# CAUSE |  这个方法最慢,建议多任务
############################################################################
library(cause)
library(dplyr)
df1_raw = read.csv('/home/yanglab_data/user/xiuxh/gwas/ecg/batch3/bolt/out/I_214.tab', sep='\t')
df2_raw = read.csv('/home/yanglab_data/user/zhanghy/gwas/summary/finngen/clean/finr5_348.txt.gz', sep='\t')
path_cause_ref = '/home/yanglab_data/user/zhanghy/project/mr_server/db/cause_ref/eur/chr'

df1 = df1_raw%>%rename(POS=BP, A1=ALLELE1, A2=ALLELE0, FRQ=A1FREQ, P=P_BOLT_LMM_INF) # 重命名，XXH的GWAS都要重命名。
df <- gwas_merge(df1, df2_raw, snp_name_cols = c("SNP", "SNP"), beta_hat_cols = c("BETA", "BETA"), 
    se_cols = c("SE", "SE"), A1_cols = c("A1", "A1"), A2_cols = c("A2", "A2"))

## 任务1，选择SNP做CAUSE
variants <- df %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
pruned = c()
for (chr in 1:22){ # 这里也可以多任务
    print(paste0('chr ', chr))
    ld <- readRDS(paste0(path_cause_ref, chr, '_AF0.05_0.1.RDS'))
    snp_info <- readRDS(paste0(path_cause_ref, chr, '_AF0.05_snpdata.RDS'))
    out <- ld_prune(variants = variants, ld = ld, total_ld_variants = snp_info$SNP, pval_cols = c("pval1"),pval_thresh = c(1e-3))
    pruned = c(pruned, out)
    print(length(pruned))
}

## 任务2，估计初始参数
set.seed(100)
varlist <- with(df, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(df, varlist)

## 最后做cause
res = try(cause(X=df, variants = pruned, param_ests = params))


############################################################################
# LDSC
############################################################################
### 第一步，修改心电数据列名, R 环境
library(dplyr)
df_raw = read.csv('/home/yanglab_data/user/xiuxh/gwas/ecg/batch3/bolt/out/I_214.tab', sep='\t')
df = df_raw%>%rename(POS=BP, A1=ALLELE1, A2=ALLELE0, FRQ=A1FREQ, P=P_BOLT_LMM_INF)%>%mutate(N=34519)%>%
    select(SNP, CHR, POS, FRQ, A1, A2, BETA, SE, P, N)

write.table(df, '/home/zhanghy/wangxf/I_214.txt', quote=F, sep='\t', row.names=F)

### 第二步，munge, Linux 环境
export PATH=/home/yanglab/zhanghy/miniconda2/bin/:$PATH 

file_in="/home/zhanghy/wangxf/I_214.txt"
file_out="/home/zhanghy/wangxf/I_214"
/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/munge_sumstats.py \
--sumstats $file_in \
--out $file_out \
--merge-alleles /home/yanglab_data/user/zhanghy/soft_slurm/ldsc/w_hm3.snplist

### 第三步，LDSC，Linux 环境
path_ld='/home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr/'
file_in_1="/home/zhanghy/wangxf/I_214"
file_in_2="/home/yanglab_data/user/zhanghy/gwas/summary/finngen/munge/finr5_348"
file_out="/home/zhanghy/wangxf/I_214_and_finr5_348"

/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/ldsc.py \
--rg  ${file_in_1}.sumstats.gz,${file_in_2}.sumstats.gz \
--ref-ld-chr $path_ld \
--w-ld-chr $path_ld \
--out $file_out




### res 
/home/yanglab_data/user/wangxf/ecg_gwas/dat_mr

datas = list.files('./')
datas = files[grepl('tab', files)]

collect = data.frame()
for (data in datas){
    df = read.table(data)
}