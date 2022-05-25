# eqtl 
#=====================================================================================
## clean GTEx FDXR eqtl | hg38
#=====================================================================================
path="/bigdat1/user/zhanghy/gwas/summary/eqtl/raw/"
zcat ${path}GTEx_Analysis_v8_QTLs_GTEx_Analysis_v8_eQTL_all_associations_Whole_Blood.allpairs.txt.gz | \
egrep "ENSG00000161513|gene_id" > ${path}GTEx_v8_FDXR.txt

library(dplyr)
library(stringr)

path = '/bigdat1/user/zhanghy/gwas/summary/eqtl/'
file_bim = '/bigdat1/user/zhanghy/gwas/bfile/1000g/no_mhc/eur/eur_nomhc.bim'
file_map = '/bigdat1/user/zhanghy/gwas/summary/finngen/para/map_hg38_to_hg19'
file_out = '/bigdat1/user/zhanghy/gwas/summary/eqtl/clean/GTEx_v8_FDXR.txt'

bim = read.table(file_bim)[, c(2, 1, 4)]
map = read.table(file_map, header=1)
df = read.table(paste0(path, 'raw/GTEx/GTEx_v8_FDXR.txt'), header=1)

col = str_split_fixed(df$variant_id, "_", 5)
df$CHR_38 = gsub('chr', '', col[,1])
df$POS_38 = as.numeric(col[,2])
df$A1 = col[,4]
df$A2 = col[,3]

df = df%>%merge(map, by=c('CHR_38', 'POS_38'))
df = df%>%rename(BETA=slope, SE=slope_se, P=pval_nominal, FRQ=maf)
df$N = round(df$ma_count/df$FRQ)

df = df[,c('SNP', 'CHR', 'POS', 'A1', 'A2', 'FRQ', 'BETA', 'SE', 'P', 'N')]

write.table(df, file_out, sep='\t', quote=F, row.names=F)