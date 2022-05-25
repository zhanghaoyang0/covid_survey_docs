
#=====================================================================================
# prs 
#=====================================================================================
### select gwas in bfile | map a1, a2
### Most PRS software will perform strand-flipping automatically, so not need to filp beta
library(dplyr)
complement = function(string){
    string = strsplit(string, '')[[1]]
    if (sum(string%in%c('A', 'T', 'C', 'G')==F)>0){return(NA); break}
    comp = sapply(string, function(x){switch(x, "A" = "T", "C" = "G", "T" = "A", "G" = "C", return(NA))})
    comp = paste(comp, collapse='')
    return(comp)
}
data = 'spracklen_nature2020_t2d'
sex = 'f'
file_in = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/clean/', data, '_', sex, '.txt.gz')
file_bfile = paste0('/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/chn_', sex, '.bim')
file_out = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/prs/', data, '_', sex, '.txt')

bim = read.table(file_bfile)%>%rename(SNP=V2, A1=V5, A2=V6)
df_raw = read.table(file_in, header=1)%>%filter(!((A1=='T'&A2=='A')|(A1=='A'&A2=='T')|(A1=='C'&A2=='G')|(A1=='G'&A2=='C')))

info = df_raw[,c('SNP', 'A1', 'A2')]%>%merge(bim[,c('SNP', 'A1', 'A2')], 'SNP')
info = info%>%mutate(A1.z=sapply(A1.x, complement), A2.z=sapply(A2.x, complement))

match = info%>%filter(A1.x == A1.y & A2.x== A2.y)%>%pull(SNP)
rmatch = info%>%filter(A1.x == A2.y & A2.x== A1.y)%>%pull(SNP)
cmatch = info%>%filter(A1.z == A1.y & A2.z== A2.y)%>%pull(SNP)

df = df_raw%>%filter(SNP%in%c(match, rmatch, cmatch))%>%select(SNP, A1, BETA, P)
write.table(df, file_out, row.names=F, quote=F, sep='\t')

### clump
data='spracklen_nature2020_t2d'
sex='m'
file_gwas="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/prs/${data}_${sex}.txt"
file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/chn_${sex}"
file_clumped="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/prs/${data}_${sex}"

#--clump-p1 1 --clump-r2 0.1 --clump-kb 250 \

plink \
    --bfile $file_bfile \
    --clump $file_gwas \
    --clump-p1 1 --clump-r2 0.05 --clump-kb 1000 \
    --clump-snp-field SNP --clump-field P \
    --out $file_clumped

# 126738

awk 'NR!=1{print $3}' ${file_clumped}.clumped >  "${file_clumped}_clumped-snp"
awk '{print $1,$4}' $file_gwas > "${file_clumped}_pval"

file_range='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/prs/range_list'

echo "0.00000001 0 0.00000001" > $file_range 
echo "0.00001 0 0.00001" >> $file_range 
echo "0.001 0 0.001" >> $file_range 
echo "0.01 0 0.01" >> $file_range 
echo "0.05 0 0.05" >> $file_range
echo "0.1 0 0.1" >> $file_range
echo "0.5 0 0.5" >> $file_range

plink \
    --bfile $file_bfile \
    --score $file_gwas 1 2 3 header \
    --q-score-range $file_range "${file_clumped}_pval" \
    --extract "${file_clumped}_clumped-snp" \
    --out "${file_clumped}_prs"

#=====================================================================================
# prs reg
# strange, the beta of t2d prs is negative
#=====================================================================================
library(dplyr)
h2l_R2 <- function(k, r2, p) {
    x = qnorm(1-k); z = dnorm(x); i = z/k
    C = k*(1-k)*k*(1-k)/(z^2*p*(1-p))
    theta = i*((p-k)/(1-k))*(i*((p-k)/(1-k))-x)
    h2l_R2 = C*r2 / (1 + C*theta*r2)
}
path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/prs/'
df_raw = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/data/icd10_clean.txt', header=1)%>%mutate(cat = ifelse(grepl('H25|H26|H28', icd), 1, 0))


out = c()
for (sex in c('m', 'f')){
    for (thres in c('0.5', '0.1', '0.05', '0.001', '0.00001', '0.00000001')){
        prs = read.table(paste0(path, 'spracklen_nature2020_t2d_', sex, '_prs.', thres, '.profile'), header=1)%>%select(FID, SCORE)
        df = df_raw%>%select(-IID, icd)%>%merge(prs, 'FID')%>%mutate(batch=as.factor(batch), center=as.factor(center))%>%filter(age>50) # temp
        lr_full <- glm(cat ~ SCORE + age + batch + center , data = df, family = binomial())
        lr_base <- glm(cat ~  age + batch + center, data = df, family = binomial())
        out = c(out, sex, thres, summary(lr_full)$coefficients[2, c(1, 2, 4)])
    }
}
res = data.frame(matrix(out, ncol=5, byrow=T))
names(res) = c('sex', 'prs_thres', 'b', 'se', 'p')
res

