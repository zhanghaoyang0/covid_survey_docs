## prs 
#=====================================================================================
# qc
#=====================================================================================
### base data
path = '/bigdat1/user/zhanghy/gwas/summary/eqtl/'

df = read.table(paste0(path, 'clean/GTEx_v8_FDXR.txt'), header=1)
# remove ambiguous SNPs
drop = (df$A1=='A'&df$A2=='T'|
          df$A1=='T'&df$A2=='A'|
          df$A1=='C'&df$A2=='G'|
          df$A1=='G'&df$A2=='C' )
df = df[!drop,]

write.table(df, paste0(path, 'prs/for_prs/GTEx_v8_FDXR.txt'), sep='\t', row.names=F, quote=F)
write.table(df$SNP, paste0(path, 'prs/for_prs/temp.snp'), sep='\t', row.names=F, col.names=F, quote=F)
write.table(df[,c('SNP', 'A1')], paste0(path, 'prs/for_prs/GTEx_v8_FDXR.a1'), sep='\t', row.names=F, col.names=F, quote=F)

### target
## get bim for map
file_bfile="/bigdat1/user/zhanghy/gwas/ukb/bfile/brt"
file_extract_snp="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/temp.snp"
file_out="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/temp"

plink \
--bfile $file_bfile \
--make-just-bim \
--extract $file_extract_snp \
--out $file_out 

wc -l ${file_out}.bim

## get mapped snp 
file_gwas = '/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/GTEx_v8_FDXR.txt'
file_bim = "/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/temp.bim"
file_snp = "/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/ukb-brt_map_GTEx-v8-FDXR.snp"

df = read.table(file_gwas, header=1)
bim = read.table(file_bim)

keep = c()
for (i in 1:nrow(bim)){
  alleles1 = df[df$SNP==bim[1, 2], c('A1', 'A2')]
  alleles2 = bim[i, c(5,6)]
  if (length(intersect(alleles2, alleles1))==2){
    keep = c(keep, i)
  }
}
bim = bim[keep,]
dim(bim)

write.table(bim[,2], file_snp, row.names=F, col.names=F, quote=F, sep='\t')

## get bfile 
file_a1="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/GTEx_v8_FDXR.a1"
file_bfile="/bigdat1/user/zhanghy/gwas/ukb/bfile/brt"
file_extract_snp="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/ukb-brt_map_GTEx-v8-FDXR.snp"
file_out="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/ukb-brt_map_GTEx-v8-FDXR"

plink \
--bfile $file_bfile \
--a1-allele $file_a1 \
--make-bed \
--extract $file_extract_snp \
--out $file_out 

wc -l ${file_out}.bim

#=====================================================================================
# prs
#=====================================================================================
## clump
file_bfile="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/ukb-brt_map_GTEx-v8-FDXR"
file_gwas="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/GTEx_v8_FDXR.txt"
file_out="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/GTEx_v8_FDXR"

plink \
--bfile $file_bfile \
--clump-p1 0.5 \
--clump-r2 0.5 \
--clump-kb 250 \
--clump $file_gwas \
--clump-snp-field SNP \
--clump-field P \
--out $file_out

# rs690514	17	72862593 3.08346e-05
# rs541242	17	72871095 6.29737e-06

# r2 matrix
file_bfile="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/ukb-brt_map_GTEx-v8-FDXR"
file_out="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/ukb-brt_map_GTEx-v8-FDXR"

plink --bfile $file_bfile \
--r2 \
--out $file_out

awk 'NR==1||$3=="rs690514"' ukb-brt_map_GTEx-v8-FDXR.ld

### range
file_gwas="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/GTEx_v8_FDXR.txt"
file_list="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/range_list"
file_pval="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/GTEx_v8_FDXR.pval"

echo "0.001 0 0.001" > $file_list 
echo "0.01 0 0.01" > $file_list 
echo "0.05 0 0.05" >> $file_list 
echo "0.1 0 0.1" >> $file_list 
echo "0.3 0 0.3" >> $file_list 
echo "0.5 0 0.5" >> $file_list 

awk '{print $1,$9}' $file_gwas > $file_pval

## prs
# 1,4,7: snp id, a1, beta | header: with header
file_bfile="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/ukb-brt_map_GTEx-v8-FDXR"
file_gwas="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/GTEx_v8_FDXR.txt"
file_list="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/range_list"
file_pval="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/for_prs/GTEx_v8_FDXR.pval"
file_out="/bigdat1/user/zhanghy/gwas/summary/eqtl/prs/prs_out/GTEx_v8_FDXR.prs"

plink \
--bfile $file_bfile \
--score $file_gwas 1 4 7 header \
--q-score-range $file_list $file_pval \
--out $file_out
#=====================================================================================
# prs regression | chn 
#=====================================================================================
## prs_mod.r
library(dplyr)
library(pROC)

args = commandArgs(T)
i = as.numeric(args[1])

h2l_R2 <- function(k, r2, p) {
  x = qnorm(1-k)
  z = dnorm(x)
  i = z/K
  C = k*(1-k)*k*(1-k)/(z^2*p*(1-p))
  theta = i*((p-k)/(1-k))*(i*((p-k)/(1-k))-x)
  h2l_R2 = C*r2 / (1 + C*theta*r2)
} # k is pop prev, p is sample p

# #Nagelkerke's R2 
# LLf = logLik(lr_full)      # full model
# LLc = logLik(lr_base)    # base model with covs
# Ng_R2v = (1-exp((2/N)*(LLc[1]-LLf[1])))/(1-exp((2/N)*LLc[1]))
# h2l_r2v=h2l_R2(K,Ng_R2v,P) 

path = '/bigdat1/user/zhanghy/gwas/ukb/fdxr/'
load(paste0(path, 'temp_phewas.rdata'))
df = pheno%>%merge(covar[,c('id', paste0('pc', 1:10))])

# thress = c(0.001, 0.05, seq(0.1, 0.5, 0.1))
traits1 = traits[grepl('I', traits)]
# traits1 = traits

for (trait in traits1){
  print(trait)
  if (trait %in% info[,1]){
    df$y = grepl(trait, df$traits)
  } else { # numeric traits
    df$y = df[,trait]
  }
  res = data.frame()
  
  for (sex in c('male', 'female')){
    if (sex=='both_sex') {
      df1 = df
      base_formula = as.formula(paste('y~sex+age+batch+center+', paste0('pc', c(1:10), collapse='+')))
      full_formula = as.formula(paste('y~prs+prs*sex+sex+age+batch+center+', paste0('pc', c(1:10), collapse='+')))
    } else {
      base_formula = as.formula(paste('y~age+batch+center+', paste0('pc', c(1:10), collapse='+')))
      full_formula = as.formula(paste('y~prs+age+batch+center+', paste0('pc', c(1:10), collapse='+')))
      if (sex=='male'){
        df1 = df[df$sex==1,]
      } else {
        df1 = df[df$sex==0,]
      }
    }
    if (class(df$y)=='logical'){
      mod_base = glm(base_formula, data = df1, family = "binomial")
    } else {
      mod_base = lm(base_formula, data = df1)
    }
    for (thres in c(0.01, 0.05, 0.1)){
      file_out = paste0(path, 'result/prs/', trait, '_', sex, '_', thres, '.rdata')
      if (file.exists(file_out)) {
        next
      }
      print(thres)
      file_score = paste0(path, '../../summary/eqtl/prs/prs_out/GTEx_v8_FDXR.prs.', thres, '.profile')
      score = read.table(file_score, header=1)[,c(1,6)]%>%rename(id=FID, prs=SCORE)
      df2 = df1%>%merge(score, by='id')
      
      if (class(df$y)=='logical'){
        mod_full = glm(full_formula, data = df2, family = "binomial")
      } else {
        mod_full = lm(full_formula, data = df2)
      }
      # log likelihood
      ll_base = logLik(mod_base)  
      ll_full = logLik(mod_full)    
      # auc
      pred_base = predict(mod_base, df2, type="response")
      pred_full = predict(mod_full, df2, type="response")
      auc_base = auc(df2$y, pred_base)[1]
      auc_full = auc(df2$y, pred_full)[1]
      
      coef = as.data.frame(summary(mod_full)$coefficients)
      coef = coef[!grepl('batch|center|pc|Inter|sex', rownames(coef)), c(1,2,4)]
      
      if (sum(grepl('binomial', summary(mod_full)$call))>0){
        reg_class = 'logistic'} else {reg_class = 'linear'}
      
      res = c(trait, thres, reg_class, auc_base, auc_full, ll_base, ll_full, as.vector(t(coef)))
      names(res) = c('trait', 'prs_pthres', 'reg', 'auc_base', 'auc_full', 'log_likhood_base', 'log_likhood_full',
                     as.vector(sapply(c('prs', 'age'), function(x){paste0(x, c('_b', '_se', '_p'))})))
      save(res, file=file_out)
    }
  }
}


      

## prs_mod.sh
#!/usr/bin/bash
i="$1"
Rscript prs_mod.r $i

## run
for i in {36..39}  
do
node=8
sbatch -N1 -n1 -c2 -w compute$node prs_mod.sh $i
done


## collect
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/gwas/ukb/fdxr/result/'

file_code = paste0(path, '../../para/phewas_pheno_code.csv')
code = read.csv(file_code)%>%rename(trait=coding)

files = list.files(paste0(path, 'prs'))
# files = files[grepl('I', files)]

collect = c()
for (file in files){
  load(paste0(path, 'prs/', file))
  sex = strsplit(file, '_')[[1]][2]
  collect = c(collect, res, sex)
}
collect = data.frame(matrix(collect, ncol=(length(res)+1), byrow=T))
names(collect) = c(names(res), 'sex')

collect[,4:11] = sapply(collect[,4:11], as.numeric)
collect = collect%>%merge(code, by='trait')

collect[collect$prs_p<0.05,]

write.csv(collect, paste0(path, 'collect/prs.csv'), row.names = F)

        
#=====================================================================================
# plot
#=====================================================================================
library(forestplot)
library(grid)

path = 'D:/nutstore/project/gwas/FDXR/data/'

df = read.csv(paste0(path, 'prs.csv'))
df = df[nchar(df$trait)!=3,]
df = df[order(df$prs_b),]

df$prs_b_l = df$prs_b-1.96*df$prs_se
df$prs_b_u = df$prs_b+1.96*df$prs_se

keep = 'G|F'
pthres = 0.01 # 0.1 0.05 0.01

temp = df[grepl(keep, df$trait)&df$prs_pthres==pthres,];temp$trait

sub1 = temp[temp$sex=='male', c('trait', 'meaning', 'prs_b', 'prs_b_l', 'prs_b_u', 'prs_p')]
sub2 = temp[temp$sex=='female', c('trait', 'meaning', 'prs_b', 'prs_b_l', 'prs_b_u', 'prs_p')]
names(sub1)[3:6] = paste0(names(sub1)[3:6], '.m')
names(sub2)[3:6] = paste0(names(sub2)[3:6], '.f')
sub = sub1%>%merge(sub2, by=c('trait', 'meaning'))
# sub = sub[order(sub$trait),]

tabletext = as.list(c('Disease', unlist(sub$meaning)))
legend = c('male', 'female')


if (pthres == 0.01){
  x_min = -8
  x_max = 8
  step = 2
} else if (pthres == 0.05){
  x_min = -15
  x_max = 15
  step = 5
} else if (pthres == 0.1){
  x_min = -20
  x_max = 25
  step = 5
}

forestplot(tabletext,
           mean  = cbind(c(NA, sub$prs_b.m), c(NA, sub$prs_b.f)),
           lower = cbind(c(NA, sub$prs_b_l.m), c(NA, sub$prs_b_l.f)),
           upper = cbind(c(NA, sub$prs_b_u.m), c(NA, sub$prs_b_u.f)),
           new_page = T, boxsize=0.25, line.margin=0.6, 
           legend = legend,
           legend_args = fpLegend(pos = list(x=1, y=0.98), gp=gpar(col="#CCCCCC", fill="#F9F9F9")),
           is.summary=c(T, rep(F, nrow(sub))), 
           cex=15.5, vertices = TRUE, xlog=FALSE, zero=0, 
           #clip=c(x_min, x_max), 
           xlab="Effect of PRS for FDXR expression on diseases", 
           col=fpColors(box=c("blue", "red"), line="darkblue", summary="royalblue"),
           #xticks = seq(from = x_min, to = x_max, by = step),
           txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", cex=1),
                            xlab  = gpar(fontfamily = "", cex = 1),
                            title = gpar(fontfamily = "", cex=1, align='l'),
                            legend = gpar(fontfamily = "", cex = 1)))


t = df[grepl('G|F', df$trait)&df$prs_p<0.05, c('meaning', 'sex', 'prs_pthres', 'prs_b', 'prs_se', 'prs_p', 'auc_base', 'auc_full')]
write.csv(t, 'C:/Users/Admin/Desktop/t.csv')
