## phewas 
#=====================================================================================
## chr17:72,858,619-72,869,156
## key snp: 
## 17	rs690514	72862593	T	C
#=====================================================================================
## get snp list of FDXR 
file_bim = '/bigdat1/user/zhanghy/gwas/ukb/bfile/brt.bim'
file_eqtl = '/bigdat1/user/zhanghy/gwas/summary/eqtl/clean/GTEx_v8_FDXR.txt'
file_snp = '/bigdat1/user/zhanghy/gwas/ukb/fdxr/FDXR_snplist'

eqtl = read.table(file_eqtl, header=1)
eqtl$P_FDR = p.adjust(eqtl$P,method="fdr",n=nrow(eqtl))
bim = read.table(file_bim)
snp1 = bim[bim[,1]==17&bim[,4]>=72858619&bim[,4]<=72869156, 2]
snp2 = eqtl[eqtl$P_FDR<0.05, 'SNP']

snp = unique(c(snp1, snp2))

write.table(snp, file_snp, row.names=F, col.names=F, quote=F)

## check
file_check = '/bigdat1/user/zhanghy/gwas/ukb/fdxr/check_snplist'
ldl = read.table('/bigdat1/user/zhanghy/gwas/summary/gwascatalog/raw/ldl/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz', header=1)

ldl = ldl[nchar(ldl$REF)==1&nchar(ldl$ALT)==1,]
ldl = ldl[ldl$pvalue<1e-5&ldl$EFFECT_SIZE>0,]
ldl = ldl[order(ldl$EFFECT_SIZE, decreasing=T),]

snp3 = ldl$rsID[1:5000]
write.table(snp3, file_check, row.names=F, col.names=F, quote=F)


## prepare genotype
# recodeA make genotype to 0,1,2 (number of MA)
# recodeD make genotype to 0,1 (het or not)
# recodeAD make genotype to both col
file_bfile="/bigdat1/user/zhanghy/gwas/ukb/bfile/brt"
file_snp="/bigdat1/user/zhanghy/gwas/ukb/fdxr/check_snplist"
file_out="/bigdat1/user/zhanghy/gwas/ukb/fdxr/check_recodeA"

plink \
--recode A \
--bfile $file_bfile \
--extract $file_snp \
--out $file_out


### phewas
# nominal 
# genotypes2[, j] = as.character(genotypes2[, j])
# out = phewas(phenotypes=phenotypes, genotypes=genotypes2,
# covariates=covar, additive.genotypes = F) 

library(PheWAS)
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/gwas/ukb/fdxr/'

args = commandArgs(T)
i = as.numeric(args[1])

# prepare rdata for phewas, done
if (F){
  file_pc  = paste0(path, '../bfile/brt.eigenvec')
  file_info = paste0(path, '../data/brt_pheno_freq>0.001.csv')
  file_allinfo = paste0(path, '../para/phewas_pheno_code.csv')
  
  info = read.csv(file_info)
  allinfo = read.csv(file_allinfo)
  pc = read.table(file_pc)[,2:12]
  names(pc) = c('id', paste0('pc', 1:10))
  
  samples_keep = pc$id
  
  geno = read.table(paste0(path, 'FDXR_recodeA.raw'),header=TRUE)
  geno = geno[,c(1,7:ncol(geno))]
  
  pheno = read.table(paste0(path, '../data/icd10_body_self-report_blood_clean.txt'), header=1)
  pheno$traits = paste0(pheno$icd, ';', pheno$selfreport)
  names(geno)[1] = names(pheno)[1] = "id"
  geno = geno[geno$id%in%samples_keep,]
  pheno = pheno[pheno$id%in%samples_keep,]
  pheno = pheno%>%merge(pc, by='id')
  
  # add self-define traits, done
  # self-define traits
  traits_add = c('pad', 'cataract', 'ihd')
  if (F){
    # cataract
    cataract_code = c('H25', 'H26', 'H28')
    pheno$cataract = grepl('H25|H26|H28', pheno$traits)
    # pad
    pad_code = c('I7020','I7021','I7080','I7081','I7090','I7091','I721','I724','I730',
                 'I731','I738','I739', 'I742','I743','I744','I745','I748','I749','I770',
                 'I776','I778','I779','I792')
    pheno$pad = grepl(paste0(pad_code, collapse='|'), pheno$traits)
    # ihd, Ischaemic heart diseases
    ihd_code = paste0('I', 20:25)
    pheno$ihd = grepl(paste0(ihd_code, collapse='|'), pheno$traits)
  }
  pheno[,c('icd', 'selfreport')] = NULL
  
  traits_numeric = c('height', 'sitting_height', 'weight', 'bmi', 'hc', 'wc', names(pheno)[grepl('blood', names(pheno))])
  
  # all traits
  traits = c(info[,1], traits_numeric, traits_add)
  save(pheno, geno, info, traits, file=paste0(path, 'temp_phewas.rdata'))
}

load(paste0(path, 'temp_phewas.rdata'))
covar = pheno[,c('id', 'sex', 'age', 'center', 'batch', paste0('pc', 1:10))]
covar[,c('center', 'batch')] = sapply(covar[,c('center', 'batch')], as.character)

traits1 = traits[grepl('F', traits)]
# traits1 = traits

for (trait in traits1){
  file_out = paste0(path, 'result/phewas/', trait, '.csv')
  if (!file.exists(file_out)){
    cols = c("phenotype", "snp", "adjustment", "beta", "SE", "OR", "p", "type", "n_total", "n_cases", "n_controls", "HWE_p","allele_freq", "n_no_snp", 
             "note", "model", "sex")
    cat(cols, sep=',', file=file_out, append=F)
    cat('\n', file=file_out, append=T)
  }
  if (trait %in% info[,1]){
    sexs = c('male', 'female')
    # sexs = sexs[info[info[,1]==trait, c('prop_case_b')]>0.01] # sex is intersectial with geno
    pheno$y =  grepl(trait, pheno$traits)
  } else { # numeric traits
    sexs = c('male', 'female')
    pheno$y = pheno[,trait]
  } 
  for (sex in sexs){
    print(sex)
    if (sex=='male'){
      id_keep = covar[covar$sex==1, 'id']
    } else if (sex=='female'){
      id_keep = covar[covar$sex==0, 'id']
    } else{
      id_keep = covar[, 'id']
    }
    phenotypes = pheno[pheno$id%in%id_keep, c('id', 'y', 'traits')]
    if (grepl('G|F', trait)){
      drop = phenotypes$y==F&grepl('G|F',phenotypes$traits)
      phenotypes = phenotypes[!drop,]
    }
    phenotypes = phenotypes[,c('id', 'y')]
    for (snp in names(geno)[2:ncol(geno)]){
      print(snp)
      genotypes_raw = geno[, c('id', snp)]
      for (model in c('additive', 'recessive')){
        print(model)
        if (sum(grepl(snp, readLines(file_out))& 
                grepl(sex, readLines(file_out))& 
                grepl(model, readLines(file_out)))>0) {next}
        
        genotypes = genotypes_raw
        if (model == 'recessive'){
          genotypes[!is.na(genotypes[,2]) & genotypes[,2]==1, 2] = 0
          genotypes[!is.na(genotypes[,2]) & genotypes[,2]==2, 2] = 1
        }else if (model == 'dominant'){
          genotypes[!is.na(genotypes[,2]) & genotypes[,2]==1, 2] = 1
          genotypes[!is.na(genotypes[,2]) & genotypes[,2]==2, 2] = 1
        }
        
        covariates = covar
        if (sex!='both_sex'){
          covariates$sex = NULL
        } 
        # interact term
        if (sex=='both_sex'){ 
          covariates = covariates%>%merge(genotypes)
          covariates[,snp] = covariates[,snp]*covariates$sex
          covariates = covariates%>%rename(geno_sex=paste0(snp))
        }
        
        out = phewas(phenotypes=phenotypes, genotypes=genotypes,
                      covariates=covar, additive.genotypes = T) 
        
        out$model = model
        out$phenotype = trait
        out$sex = sex
        cat(unlist(out[1,]), sep=',', file=file_out, append=T)
        cat('\n', file=file_out, append=T)
      }
    }
  }
}


### phewas.sh
#!/usr/bin/bash
i="$1"
Rscript phewas.r $i


for i in {31..39}  
do
node=9
sbatch -N1 -n1 -c2 -w compute$node phewas.sh $i
done



## collect
library(dplyr)
library(stringr)
path = '/bigdat1/user/zhanghy/gwas/ukb/fdxr/'
file_bim = '/bigdat1/user/zhanghy/gwas/ukb/bfile/brt_qc/qc_1_bychr/chr17_1.bim'
file_eqtl = '/bigdat1/user/zhanghy/gwas/summary/eqtl/clean/GTEx_v8_FDXR.txt'
file_code = paste0(path, '../para/phewas_pheno_code.csv')

code = read.csv(file_code)
eqtl = read.table(file_eqtl, header=1)[,c('SNP', 'A1', 'A2', 'BETA', 'P')]%>%rename(snp=SNP, eqtl_beta=BETA, eqtl_p=P)
temp = eqtl
t = temp$A1
temp$A1 = temp$A2
temp$A2 = t
temp$eqtl_beta = -temp$eqtl_beta
eqtl = rbind(eqtl, temp)

bim = read.table(file_bim)[,c(2,5,6)]
names(bim) = c('snp', 'A1', 'A2')

traits = list.files(paste0(path, 'result/phewas'))
# traits = traits[grepl('check', traits)]

df = data.frame()
for (trait in traits){
  temp = read.csv(paste0(path, 'result/phewas/', trait))
  df = rbind(df, temp)
}

t = df[grepl('ihd', df$phenotype),]

col = str_split_fixed(df$snp, '_', 2)
df$snp = col[,1]
df$A1 = col[,2]
dim(df)

df$A2 = NA
for (i in unique(df$snp)){
  A1 = df[df$snp==i, 'A1'][1]
  A2 = setdiff(unlist(bim[bim$snp==i, c('A1', 'A2')]), A1)
  df[df$snp==i, 'A2'] = A2
}
sum(is.na(df$A2))
dim(df)

df = df%>%merge(eqtl)
df = df%>%rename(coding=phenotype)
df = df%>%merge(code, by='coding', all.x=T)
dim(df)

df$t = paste0(df$coding, '_', df$sex, '_', df$model)

df1 = data.frame()
for (i in unique(df$t)){
  sub = df[df$t==i,]
  sub$fdr.p = p.adjust(df[df$t==i, 'p'], 'fdr')
  df1 = rbind(df1, sub)
}
df1$t = NULL

t = df[df$p<0.05,]
t[t$snp=='rs690514',]

write.csv(df1, paste0(path, 'result/collect/phewas.csv'), row.names = F)
