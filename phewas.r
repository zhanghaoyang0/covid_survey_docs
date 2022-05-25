## phewas 
#=====================================================================================
# genecard hg19 Entrez Gene
info = setNames(data.frame(matrix(c('TSC2', 16, 2097986, 2139492, 'HUWE1', 23, 53559057, 53713664, 'USP7', 16, 8985954, 9057763, 'TRIP12', 2, 230628553, 230787902,'FMR1', 23, 146993437, 147032645, 
'HACE1', 6, 105175969, 105307794, 'BRCA1', 17, 41196312, 41277381,'MKRN2', 3, 12598586, 12625212, 'CHD8', 14, 21853358, 21924282, 'FDXR', 17, 72858619, 72869119), ncol=4, byrow=T)), c('gene', 'chr', 'start', 'end'))
# gene chr     start       end
# TSC2  16   2097986   2139492
# HUWE1  23  53559057  53713664
# USP7  16   8985954   9057763
# TRIP12   2 230628553 230787902
# FMR1  23 146993437 147032645
# HACE1   6 105175969 105307794
# BRCA1  17  41196312  41277381
# MKRN2   3  12598586  12625212
# CHD8  14  21853358  21924282
# FDXR  17  72858619  72869119
# 没有TRP12, 只有TRIP12
# 没有FMRP, 只有FMR1 (FMRP translational regulator 1)
#=====================================================================================
## get snp list of FDXR 
source('/home/yanglab_data/user/zhanghy/project/temp/code/source_phewas.r')

file_bim = '/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt.bim'
file_snp = '/home/yanglab_data/user/zhanghy/gwas/ukb/phewas/data/snp'

bim = read.table(file_bim)

snp_list = list()
for (i in 1:nrow(info)){
  chr = info[i, 'chr']; start = info[i, 'start']; end = info[i, 'end']
  snp = bim%>%filter(V1==chr&V4>=start&V4<=end)%>%pull(V2)
  snp_list[[info[i, 'gene']]] = snp
  print(i); print(length(snp))
}
save(snp_list, file='/home/yanglab_data/user/zhanghy/gwas/ukb/phewas/data/snp_list.rdata')

# file_eqtl = '/home/yanglab_data/user/zhanghy/gwas/summary/eqtl/clean/GTEx_v8_FDXR.txt'
# eqtl = read.table(file_eqtl, header=1)
# eqtl$P_FDR = p.adjust(eqtl$P,method="fdr",n=nrow(eqtl))
# snp2 = eqtl[eqtl$P_FDR<0.05, 'SNP']
# snp = unique(c(snp1, snp2))

write.table(unique(unlist(snp_list)), file_snp, row.names=F, col.names=F, quote=F)


## prepare genotype
# recodeA make genotype to 0,1,2 (number of MA)
# recodeD make genotype to 0,1 (het or not)
# recodeAD make genotype to both col
file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt"
file_snp="/home/yanglab_data/user/zhanghy/gwas/ukb/phewas/data/snp"
file_out="/home/yanglab_data/user/zhanghy/gwas/ukb/phewas/data/snp_recodeA"

plink \
--recode A \
--bfile $file_bfile \
--extract $file_snp \
--out $file_out


### phewas
# nominal 
# genotypes2[, j] = as.character(genotypes2[, j])
# out = phewas(phenotypes=phenotypes, genotypes=genotypes2, covariates=covar, additive.genotypes = F) 

library(PheWAS)
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/gwas/ukb/phewas/'

args = commandArgs(T)
i = as.numeric(args[1])

# prepare rdata for phewas, done
if (F){
  file_pc  = paste0(path, '../bfile/brt.eigenvec')
  file_info1 = paste0(path, '../data/brt_pheno_freq>0.001.csv')
  file_info2 = paste0(path, '../para/phewas_pheno_code.csv')
  
  info1 = read.csv(file_info1)
  info2 = read.csv(file_info2)
  pc = read.table(file_pc)[,2:12]; names(pc) = c('id', paste0('pc', 1:10))
  samples_keep = pc$id
  
  geno_raw = read.table(paste0(path, 'data/snp_recodeA.raw'),header=TRUE)
  geno = geno_raw%>%select(-IID, -PAT, -MAT, -SEX, -PHENOTYPE)%>%rename(id=FID)%>%filter(id%in%samples_keep)

  pheno_raw = read.table(paste0(path, '../data/icd10_body_self-report_blood_clean.txt'), header=1)
  pheno = pheno_raw%>%mutate(id=FID, traits=paste0(icd, ';', selfreport))%>%filter(id%in%samples_keep)
  pheno = pheno%>%select(-contains('blood'), -icd, -selfreport)%>%merge(pc, by='id') # we drop blood trait now
  
  # # self-define traits
  traits_add = c()
  # traits_add = c('pad', 'cataract', 'ihd')
  # if (F){
  #   # cataract
  #   cataract_code = c('H25', 'H26', 'H28')
  #   pheno$cataract = grepl('H25|H26|H28', pheno$traits)
  #   # pad
  #   pad_code = c('I7020','I7021','I7080','I7081','I7090','I7091','I721','I724','I730',
  #                'I731','I738','I739', 'I742','I743','I744','I745','I748','I749','I770',
  #                'I776','I778','I779','I792')
  #   pheno$pad = grepl(paste0(pad_code, collapse='|'), pheno$traits)
  #   # ihd, Ischaemic heart diseases
  #   ihd_code = paste0('I', 20:25)
  #   pheno$ihd = grepl(paste0(ihd_code, collapse='|'), pheno$traits)
  # }
  # pheno[,c('icd', 'selfreport')] = NULL
  
  traits_numeric = c('height', 'sitting_height', 'weight', 'bmi', 'hc', 'wc')
  
  # all traits
  traits = c(info1%>%filter(prop_case_b>0.001)%>%pull(coding), traits_numeric, traits_add)
  save(pheno, geno, info1, traits, file=paste0(path, 'data/for_phewas.rdata'))
}

load(paste0(path, 'data/for_phewas.rdata'))
load(paste0(path, 'data/snp_list.rdata'))
covar = pheno[,c('id', 'sex', 'age', 'center', 'batch', paste0('pc', 1:10))]
covar[,c('center', 'batch')] = sapply(covar[,c('center', 'batch')], as.character)

# E Endocrine, nutritional and metabolic diseases
# F Mental and behavioural disorders
# G Diseases of the nervous system
# H Diseases of the ear and mastoid process
# I Diseases of the circulatory system
# M Diseases of the musculoskeletal system and connective tissue
# Q Congenital malformations, deformations and chromosomal abnormalities

traits1 = traits[grepl('F|G|H|I', traits)]
traits1 = traits1[nchar(traits1)==3]
traits1 = traits1[i]

# sum(grepl('F', traits1)); sum(grepl('G', traits1)); sum(grepl('H', traits1)); sum(grepl('I', traits1))
# file_code = paste0(path, '../para/phewas_pheno_code.csv')
# map_pheno = read.csv(file_code)
# map_pheno%>%filter(coding%in%traits1)%>%pull(meaning)

for (trait in traits1){
  file_out = paste0(path, 'result/phewas/', trait, '.csv')
  if (!file.exists(file_out)){
    cols = c("phenotype", "snp", "adjustment", "beta", "SE", "OR", "p", "type", "n_total", "n_cases", "n_controls", "HWE_p","allele_freq", "n_no_snp", 
             "note", "model", "sex")
    cat(cols, sep=',', file=file_out, append=F)
    cat('\n', file=file_out, append=T)
  }
  if (trait %in% info1[,1]){
    pheno$y =  grepl(trait, pheno$traits)
  } else { # numeric traits
    pheno$y = pheno[,trait]
  }
  # sexs = c('male', 'female', 'both')
  sexs = c('both')
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
      if (!strsplit(snp, '_')[[1]][1]%in%snp_list[['FDXR']]) {next}
      print(snp)
      genotypes_raw = geno[, c('id', snp)]
      for (model in c('additive')){
        print(model)
        if (sum(grepl(snp, readLines(file_out))& 
                grepl(sex, readLines(file_out))& 
                grepl(model, readLines(file_out)))>0) {next}
        genotypes = genotypes_raw
        if (model == 'recessive'){
          genotypes[!is.na(genotypes[,2]) & genotypes[,2]==1, 2] = 0; genotypes[!is.na(genotypes[,2]) & genotypes[,2]==2, 2] = 1
        }else if (model == 'dominant'){
          genotypes[!is.na(genotypes[,2]) & genotypes[,2]==1, 2] = 1; genotypes[!is.na(genotypes[,2]) & genotypes[,2]==2, 2] = 1
        }
        covariates = covar
        if (sex!='both_sex'){
          covariates$sex = NULL
        } 
        # # interact term
        # if (sex=='both_sex'){ 
        #   covariates = covariates%>%merge(genotypes)
        #   covariates[,snp] = covariates[,snp]*covariates$sex
        #   covariates = covariates%>%rename(geno_sex=paste0(snp))
        # }
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

# 121
for i in {61..80}  
do
node=5
sbatch -N1 -n1 -c1 -w compute$node phewas.sh $i
done


## collect
library(dplyr)
library(stringr)
path = '/home/yanglab_data/user/zhanghy/gwas/ukb/phewas/'

load(paste0(path, 'data/snp_list.rdata'))
map_snp = data.frame()
for (i in names(snp_list)){
  map_snp = rbind(map_snp, data.frame(snp=snp_list[[i]], gene=i))
}

file_code = paste0(path, '../para/phewas_pheno_code.csv')
map_pheno = read.csv(file_code)%>%rename(phenotype=coding)

res = data.frame()
for (file in list.files(paste0(path, 'result/phewas/age60/'))){
  df = read.csv(paste0(path, 'result/phewas/age60/', file))
  res = rbind(res, df)
}

res$snp = sapply(res$snp, function(x){strsplit(x, '_')[[1]][1]})
res = res%>%filter(sex=='both')%>%merge(map_pheno)%>%merge(map_snp)

t = res%>%filter(p<0.005)
unique(res$gene); unique(t$gene)

write.csv(t, paste0(path, 'result/t60.csv'), row.names=F)

#=====================================================================================
# age>60
#=====================================================================================
# phewas60.r
library(PheWAS)
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/gwas/ukb/phewas/'

# args = commandArgs(T)
# i = as.numeric(args[1])

load(paste0(path, 'data/for_phewas.rdata'))
load(paste0(path, 'data/snp_list.rdata'))
pheno = pheno[pheno$age>60,]
covar = pheno[,c('id', 'sex', 'age', 'center', 'batch', paste0('pc', 1:10))]
covar[,c('center', 'batch')] = sapply(covar[,c('center', 'batch')], as.character)

traits1 = traits[grepl('F|G|H|I|E|M|Q', traits)]
traits1 = traits1[grepl('Q', traits1)]
traits1 = traits1[nchar(traits1)==3]
# traits1 = traits1[i]

for (trait in traits1){
  print(trait)
  file_out = paste0(path, 'result/phewas/age60/', trait, '.csv')
  if (!file.exists(file_out)){
    cols = c("phenotype", "snp", "adjustment", "beta", "SE", "OR", "p", "type", "n_total", "n_cases", "n_controls", "HWE_p","allele_freq", "n_no_snp", 
             "note", "model", "sex")
    cat(cols, sep=',', file=file_out, append=F)
    cat('\n', file=file_out, append=T)
  }
  if (trait %in% info1[,1]){pheno$y =  grepl(trait, pheno$traits)
  } else {pheno$y = pheno[,trait]}
  sexs = c('both')
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
      if (!strsplit(snp, '_')[[1]][1]%in%snp_list[['FDXR']]) {next}
      print(snp)
      genotypes_raw = geno[, c('id', snp)]
      for (model in c('additive')){
        if (sum(grepl(snp, readLines(file_out))&grepl(sex, readLines(file_out))&grepl(model, readLines(file_out)))>0) {next}
        genotypes = genotypes_raw
        if (model == 'recessive'){
          genotypes[!is.na(genotypes[,2]) & genotypes[,2]==1, 2] = 0; genotypes[!is.na(genotypes[,2]) & genotypes[,2]==2, 2] = 1
        }else if (model == 'dominant'){
          genotypes[!is.na(genotypes[,2]) & genotypes[,2]==1, 2] = 1; genotypes[!is.na(genotypes[,2]) & genotypes[,2]==2, 2] = 1
        }
        covariates = covar
        if (sex!='both_sex'){
          covariates$sex = NULL
        } 
        out = phewas(phenotypes=phenotypes, genotypes=genotypes,covariates=covar, additive.genotypes = T) 
        out$model = model
        out$phenotype = trait
        out$sex = sex
        cat(unlist(out[1,]), sep=',', file=file_out, append=T)
        cat('\n', file=file_out, append=T)
      }
    }
  }
}

# 147
for i in {131..147}  
do
node=3
sbatch -N1 -n1 -c1 -w compute$node phewas60.sh $i
done



## old
file_bim = '/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt_qc/qc_1_bychr/chr17_1.bim'
file_eqtl = '/home/yanglab_data/user/zhanghy/gwas/summary/eqtl/clean/GTEx_v8_FDXR.txt'
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
