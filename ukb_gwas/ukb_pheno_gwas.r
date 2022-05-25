## ukb_pheno_gwas
#=====================================================================================
# extract base character and icd col 
# 21022 age
# 31 sex
# 21000 ethnicity  
# 22000 batch
# 54 center 
# 26412 employment score https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=26412
# 26414 education score https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=26414
# 738 income before tax https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=738

# 20001	Cancer code, self-reported
# 20002	Non-cancer illness code, self-reported
# 41270 icd-10

# 21001 bmi
# 12144	Height | only has 2 and 3
# 21002	Weight
# 48	Waist circumference
# 49	Hip circumference
# 20015	Sitting height
#=====================================================================================
## get field_col
library(dplyr)
all_field = as.vector(unlist(read.table('/home/yanglab_data/pub/DownloadUKB/Converted/Main1.tab', nrows=1)))
code1 = read.csv('/home/yanglab_data/user/zhanghy/gwas/ukb/info/Data_Dictionary_Showcase.csv')
code2 = read.csv('/home/yanglab_data/user/zhanghy/gwas/ukb/info/Codings_Showcase.csv')

## diet
t = code1%>%filter(grepl('Diet', code1$Path))
keep_field = c('f.eid2', paste0('f.', t$FieldID, '.0.0'))
n_field = unlist(sapply(keep_field, function(x){which(all_field==x)}))
paste0(n_field, collapse=',')

for (i in unique(t$Coding)){
  print(i)
  print(code2%>%filter(Coding==i))
}
105030

## base_icd10
keep_field = c('f.eid2', 'f.21022.0.0', 'f.31.0.0', 'f.21000.0.0', 'f.22000.0.0', 'f.54.0.0', 'f.26412.0.0', 'f.26414.0.0', 'f.738.0.0')
n_field = sapply(keep_field, function(x){which(all_field==x)})
n_field = c(n_field,  which(grepl('41270.0', all_field)))
paste0(n_field, collapse=',')

## edu, employ, income
keep_field = c('f.eid2', 'f.26412.0.0', 'f.26414.0.0', 'f.738.0.0')
n_field = sapply(keep_field, function(x){which(all_field==x)})
n_field = c(n_field,  which(grepl('41270.0', all_field)))
paste0(n_field, collapse=',')

## body_selfreported-disease
keep_field = c('f.eid2', 'f.21001.0.0', 'f.12144.2.0', 'f.21002.0.0', 'f.48.0.0', 'f.49.0.0', 'f.20015.0.0')
n_field = sapply(keep_field, function(x){which(all_field==x)})
n_field = c(n_field, which(grepl('20001.0', all_field)), which(grepl('20002.0', all_field)))
paste0(n_field, collapse=',')

## age cataract diagnosed
keep_field = c('f.eid2', 'f.4700.0.0')
n_field = sapply(keep_field, function(x){which(all_field==x)})
paste0(n_field, collapse=',')

## blood
# Category 2403: Blood, blood-forming organs and certain immune disorders | missing
# Category 17518: Blood biochemistry
# Category 100081: Blood count
code = read.csv('/home/yanglab_data/user/zhanghy/gwas/ukb/info/Data_Dictionary_Showcase.csv')[,c(1:4,8)]
t = code[code$Category%in%c(17518, 100081),]

keep_field = code[code$Category%in%c(17518, 100081), 'FieldID']
keep_field = c('f.eid2', paste0('f.', keep_field, '.0.0'))

n_field = sapply(keep_field, function(x){which(all_field==x)})
n_field
paste0(n_field, collapse=',')

## extract
cat /home/yanglab_data/pub/DownloadUKB/Converted/Main1.tab | cut -f 1,9833,24,9782,9992,95,15741-15953 \
> /home/yanglab_data/user/zhanghy/gwas/ukb/data/icd10.tab

cat /home/yanglab_data/pub/DownloadUKB/Converted/Main1.tab | cut -f 1,9785,6461,9789,74,78,8053,6708-6713,\
6732-6765 > /home/yanglab_data/user/zhanghy/gwas/ukb/data/body_selfreport.tab

cat /home/yanglab_data/pub/DownloadUKB/Converted/Main1.tab | cut -f 1,13764,13777,13790,13803,13816,13829,13842,13855,13868,13881,13894,13907,13920,13933,13946,13959,13972,13985,13998,14011,14024,14037,14050,14063,14076,14089,14102,14115,14128,14141,14154,14222,14236,14250,14264,14278,14292,14306,14320,14334,14348,14362,14376,14390,14404,14418,14432,14444,14458,14472,14486,14500,14514,14528,14542,14556,14570,14584,14598,14612,14626 \
> /home/yanglab_data/user/zhanghy/gwas/ukb/data/blood.tab

cat /home/yanglab_data/pub/DownloadUKB/Converted/Main1.tab | cut -f 1,13742,13744,485 \
> /home/yanglab_data/user/zhanghy/gwas/ukb/data/edu_etc.tab

cat /home/yanglab_data/pub/DownloadUKB/Converted/Main1.tab | cut -f 1,3021 \
> /home/yanglab_data/user/zhanghy/gwas/ukb/data/cataract_age.tab

cat /home/yanglab_data/pub/DownloadUKB/Converted/Main1.tab | cut -f 1,672,676,680,684,688,692,696,700,704,708,712,716,720,724,728,732,736,740,744,748,752,756,760,764,768,1083,1437,5913,6426,6427,6439,6449,8110,8111,8116,8121,8126,8131,8136,8141,8246,8266,8296,8361,8386,8411,8436,8456,8481,8496,8506,8516,8531,8541,8551,8561,8571,8581,8591,8601,8651,8676,16757,16762,16767,16772,16777,16782,16787,16792,16797,16802,16807,16812,16817,16822,16827,16832,16837,16842,16847,16852,16857,16862,16867,16872,16877,16882,16887,16892,16897,16902,16907,16912,16917,16922,16927,16932,16937,16942,16947,16952,16957,16962,16967,16972,16977,16982,16987,16992,16997,17002,17007,17012,17017,17022,17027,17032,17037,17042,17047,17052,17057,17062,17067,17072,17077,17082,17087,17092,17097,17102,17107,17112,17117,17122,17127,17132,17137,17142,17147,17152,17157,17162,17167,17172,17177,17182,17187,17192,17197,17202,17207,17212,17217,17222,17227,17232,17237,17242,17247,17252,17257,17262,17267,17272,17277,17282,17287,17292,17297,17302,17307,17312,17317,17322,17327,17332,17337,17342,17347,17352,17357,17362,17367,17372,17377,17382,17387,17392,17397,17402,17407,17412,17417,17422,17427,17432,17437,17442,17447,17452,17457,17462,17467,17472,17477,17482,17487,17492,17497,17502,17507,17512,17517,17522,17527,17532,17537,17542,17547,17552,17557,17562,17567,17572,17577,17582,17587,17592,17597,17602,17607,17612,17617,17622,17627,17632,17637,17642,17647,17652,17657,17662,17667,17672,17677,17682,17687,17692,17697,17702,17707,17712,17717,17722,17727,17732,17737,17742,17747,17752,17757,17762,17767,17772,17777,17782,17787,17792,17797,17802,17807,17812,17817,17822,17827,17832,17837,17842,17847,17852,17857,17862,17867,17872,17877,17882,17887,17892,17897,17902,17907,17912,17917,17922,17927,17932,17937,17942,17947,17952,17957,17962,17967,17972,17977,17982,17987,17992,17997,18002,18007,18012,18017,18022,18027,18032,18037,18042,18047,18052,18057,18062,18067,18072,18077,18082,18087,18092,18097,18102,18107,18112,18117,18122,18127,18132,18137,18142 \
> /home/yanglab_data/user/zhanghy/gwas/ukb/data/diet.tab
#=====================================================================================
# clean icd-10
#=====================================================================================
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/gwas/ukb/data/'

df_raw = read.table(paste0(path, 'icd10.tab'), header=1)
df = df_raw
icd_col = colnames(df)[grepl('41270', colnames(df))]

df = cbind(df$f.eid2, df)
colnames(df)[1:7] = c('FID', 'IID', 'sex', 'center', 'ethnicity', 'age', 'batch')
df = df[,c('FID', 'IID', 'sex', 'age', 'ethnicity', 'center', 'batch', colnames(df)[8:ncol(df)])]
# drop na 
df = df[rowSums(is.na(df[, 1:7]))==0,]

# collapse icd
df$icd = NA
for (i in 1:nrow(df)){
  if (i%%500==0){
    print(paste0(i, '/', nrow(df)))
  }
  t = df[i, icd_col]
  t = t[!is.na(t)]
  df[i, 'icd'] = paste(t, collapse=";")
}
df[,icd_col] = NULL


# add NONE on empty cell, nor bolt will have parse error
df[df$icd=='', 'icd'] = 'NONE'

write.table(df, paste0(path, 'icd10_clean.txt'), row.names=F, sep='\t', quote=F)
#=====================================================================================
# clean body_self-report
# 20001	Cancer code, self-report
# 20002	Non-cancer illness code, self-report
#=====================================================================================
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/gwas/ukb/data/'

df_icd10 = read.table(paste0(path, 'icd10_clean.txt'), header=1)

# body_selfreport
df = read.table(paste0(path, 'body_selfreport.tab'), header=1)
df = df[df$f.eid2%in%df_icd10$IID,] # drop sample with missing base info 
cancer_col = colnames(df)[grepl('20001', colnames(df))]
noncancer_col = colnames(df)[grepl('20002', colnames(df))]

df[,cancer_col] = sapply(df[,cancer_col], function(x){paste0('c_', x, '_c')})
df[,noncancer_col] = sapply(df[,noncancer_col], function(x){paste0('nc_', x, '_nc')})
df[df=='nc_NA_nc'|df=='c_NA_c'] = NA

df$selfreport = NA

for (i in 1:nrow(df)){
  if (i%%500==0){
    print(paste0(i, '/', nrow(df)))
  }
  traits = df[i, c(cancer_col, noncancer_col)]
  traits = traits[!is.na(traits)]
  if (length(traits)==0){
    traits = 'NONE'
  } else{
    traits = paste(traits, collapse=";")
  }
  df[i, 'selfreport'] = traits
}

df = df%>%rename(FID=f.eid2, wc=f.48.0.0, hc=f.49.0.0, height=f.12144.2.0, 
                 weight=f.21002.0.0,sitting_height=f.20015.0.0, bmi=f.21001.0.0)
df = df[,c('FID', 'selfreport', 'height', 'sitting_height', 'weight', 'bmi', 'hc', 'wc')]

res = df_icd10%>%merge(df, by='FID')

# blood
df = read.table(paste0(path, 'blood.tab'), header=1)
names(df) = c('FID', gsub('.0.0', '', gsub('f.', 'blood_', names(df)[2:ncol(df)]), fixed=T))
res = res%>%merge(df, by='FID')

write.table(res, paste0(path, 'icd10_body_self-report_blood_clean.txt'), row.names=F, sep='\t', quote=F)
#=====================================================================================
# description of trait and nsample
#=====================================================================================
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/gwas/ukb/data/'
df = read.table(paste0(path, 'icd10_body_self-report_blood_clean.txt'), header=1)

## select sample
file_fam = '/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt.fam'
file_out = '/home/yanglab_data/user/zhanghy/gwas/ukb/data/brt_pheno_freq.csv'

samples_keep = read.table(file_fam)[,2]
df = df[df$FID%in%samples_keep,]
df$traits = paste0(df$icd, ';', df$selfreport)

# get phenofreq, done!
if (F){ 
  traits_sr = c() # sr = selfreport
  traits_icd = c()
  for (i in 1:nrow(df)){
    if (i%%500==0){
      print(paste0(i, '/', nrow(df)))
    }
    traits_sr = c(traits_sr, strsplit(df[i, 'selfreport'], ';')[[1]])
    traits_sr = unique(traits_sr)
    traits_icd = c(traits_icd, strsplit(df[i, 'icd'], ';')[[1]])
    traits_icd = unique(traits_icd)
  }
  
  traits_icd = c(traits_icd, sapply(traits_icd, function(x){substr(x, 1, 3)[1]})) # add main group
  traits = c(traits_icd, traits_sr)
  traits = unique(traits)
  traits = traits[!traits%in%c('NONE', 'NON')]
  
  save(traits, file=paste0(path, 'brt_all_traits.rdata'))
  
  load(paste0(path, 'brt_all_traits.rdata'))
  
  cols = c('trait', 'n_case_m', 'n_ctrl_m', 'prop_case_m', 'n_case_f',
           'n_ctrl_f', 'prop_case_f', 'n_case_b', 'n_ctrl_b', 'prop_case_b')
  cat(cols, sep=',', file=file_out, append=F)
  cat('\n', file=file_out, append=T)
  
  for (trait in traits){
    if (sum(grepl(paste0(trait, ','), readLines(file_out)))>0) {next}
    
    print(which(traits==trait))
    n_case_m = sum(df$sex==1& grepl(trait, df$traits))
    n_ctrl_m = sum(df$sex==1& !grepl(trait, df$traits))
    prop_case_m = n_case_m/(n_case_m+n_ctrl_m)
    
    n_case_f = sum(df$sex==0& grepl(trait, df$traits))
    n_ctrl_f = sum(df$sex==0& !grepl(trait, df$traits))
    prop_case_f = n_case_f/(n_case_f+n_ctrl_f)
    
    n_case_b = sum(grepl(trait, df$traits))
    n_ctrl_b = sum(!grepl(trait, df$traits))
    prop_case_b = n_case_b/(n_case_b+n_ctrl_b)
    
    out = c(trait, n_case_m, n_ctrl_m, prop_case_m, n_case_f,
            n_ctrl_f, prop_case_f, n_case_b, n_ctrl_b, prop_case_b)
    
    cat(out, sep=',', file=file_out, append=T)
    cat('\n', file=file_out, append=T)
  }
}

res = read.csv(file_out)

sum(res$n_case_m+res$n_case_f != res$n_case_b)
sum(res$n_ctrl_m+res$n_ctrl_f != res$n_ctrl_b)

# # combine coding
# code1 = read.csv(paste0(path, '../info/Data_Dictionary_Showcase.csv'))
# code2 = read.csv(paste0(path, '../info/Codings_Showcase.csv'))
# self_cancer_code = code2[code2$Coding==code1[code1$FieldID=='20001', 'Coding'],][,2:3]
# self_noncancer_code = code2[code2$Coding==code1[code1$FieldID=='20002', 'Coding'],][,2:3]
# code1[code1$Category%in%c(17518, 100081), c('FieldID', 'Field')]
# blood_code = code1[code1$Category%in%c(17518, 100081), c('FieldID', 'Field')]
# names(blood_code) = c('coding', 'meaning')
# 
# self_cancer_code[,1] = paste0('c_', self_cancer_code[,1], '_c')
# self_noncancer_code[,1] = paste0('nc_', self_noncancer_code[,1], '_nc')
# blood_code[,1] = paste0('blood_', blood_code[,1])
# 
# icd_code = read.csv(paste0(path, '../info/icd10_all.csv'))[,1:2]
# names(self_cancer_code) = names(self_noncancer_code) = names(icd_code)
# 
# code = rbind(icd_code, self_noncancer_code, self_cancer_code, blood_code)
# code[code$coding=='F171', 2] = 'F17.1 Harmful use (tobacoo)'
# code[code$coding=='F101', 2] = 'F10.1 Harmful use (alcohol)'
# 
# write.csv(code, file_code, row.names=F)

file_code = paste0(path, '../para/phewas_pheno_code.csv')
code = read.csv(file_code)

# add pheno name
names(res)[1] = 'coding'
res[!res$coding%in%code[,1],]
res = code%>%merge(res, by='coding')
res = rbind(res[!grepl('_c', res[,1]),], res[grepl('_c', res[,1]),])

res1 = res[res$prop_case_m>0.001|res$prop_case_f>0.001,]

file_out1 = '/home/yanglab_data/user/zhanghy/gwas/ukb/data/brt_pheno_freq>0.001.csv'
write.csv(res1, file_out1, row.names=F)

#=====================================================================================
# run
#=====================================================================================
## make 1,0 pheno 
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/gwas/ukb/'
file_info = paste0(path, '/data/brt_pheno_freq>0.01.csv')
file_pheno = paste0(path, '/data/icd10_body_self-report_blood_clean.txt')
file_fam = paste0(path,'/bfile/brt.fam')
file_out1 = paste0(path, '/data/brt_pheno_for_bolt.txt')
file_out2 = paste0(path, '/data/brt_m_pheno_for_bolt.txt')
file_out3 = paste0(path, '/data/brt_f_pheno_for_bolt.txt')
file_code = paste0(path, 'para/phewas_pheno_code.csv')
file_cat_age = paste0(path, 'data/cataract_age.tab')

code = read.csv(file_code)
info = read.csv(file_info)

df = read.table(file_pheno, header=1)
samples_keep = read.table(file_fam)[,2]
df = df[df$FID%in%samples_keep,]
df$traits = paste0(df$icd, ';', df$selfreport)
df = df[,c(names(df)[1:7],'traits', 'sitting_height')]

## make 0,1 phenotypes
traits = info[info$prop_case_b>0.01, 'coding']
traits = traits[nchar(traits)==3]
traits = c(traits, 'E10') #E10=T1D 

for (trait in traits){
  df[,trait] = as.numeric(grepl(trait, df$traits))
}

## cataract
cataract_code = c('H25', 'H26', 'H28')
df$cataract = 0
df[grepl('H25|H26|H28', df$traits), 'cataract'] = 1
# cataract age
df_add = read.table(file_cat_age, header=1)%>%rename(FID=f.eid2, cataract_age=f.4700.0.0)
df = df%>%merge(df_add, by='FID')

# pad
pad_code = c('I7020','I7021','I7080','I7081','I7090','I7091','I721','I724','I730',
             'I731','I738','I739', 'I742','I743','I744','I745','I748','I749','I770',
             'I776','I778','I779','I792')
df$pad = 0
df[grepl(paste0(pad_code, collapse='|'), df$traits), 'pad'] = 1

df$traits = NULL

write.table(df, file_out1, quote=F, row.names=F)
df2 = df[df$sex==1,]; df3 = df[df$sex==0,]
write.table(df2, file_out2, quote=F, row.names=F)
write.table(df3, file_out3, quote=F, row.names=F)

# 80%sample
file_out4 = paste0(path, '/data/brt_m_0.8_pheno_for_bolt.txt')
file_out5 = paste0(path, '/data/brt_f_0.8_pheno_for_bolt.txt')
set.seed(1)
df4 = df2[sample(1:nrow(df2), 0.8*nrow(df2)),]
df5 = df3[sample(1:nrow(df3), 0.8*nrow(df3)),]
write.table(df4, file_out4, quote=F, row.names=F)
write.table(df5, file_out5, quote=F, row.names=F)


# cat-age | ~5000, limited sample size
df1 = df%>%filter(cataract_age>=40)


## bolt.sh
#!/usr/bin/bash
trait="$1"
nthread="$2"

file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt"
file_out="/home/yanglab_data/user/zhanghy/gwas/ukb/code/$trait.txt"
file_ld="/home/yanglab_data/user/zhanghy/soft_slurm/bolt/tables/LDSCORE.1000G_EUR.tab.gz"
file_pheno="/home/yanglab_data/user/zhanghy/gwas/ukb/data/brt_pheno_for_bolt.txt"
file_modsnp="/home/yanglab_data/user/zhanghy/soft_slurm/bolt/hm3/ukbEURu_imp_v2_HM3_QC_R2_09.snplist"

/home/yanglab_data/user/zhanghy/soft_slurm/bolt/BOLT-LMM_v2.3.5/bolt \
--bfile=$file_bfile \
--LDscoresFile=$file_ld \
--lmm \
--phenoFile=$file_pheno \
--phenoCol=$trait \
--covarFile=$file_pheno \
--qCovarCol=age \
--covarCol=sex \
--covarCol=center \
--covarCol=batch \
--covarMaxLevels=120 \
--modelSnps=$file_modsnp \
--numThreads=$nthread \
--statsFile=$file_out

## bolt_sex80.sh
#!/usr/bin/bash
trait="$1"
nthread="$2"
sex="$3"

file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt"
file_out="/home/yanglab_data/user/zhanghy/gwas/ukb/code/${trait}_${sex}_0.8.txt"
file_ld="/home/yanglab_data/user/zhanghy/soft_slurm/bolt/tables/LDSCORE.1000G_EUR.tab.gz"
file_pheno="/home/yanglab_data/user/zhanghy/gwas/ukb/data/brt_${sex}_0.8_pheno_for_bolt.txt"
file_modsnp="/home/yanglab_data/user/zhanghy/soft_slurm/bolt/hm3/ukbEURu_imp_v2_HM3_QC_R2_09.snplist"

/home/yanglab_data/user/zhanghy/soft_slurm/bolt/BOLT-LMM_v2.3.5/bolt \
--bfile=$file_bfile \
--LDscoresFile=$file_ld \
--lmm \
--phenoFile=$file_pheno \
--phenoCol=$trait \
--covarFile=$file_pheno \
--qCovarCol=age \
--covarCol=center \
--covarCol=batch \
--covarMaxLevels=120 \
--modelSnps=$file_modsnp \
--numThreads=$nthread \
--statsFile=$file_out


# phe
[1] "A09" "A41" "B95" "B96" "C44" "C50" "C61" "C77" "C78" "C79" "D12" "D17"
[13] "D22" "D23" "D25" "D50" "D64" "E03" "E11" "E66" "E78" "E87" "F10" "F17"
[25] "F32" "F41" "G40" "G47" "G55" "G56" "H02" "H25" "H26" "H35" "H40" "I10"
[37] "I20" "I21" "I25" "I26" "I44" "I48" "I50" "I51" "I63" "I73" "I80" "I83"
[49] "I84" "I95" "J18" "J22" "J34" "J44" "J45" "J90" "J98" "K20" "K21" "K22"
[61] "K25" "K29" "K30" "K31" "K40" "K42" "K44" "K52" "K56" "K57" "K58" "K59"
[73] "K62" "K63" "K66" "K76" "K80" "K92" "L03" "L72" "L98" "M06" "M10" "M13"
[85] "M15" "M16" "M17" "M19" "M20" "M23" "M25" "M47" "M48" "M51" "M54" "M65"
[97] "M67" "M72" "M75" "M79" "M81" "N17" "N18" "N20" "N32" "N39" "N40" "N81"
[109] "N83" "N84" "N85" "N92" "N93" "N95" "O70" "R00" "R04" "R06" "R07" "R10"
[121] "R11" "R13" "R19" "R29" "R31" "R33" "R35" "R39" "R41" "R42" "R50" "R51"
[133] "R55" "R63" "R69" "R79" "R94" "S52" "S82" "T81" "T84" "W01" "W19" "X59"
[145] "Y83" "Z03" "Z08" "Z09" "Z12" "Z13" "Z30" "Z37" "Z46" "Z50" "Z51" "Z53"
[157] "Z72" "Z80" "Z82" "Z85" "Z86" "Z87" "Z88" "Z90" "Z91" "Z92" "Z95" "Z96"
[169] "Z98"


## run
# both sex
trait=K92
nthread=50
node=10
file_log="/home/yanglab_data/user/zhanghy/gwas/ukb/code/$trait.log"
sbatch -N1 -n1 -c$nthread -o $file_log -w compute$node bolt.sh $trait $nthread

# single sex
trait=cataract
sex=f
nthread=20
node=7
# file_log="/home/yanglab_data/user/zhanghy/gwas/ukb/code/${trait}_${sex}.log"
file_log="/home/yanglab_data/user/zhanghy/gwas/ukb/code/${trait}-adjE11_${sex}_0.8.log"
# file_log="/home/yanglab_data/user/zhanghy/gwas/ukb/code/${trait}_${sex}_0.8.log"
sbatch -N1 -n1 -c$nthread -o $file_log -w compute$node bolt_sex80.sh $trait $nthread $sex
# sbatch -N1 -n1 -c$nthread -o $file_log -w compute$node bolt_sex80_noadj.sh $trait $nthread $sex


# 2, 5, 6, 7 can be used