
#=====================================================================================
# prs: clump
#=====================================================================================
pop='eas'
path="/home/yanglab_data/user/zhanghy/gwas/summary/microbiome/${pop}/"
if [ $pop == "eas" ]
then file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/chn"
else file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt"
fi

datas=(`ls ${path}result/prs/clump/for_clump`)
start=0; end="${#datas[*]}"
for ((i=$start; i<$end; i++))
do
data=${datas[$i]/.txt/}
echo $data
file_gwas="${path}result/prs/clump/for_clump/$data.txt"
file_out1="${path}result/prs/clump/clump_out/$data"
file_out2="${path}result/prs/clump/clump_out/snp/$data"
file_out3="${path}result/prs/clump/clump_out/pval/$data"

if [ ! -f ${file_out1}.clumped ]
then
plink \
    --bfile $file_bfile \
    --clump $file_gwas \
    --clump-p1 1 --clump-r2 0.05 --clump-kb 1000 \
    --clump-snp-field SNP --clump-field P \
    --out $file_out1
awk 'NR!=1{print $3}' ${file_out1}.clumped >  $file_out2
awk '{print $1,$4}' $file_gwas > $file_out3
fi done

#=====================================================================================
# prs: main
#=====================================================================================
pop='eas'
path="/home/yanglab_data/user/zhanghy/gwas/summary/microbiome/${pop}/result/prs/"

if [ $pop == "eas" ]
then file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/chn"
else file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt"
fi

file_range='${path}range_list'
echo "0.00000001 0 0.00000001" > $file_range 
echo "0.00001 0 0.00001" >> $file_range 
echo "0.001 0 0.001" >> $file_range 
echo "0.01 0 0.01" >> $file_range 
echo "0.05 0 0.05" >> $file_range
echo "0.1 0 0.1" >> $file_range
echo "0.5 0 0.5" >> $file_range

datas=(`ls ${path}clump/clump_out | grep clumped`)
start=0; end="${#datas[*]}"
for ((i=$start; i<$end; i++))
do
data=${datas[$i]/.clumped/}
file_gwas="${path}clump/for_clump/${data}.txt"
file_out="${path}prs_out/$data"
echo $data
if [ ! -f ${file_out}.0.1.profile ]
then
plink \
  --bfile $file_bfile \
  --score $file_gwas 1 2 3 header \
  --q-score-range $file_range "${path}clump/clump_out/pval/${data}" \
  --extract "${path}clump/clump_out/snp/${data}" \
  --out $file_out
fi done


#=====================================================================================
# prs reg: icd10
#=====================================================================================
source('/home/yanglab_data/user/zhanghy/project/temp/code/microbio/source_microbio.r')
pop = 'eas'
path_prs = paste0(path_list[[pop]][['data']], 'result/prs/prs_out/')
df_raw = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/data/icd10_clean.txt', header=1)

## traits
traits = paste(df$icd, collapse=";") 
traits = unique(strsplit(traits, ';')[[1]])
traits = traits[traits!='NONE']
# prs
datas = list.files(path_prs)
datas = unique(sapply(datas, function(x){strsplit(x, '.', fixed=T)[[1]][1]}))

for (data in datas){
  print(data)
  file_out = paste0(path_prs, '../reg/icd/', data, '.rdata')
  if (file.exists(file_out)){next}
  out = c()
  for (thres in c('0.5', '0.1', '0.05', '0.001', '0.00001', '0.00000001')){
    print(thres)
    file_prs = paste0(path_prs, data, '.', thres, '.profile')
    if (!file.exists(file_prs)){next}
    prs = read.table(file_prs, header=1)%>%select(FID, SCORE)
    df = df_raw%>%select(-IID, icd)%>%merge(prs, 'FID')%>%mutate(batch=as.factor(batch), center=as.factor(center))
    for (trait in traits){
      df$y = as.numeric(grepl(trait, df$icd))
      if (length(unique(df$y))==1){next}
      if (table(df$y)[2] < (nrow(df)*0.01)){next}
      mod_full <- glm(y ~ SCORE + sex + age + batch + center , data = df, family = binomial())
      mod_base <- glm(y ~ age + sex + batch + center, data = df, family = binomial())
      LLf = logLik(mod_full)[1] 
      LLc = logLik(mod_base)[1]   
      out = c(out, trait, data, thres, summary(mod_full)$coefficients[2, c(1, 2, 4)], LLf, LLc)
    }
  }
  save(out, file=file_out)
}

### collect
collect= c()
for (data in datas[1:5]){
  print(data)
  file_out = paste0(path_prs, '../reg/icd/', data, '.rdata')
  if (!file.exists(file_out)){next}
  load(file_out)
  collect = c(collect, out)
}

res = data.frame(matrix(collect, ncol=8, byrow=T))
names(res) = c('coding', 'prs_trait', 'prs_thres', 'b', 'se', 'p', 'llf', 'llc')
res = res%>%filter(p<0.001)

## add icd label
file_code = paste0('/home/yanglab_data/user/zhanghy/gwas/ukb/para/phewas_pheno_code.csv')
code = read.csv(file_code)
res = res%>%merge(code, by='coding')

res = res%>%filter(!grepl('perineal laceration', 'meaning'))

write.csv(res, paste0(path_prs, '../../collect/', pop, '_prs', '.csv'), row.names = F)


#=====================================================================================
# prs reg: diet
#=====================================================================================
source('/home/yanglab_data/user/zhanghy/project/temp/code/microbio/source_microbio.r')
code1 = read.csv('/home/yanglab_data/user/zhanghy/gwas/ukb/info/Data_Dictionary_Showcase.csv')
code2 = read.csv('/home/yanglab_data/user/zhanghy/gwas/ukb/info/Codings_Showcase.csv')
pop = 'eas'
path_prs = paste0(path_list[[pop]][['data']], 'result/prs/prs_out/')
df_raw = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/data/icd10_clean.txt', header=1)
diet_raw = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/data/diet.tab', header=1, sep='\t')
keep_id = read.table('/home/yanglab_data/user/zhanghy/gwas/summary/microbiome/eas/result/prs/prs_out/c_Bacilli.0.1.profile', header=1)$FID

diet = diet_raw%>%select(-f.105030.0.0, -f.105010.0.0, -f.105010.0.0, -f.20078.0.0)%>%rename(FID=f.eid2)
df = df_raw%>%filter(FID%in%keep_id)%>%merge(diet, by='FID')%>%mutate(batch=as.factor(batch), center=as.factor(center))
df = df[,colSums(!is.na(df))>100] # colwith > 100 sample

# prs and diet vars
datas = list.files(path_prs)
datas = gsub('.log', '', datas[grepl('.log', datas)])
vars = names(df)[grep('f.', names(df))[1]:ncol(df)]

## get coding 
for (i in vars){
  field = gsub('f.', '', gsub('.0.0', '', i, fixed=T))
  code = code1%>%filter(FieldID==field)%>%pull(Coding)
  if (is.na(code)) {next} # numeric
  map = code2%>%filter(Coding==code)%>%select(-Coding)%>%mutate(Value=as.numeric(Value))
  for (j in 1:nrow(map)){ # merge lose row order
    df[!is.na(df[,i])&df[,i]==map[j, 1], i] = map[j, 2]
  }
}
## frq
for (i in vars){
  field = gsub('f.', '', gsub('.0.0', '', i, fixed=T))
  code = code1%>%filter(FieldID==field)%>%pull(Coding)
  if (is.na(code)) {next} # numeric
  print(i)
  print(table(df[,i]))
}

for (data in datas){
  print(data)
  file_out = paste0(path_prs, '../reg/diet/', data, '.rdata')
  if (file.exists(file_out)){next}
  out = data.frame()
  for (thres in c('0.5', '0.1', '0.05', '0.001', '0.00001', '0.00000001')){
    print(thres)
    file_prs = paste0(path_prs, data, '.', thres, '.profile')
    if (!file.exists(file_prs)){next}
    prs = read.table(file_prs, header=1)%>%select(FID, SCORE)
    sub = df%>%merge(prs, 'FID')
    for (var in vars){
      formula = paste0('SCORE ~ ', var, ' + sex + age')
      mod = lm(formula, data = sub)
      coef = summary(mod)$coefficients[,c(1,2,4)]; coef = coef[-1,]
      coef = data.frame(cbind(var, data, thres, gsub(var, '', rownames(coef)), coef))
      names(coef) = c('y_diet', 'x_prs', 'prs_thres', 'var', 'b', 'se', 'p')
      out = rbind(out, coef)
      rownames(out)=NULL
    }
  }
  save(out, file=file_out)
}

### collect
collect = data.frame()
for (data in datas){
  file_out = paste0(path_prs, '../reg/diet/', data, '.rdata')
  load(file_out)
  if (nrow(out)==0){print(data)} # check
  out$group = paste(out$y_diet, out$x_prs, out$prs_thres, sep=';')
  keep = unique(out%>%filter(!var%in%c('age', 'sex')&p<0.01)%>%pull(group)) # keep sig
  out = out%>%filter(group%in%keep)
  collect = rbind(collect, out)
}

t = collect%>%filter(p<0.05)


## add icd label


write.csv(res, paste0(path_prs, '../../collect/', pop, '_prs', '.csv'), row.names = F)

