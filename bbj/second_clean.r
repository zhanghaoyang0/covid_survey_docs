#=====================================================================================
# bbj second clean on R (the previous is on py)
# some data may filter with 1000g with mhc, e.g., Case_control_14
#=====================================================================================
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/'
file_bim = '/home/yanglab_data/user/zhanghy/gwas/bfile/1000g/no_mhc/eas/eas_nomhc.bim'

bim = read.table(file_bim)[, c(2, 1, 4)]
colnames(bim) = c('SNP', 'CHR', 'POS')

para = read.csv(paste0(path, 'para/para.csv'))
para[is.na(para$total), 'total'] = para[is.na(para$total), 'case'] + para[is.na(para$total), 'ctrl']

datas = para[!grepl('eur', para[,1]), 1]

for (data in datas){
  print(data)
  if (file.exists(paste0(path, 'clean/', data))){
    next()
  }
  df_raw = read.csv(paste0(path, 'raw/', data), sep = '\t')
  df = df_raw
  sub_para = para[para$id==data,]
  # replace colname
  if (is.na(sub_para[, 'FRQ'])){
    replace_cols = c('A1', 'A2', 'POS', 'SNP', 'P', 'FRQ_CASE', 'FRQ_CTRL', 'SE')
  } else {
    replace_cols = c('A1', 'A2', 'POS', 'SNP', 'P', 'FRQ', 'SE')
  }
  
  for (col in replace_cols){
    colnames(df)[colnames(df)==sub_para[,col]] = col
  }
  # for or, the se of case_control_105 not need to transform
  if (sub_para[, 'effect_type']=='OR'){
    df$BETA = log(df$OR)
  }
  # for not total frq
  if (is.na(sub_para[, 'FRQ'])){
    df$FRQ = (df$FRQ_CASE*sub_para$case+df$FRQ_CTRL*sub_para$ctrl)/sub_para$total
  } 
  # for N
  if (!'N' %in% colnames(df)){
    df$N = sub_para$total
  }
  
  df$SNP = NULL
  df = df%>%merge(bim, by=c('CHR', 'POS'))
  df = df[df$FRQ>0.01&df$FRQ<0.99,]
  df = na.omit(df)
  
  # remove duplicated snp
  df$order = 1:nrow(df)
  df = df[order(df$P),]
  df = df[!duplicated(df$SNP),]
  df = df[order(df$order),]
  df$order = NULL
  df = df[,c('CHR', 'POS', 'SNP', 'FRQ', 'A1', 'A2', 'BETA', 'SE', 'P', 'N')]
  
  write.table(df, paste0(path, 'clean/', data, '.txt'), sep = '\t', row.names = F, quote = F)
}

#=====================================================================================
# munge 
#=====================================================================================
path='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/'
datas=`ls $path/clean/ | grep txt`
datas=($datas)

len="${#datas[*]}"

for ((i=0; i<$len; i++))
  do
data=${datas[$i]}
echo $data
out=`echo $data | sed -e "s/.txt//g"`
file_out=$path/munge/$out

if [ ! -f $file_out.sumstats.gz ] 
then
/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/munge_sumstats.py \
--sumstats $path/clean/$data \
--out $path/munge/$out \
--merge-alleles /home/yanglab_data/user/zhanghy/soft_slurm/ldsc/w_hm3.snplist
fi
done

#=====================================================================================
# h2 
#=====================================================================================
## ldsc: cal h2
path="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/"
path_ld="/home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/"
files=`ls ${path}clean | grep txt | grep Case`

for file in $files
do echo $file
file_out=${path}h2/${file/.txt.gz/}

if [ ! -f $file_out ]
then
/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/ldsc.py \
--h2 ${path}munge/${file/.txt.gz/}.sumstats.gz \
--ref-ld-chr $path_ld \
--w-ld-chr $path_ld \
--out $file_out
fi
done

## collect
path = "/home/yanglab_data/user/zhanghy/gwas/summary/bbj/h2/"
files = list.files(path)
files = gsub('.log', '', files)

out = c()
for (file in files){
  print(file)
  item = scan(paste0(path, file, '.log'), character(), sep='\t')
  item = item[grepl('h2:', item)]
  item = as.numeric(strsplit(item, ' |[()]')[[1]])
  item = item[!is.na(item)]
  item[3] = 2*pnorm(abs(item[1]/item[2]), lower.tail=FALSE)
  out = c(out, file,  item)
}

collect = data.frame(matrix(out, ncol=4, byrow=T))
colnames(collect) = c('data', 'h2', 'se', 'p')
write.table(collect, paste0(path, '../para/h2.txt'), quote=F, row.names=F, sep='\t')

#=====================================================================================
# nsnp with p<5e-8
#=====================================================================================
## nsnp with p<5e-8 
cd "/home/yanglab_data/user/zhanghy/gwas/summary/bbj/clean/"
datas=`ls | grep txt`
datas=($datas)
len="${#datas[*]}"

rm -rf ../para/nsnp_5e-8.txt 
touch ../para/nsnp_5e-8.txt 

for ((i=0; i<$len; i++)) 
  do
data=${datas[$i]}
echo $data

if [ `cat ../para/nsnp_5e-8.txt |grep $data | wc -l` == 0 ]
then
nsnp=`awk -F '\t' 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} $(f["P"])<5e-8' $data | wc -l`
echo $data $nsnp >> ../para/nsnp_5e-8.txt
if [ $nsnp -ge 10 -a ! -f ../clump/for_clump/$data ]
then
cat <(head -1 $data) <(awk -F '\t' 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} $(f["P"])<5e-8' $data) \
> ../clump/for_clump/$data
fi
fi
done

#=====================================================================================
# nsnp with p<5e-8 after clump
#=====================================================================================
path="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/"
file_bfile="/home/yanglab_data/user/zhanghy/gwas/bfile/1000g/no_mhc/eas/eas_nomhc_maf_0.01"

datas=`awk '$2>=10 {print $1}' $path/para/nsnp_5e-8.txt`
datas=($datas)
len="${#datas[*]}"

for ((i=0; i<$len; i++))
do
data=${datas[$i]}

echo $data
name=${data//.txt.gz/}
file_out=$path/clump/clump_out/$name

if [ ! -f $file_out.clumped ]
then
plink --allow-no-sex \
--clump $path/clump/for_clump/$data \
--bfile $file_bfile \
--clump-kb 1000 --clump-p1 1.0 --clump-p2 1.0 --clump-r2 0.05 \
--out $file_out
fi
done

# collect.r
set = 'bbj'

res = data.frame()
path_in = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/', set, '/clump/clump_out/')
file_out = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/', set, '/para/nsnp_5e-8_clumped.txt')
file_info = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/', set, '/para/nsnp_5e-8.txt')
info = read.table(file_info)

datas = info[info[,2]>=10, 1]
datas_proxy = paste0(c('QTL_43', 'QTL_89', 'QTL_117', 'Case_control_36'), '.txt')
datas = datas[!datas%in%datas_proxy]

out = c()
for (data in datas){
  file_in = paste0(path_in, gsub('.txt', '.clumped', data))
  if (! file.exists(file_in)){
    next()
  }
  df = read.table(file_in, header=T, stringsAsFactors = F)
  out = c(out, data, nrow(df))
}

res = matrix(out, ncol=2, byrow=T)
write.table(res, file_out, row.names = F, col.names = F, quote = F, sep = '\t')
