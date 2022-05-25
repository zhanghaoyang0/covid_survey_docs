# gwascatalog 
#=====================================================================================
# clean  
#=====================================================================================
library(dplyr)
library(stringr)
library(readxl)

path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/'
file_bim = '/bigdat1/user/zhanghy/gwas/bfile/1000g/no_mhc/eur/eur_nomhc.bim'
file_frq = '/bigdat1/user/zhanghy/gwas/bfile/1000g/no_mhc/eur/eur_nomhc.frq'
file_map = '/bigdat1/user/zhanghy/gwas/summary/finngen/para/map_hg38_to_hg19'
file_info = paste0(path, 'info/gwascatalog_info.xlsx')

info = data.frame(read_excel(file_info, sheet=1))
info = info[is.na(info$drop),]
bim = read.table(file_bim)[, c(2, 1, 4)]
frq = read.table(file_frq, header=1)[, c(1:5)]
map = read.table(file_map, header=1)
names(frq)[5] = 'FRQ'

colnames(bim) = c('SNP', 'CHR', 'POS')

for (i in 1:nrow(pair)){
  print(i)
  file = info[i, 'file']
  print(file)
  
  file_out = paste0(path, 'clean/', info[i, 'label'], '.txt')
  if (file.exists(file_out)){next}
  
  # for file with miss cell
  if (file == 'urate_chr1_22_LQ_IQ06_mac10_all_741_rsid.txt.gz'){
    df = read.table(gzfile(paste0(path, 'raw/', file)), header=1, fill=T)
    df1 = df[is.na(df[,ncol(df)]),]
    df2 = df[!is.na(df[,ncol(df)]),]
    df1$n_total_sum = NULL
    df2$RSID = NULL
    names(df1) = names(df2)
    df = rbind(df1, df2)
  } else {
    df = read.table(gzfile(paste0(path, 'raw/', file)), header=1)
  }
  
  # rename col
  for (col in names(info)[8:17]){
    names(df)[names(df)==info[i, col]] = col
  }
  df[,c('A1', 'A2')] = sapply(df[,c('A1', 'A2')], toupper)
  # for miss beta
  if (!'BETA' %in% names(df)){
    df$BETA = log(df$OR)
  }
  # for missing chr, pos
  if (!'CHR' %in% names(df)){
    col = str_split_fixed(str_split_fixed(df$SNP, "_", 2)[,1], ':', 2)
    df$CHR = col[,1]
    df$POS = col[,2]
  }
  # add 1000g SNP
  df$SNP = NULL
  if (info[i, 'build']=='hg38'){
    df = df%>%rename(CHR_38=CHR, POS_38=POS)
    df = df %>% merge(map, by = c('CHR_38', 'POS_38'))
  } else {
    df = df %>% merge(bim, by = c('CHR', 'POS'))
  }
  # for missing frq
  if (!'FRQ'%in%names(df) | sum(is.na(df$FRQ))/nrow(df)>0.5){
    df$FRQ=NULL
    df = df%>%merge(frq)
  }
  # sample size | \\D+ specifies one or more non-numeric characters
  n = gsub(',', '', info[i,'sample'])
  n = sum(as.numeric(strsplit(n, "\\D+")[[1]]))
  df$N = n
  
  df = df[,c('SNP', 'CHR', 'POS', 'A1', 'A2', 'FRQ', 'BETA', 'SE', 'P', 'N')]
  if (ncol(df)<10){
    stop('ncol of df is less than 10!')
  }
  
  df$P = as.numeric(df$P)
  df = na.omit(df)
  df = df[df$FRQ>0.01&df$FRQ<0.99,]
  
  # remove duplicated snp
  df$order = 1:nrow(df)
  df = df[order(df$P),]
  df = df[!duplicated(df$SNP),]
  df = df[order(df$order),]
  df$order = NULL
  
  write.table(df, file_out, sep='\t', quote=F, row.names=F)
}

#=====================================================================================
# munge 
#=====================================================================================
path='/bigdat1/user/zhanghy/gwas/summary/gwascatalog/'
mkdir ${path}munge

files=`ls $path/clean/ | grep txt`
files=($files)

start=0
end=$[${#files[@]}-1]
  
for i in `seq $start $end` 
do
file=${files[$i]}
echo $file
out=`echo $file | sed -e "s/.txt//g"`
if [ ! -f ${path}munge/$out.sumstats.gz ]
then
~/soft1/ldsc/munge_sumstats.py \
--sumstats $path/clean/$file \
--out $path/munge/$out \
--merge-alleles ~/soft1/ldsc/w_hm3.snplist
fi
done

#=====================================================================================
# h2 
#=====================================================================================
## ldsc: cal h2
path="/bigdat1/user/zhanghy/gwas/summary/gwascatalog/"
mkdir ${path}h2

files=`ls $path/clean`
files=($files)

start=0
end=$[${#files[@]}-1]
  
for i in `seq $start $end` 
do
file=${files[$i]}
echo $file
name=${file/.txt/}

if [ `ls $path\h2 | grep $name.log | wc -l` == 0 ]
then
~/soft1/ldsc/ldsc.py \
--h2 $path\munge/$name.sumstats.gz \
--ref-ld-chr /bigdat1/user/zhanghy/gwas/plink_file/eur_w_ld_chr/ \
--w-ld-chr /bigdat1/user/zhanghy/gwas/plink_file/eur_w_ld_chr/ \
--out $path\h2/$name
fi

done

## collect
path = "/bigdat1/user/zhanghy/gwas/summary/gwascatalog/h2/"
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
  out = c(out, paste0(file, '.txt'), item)
}

collect = data.frame(matrix(out, ncol=4, byrow=T))
colnames(collect) = c('data', 'h2', 'se', 'p')
write.table(collect, paste0(path, '../para/h2.txt'), quote=F, row.names=F, sep='\t')

#=====================================================================================
# nsnp with p<5e-8
#=====================================================================================
cd /bigdat1/user/zhanghy/gwas/summary/gwascatalog/clean
ls | wc -l
files=`ls`
files=($files)

start=0
end=$[${#files[@]}-1]

rm -rf ../para/nsnp_5e-8.txt 
touch ../para/nsnp_5e-8.txt 

for i in `seq $start $end`
do
file=${files[$i]}
echo $file

if [ `cat ../para/nsnp_5e-8.txt |grep $file | wc -l` == 0 ]
then
nsnp=`awk -F '\t' 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} $(f["P"])<5e-8' $file | wc -l`
echo $file $nsnp >> ../para/nsnp_5e-8.txt
h2_p=`awk -v file="$file" '$1 == file {print $3}' ../para/h2.txt`
h2_is_sig=`expr $h2_p \< 0.05`

if [ $nsnp -ge 10 -a h2_is_sig ]
then
cat <(head -1 $file) <(awk -F '\t' 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} $(f["P"])<5e-8' $file) \
> ../clump/for_clump/5e-8/$file
fi
fi

done

#=====================================================================================
# nsnp with p<5e-8 after clump
#=====================================================================================
path_in="/bigdat1/user/zhanghy/gwas/summary/gwascatalog/clump/for_clump/5e-8/"
path_out="$path_in/../../clump_out/5e-8/"
path_bfile="/bigdat1/user/zhanghy/gwas/bfile/1000g/no_mhc/eur/eur_nomhc_maf_0.01"

files=`awk '$2>=10 {print $1}' $path_in/../../../para/nsnp_5e-8.txt`
files=($files)

start=0
end=$[${#files[@]}-1]
  
for i in `seq $start $end`
do
file=${files[$i]}
echo $file
name=${file//.txt/}
if [ ! -f $path_out$name\.clumped ]
then
plink --allow-no-sex \
--clump $path_in$file \
--bfile $path_bfile \
--clump-kb 1000 --clump-p1 1.0 --clump-p2 1.0 --clump-r2 0.05 \
--out $path_out$name
fi
done

# collect.r
set = 'gwascatalog'

res = data.frame()
path_in = paste0('/bigdat1/user/zhanghy/gwas/summary/', set, '/clump/clump_out/5e-8/')
path_out = paste0('/bigdat1/user/zhanghy/gwas/summary/', set, '/para/nsnp_5e-8_clumped.txt')
path_info = paste0('/bigdat1/user/zhanghy/gwas/summary/', set, '/para/nsnp_5e-8.txt')
info = read.table(path_info)

files = info[info[,2]>=10, 1]

out = c()
for (file in files){
  df = read.table(paste0(path_in, gsub('.txt', '.clumped', file)), header=T, stringsAsFactors = F)
  out = c(out, file, nrow(df))
}

res = matrix(out, ncol=2, byrow=T)
write.table(res, path_out, row.names = F, col.names = F, quote = F, sep = '\t')