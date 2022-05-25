# gestd
#=====================================================================================
# clean  
# 2015 is unknown
# 2018 mix ukb, drop
# we use 2015 and 2013
#=====================================================================================
library(dplyr)
library(stringr)

path = '/home/yanglab_data/user/zhanghy/gwas/summary/giant/'
bim = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt.bim')[,c(2,1,4)]
bim = bim%>%rename(SNP=V2, CHR=V1, POS=V4)%>%filter(!(CHR==6&POS<=33448354&POS>=28477797))

info = read.csv(paste0(path, 'info/info.csv'))
info = info[grepl('2015', info$label)|grepl('2013', info$label),]
traits = info[,1]

for (trait in traits[5:9]){
  print(trait)
  file_out = paste0(path, 'clean/', info[info$trait==trait, 'label'], '.txt')
  if (file.exists(file_out)){
    next()
  }
  
  # read
  if (trait%in%traits[grepl('2015', traits)&(!grepl('HIPadjBMI', traits))]){
    df_raw = read.table(paste0(path, 'raw/f/', trait), skip=8)
  } else{
    df_raw = read.table(paste0(path, 'raw/f/', trait), header=1)
  }
  # set colnames
  df = df_raw
  if (ncol(df)==8){
    colnames(df) = c('SNP', 'A1', 'A2', 'FRQ', 'BETA', 'SE', 'P', 'N')
  }  else {
    colnames(df) = c('SNP', 'Chr', 'Pos', 'A1', 'A2', 'FRQ', 'BETA', 'SE', 'P', 'N')
  }
  df = df%>%merge(bim, by=c('SNP'))
  df$A1=toupper(df$A1)
  df$A2=toupper(df$A2)
  
  # clean
  df = df%>%select(SNP, CHR, POS, A1, A2, FRQ, BETA, SE, P, N)%>%mutate(A1=toupper(A1), A2=toupper(A2), P=as.numeric(P))%>%na.omit()%>%filter(FRQ>0.01&FRQ<0.99)

  # remove duplicated snp
  df = df%>%mutate(id=1:nrow(df))%>%arrange(P)%>%filter(!duplicated(SNP))%>%arrange(id)%>%select(-id)

  write.table(df, file_out, quote = F, row.names = F, sep = '\t')
}

## gestd
# GEST_DIABETES   finr5_263       5687    117892  Gestational diabetes (for exclusion)
library(dplyr)

map = read.table('/home/yanglab_data/user/zhanghy/db/map/ukb_brt/map_hg38_to_hg19', header=1, sep='\t')

df_raw = read.table(gzfile('/home/yanglab_data/user/zhanghy/gwas/summary/finngen/raw/finngen_R5_GEST_DIABETES.gz'), header=1, sep='\t', comment.char='', fill=T)
df = df_raw%>%rename(CHR_38=X.chrom, POS_38=pos, A1=alt, A2=ref, BETA=beta, SE=sebeta, FRQ=maf, P=pval)%>%filter(FRQ>0.01&FRQ<0.99)%>%na.omit()%>%
    select(CHR_38, POS_38, A1, A2, FRQ, BETA, SE, P)%>%mutate(N=5687+117892)%>%mutate(CHR_38=ifelse(CHR_38=='X', '23', CHR_38))

df = df%>%merge(map, by = c('CHR_38', 'POS_38'))%>%select(-CHR_38, -POS_38)%>%filter(!(CHR==6&POS<=33448354&POS>=28477797))
df = df%>%mutate(id=1:nrow(df))%>%arrange(P)%>%filter(!duplicated(SNP))%>%arrange(id)%>%select(-id)
df = df[,c('A1', 'A2', 'FRQ', 'BETA', 'SE', 'P', 'CHR', 'POS', 'SNP', 'N')]

write.table(df, '/home/yanglab_data/user/zhanghy/gwas/summary/finngen/clean/finr5_263.txt', row.names=F, quote=F, sep='\t')
#=====================================================================================
# munge | munge first, since some has error:
#=====================================================================================
path="/home/yanglab_data/user/zhanghy/gwas/summary/giant/"
files=`ls $path/clean/ | grep txt`
files=($files)

for i in {0..8}
do
file=${files[$i]}
echo $file
out=`echo $file | sed -e "s/.txt//g"`
if [ ! -f $path/munge/$out.sumstats.gz ]
then
/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/munge_sumstats.py \
--sumstats $path/clean/$file \
--out $path/munge/$out \
--merge-alleles /home/yanglab_data/user/zhanghy/soft_slurm/ldsc/w_hm3.snplist_withx
fi
done

#=====================================================================================
# h2 
#=====================================================================================
## ldsc: cal h2
path="/home/yanglab_data/user/zhanghy/gwas/summary/giant/"
path_ld="/home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr/"

files=`ls $path/munge| grep sum`
for file in $files
do
echo $file
name=${file/.sumstats.gz/}

if [ ! -f $path\h2/$name.log ]
then
/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/ldsc.py \
--h2 $path\munge/$name.sumstats.gz \
--ref-ld-chr $path_ld \
--w-ld-chr $path_ld \
--out $path\h2/$name
fi

done

## collect
path = "/home/yanglab_data/user/zhanghy/gwas/summary/giant/h2/"
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
## nsnp with p<5e-8 in 
cd /home/yanglab_data/user/zhanghy/gwas/summary/giant/clean
files=`ls ../munge/ | grep sum`

rm -rf ../para/nsnp_5e-8.txt
touch ../para/nsnp_5e-8.txt

p_col=9
for file in $files
do
echo $file
file=${file/sumstats.gz/txt}

if [ `cat ../para/nsnp_5e-8.txt |grep $file | wc -l` -eq 0 ]
then
nsnp=`awk -v p_col=$p_col '$p_col<5e-8 {print $0}' $file | wc -l`
echo -e $file'\t'$nsnp >> ../para/nsnp_5e-8.txt
h2_p=`awk -v file="$file" '$1 == file {print $3}' ../para/h2.txt`
h2_is_sig=`expr $h2_p \< 0.05`

if [ $nsnp -ge 10 -a h2_is_sig ]
then
cat <(head -1 $file) <(awk -v p_col=$p_col '$p_col<5e-8 {print $0}' $file) > ../clump/for_clump/$file
fi
fi

done

#=====================================================================================
# nsnp with p<5e-8 after clump
#=====================================================================================
## clump
path_in="/home/yanglab_data/user/zhanghy/gwas/summary/giant/clump/for_clump/"
path_out="$path_in/../clump_out/"
path_bfile="/home/yanglab_data/user/zhanghy/gwas/bfile/1000g/no_mhc/eur/eur_nomhc_maf_0.01"
files=`awk '$2>=10 {print $1}' $path_in/../../para/nsnp_5e-8.txt`

for file in $files
do
echo $file
name=${file//.txt/}

if [ ! -f  $path_out$name.clumped ]
then
plink --allow-no-sex \
--clump $path_in$file \
--bfile $path_bfile \
--clump-kb 1000 --clump-p1 1.0 --clump-p2 1.0 --clump-r2 0.05 \
--out $path_out$name
fi
done

## collect.r
res = data.frame()
path_in = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/giant/clump/clump_out/')
file_out = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/giant/para/nsnp_5e-8_clumped.txt')
file_info = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/giant/para/nsnp_5e-8.txt')
info = read.table(file_info)

files = info[info[,2]>=10, 1]

out = c()
for (file in files){
  df = read.table(paste0(path_in, gsub('.txt', '.clumped', file)), header=T, stringsAsFactors = F)
  out = c(out, file, nrow(df))
}

res = matrix(out, ncol=2, byrow=T)
write.table(res, file_out, row.names = F, col.names = F, quote = F, sep = '\t')
