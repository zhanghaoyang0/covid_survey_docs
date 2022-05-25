#=====================================================================================
# download
#=====================================================================================
# sex stratified mr
# download like

wget -c http://jenger.riken.jp/104/
  tar -zxvf index.html
rm -rf index.html

#=====================================================================================
# clean
# ***.x: statistics for male;  ***.y: statistics for female
#=====================================================================================
library(dplyr)
library(gtools)

path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/'
file_bim = '/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/chn.bim'

info = read.csv(paste0(path, '/para/info.csv'))
bim = read.table(file_bim)[, c(1, 4, 2)]
bim = bim%>%rename(SNP=V2, CHR=V1, POS=V4)%>%filter(!(CHR==6&POS<=33448354&POS>=28477797))

traits = list.files(paste0(path, 'raw/'))
traits = traits[grepl('Case', traits)]

for (trait in traits){
  print(trait)
  if (sum(file.exists(paste0(path, 'clean/ukb_ref/', trait, '_', c('m', 'f'), '.txt.gz')))==2){next}
  n_x = sum(info[info$id==strsplit(trait, '_')[[1]][3], c('n_ctrl_m', 'n_case_m')]); n_y = sum(info[info$id==strsplit(trait, '_')[[1]][3], c('n_ctrl_f', 'n_case_f')])
  n_x_case = info[info$id==strsplit(trait, '_')[[1]][3], 'n_case_m']; n_y_case = info[info$id==strsplit(trait, '_')[[1]][3], 'n_case_f']
  n_x_ctrl = info[info$id==strsplit(trait, '_')[[1]][3], 'n_ctrl_m']; n_y_ctrl = info[info$id==strsplit(trait, '_')[[1]][3], 'n_ctrl_f']
  
  files = list.files(paste0('/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/raw/', trait))
  files = mixedsort(files[!grepl('READ', files)])
  print(length(files))
  df = data.frame()
  for (file in files){
    print(file)
    sub = read.table(gzfile(paste0(path, '/raw/', trait, '/', file)), header =1)
    sub = sub%>%mutate(CHR = ifelse(CHR=="X", 23, CHR))%>%merge(bim, by=c('CHR', 'POS'))
    df = dplyr::bind_rows(df, sub)
    print(nrow(df))
  }
  df = df%>%mutate(N.x=n_x, N.y=n_y, FRQ.x=(AF.Cases.x*n_x_case+AF.Controls.x*n_x_ctrl)/n_x, FRQ.y=(AF.Cases.y*n_y_case+AF.Controls.y*n_y_ctrl)/n_y)
  for (sex in c('x', 'y')){
    label=ifelse(sex=='x', 'm', 'f')
    file_out = paste0(path, '/clean/ukb_ref/', trait, '_', label, '.txt')
    sub = df[,c('SNP', 'CHR', 'POS', 'Allele1', 'Allele2', paste0(c('BETA.', 'SE.', 'p.value.', 'FRQ.', 'N.'), sex))]
    colnames(sub) = c('SNP', 'CHR', 'POS', 'A2', 'A1', 'BETA', 'SE', 'P', 'FRQ', 'N')
    # filter, remove duplicated snp
    sub = sub%>%select(CHR, POS, SNP, FRQ, A1, A2, BETA, SE, P, N)%>%na.omit()%>%filter(FRQ>0.01&FRQ<0.99)
    sub = sub%>%mutate(id=1:nrow(sub))%>%arrange(P)%>%filter(!duplicated(SNP))%>%arrange(id)%>%select(-id)
    write.table(sub, file_out, row.names=F, quote=F, sep='\t')
  }
}

#=====================================================================================
# clean single sex (not stratify gwas)
#=====================================================================================
library(dplyr)
library(stringr)
path = "/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/raw/single_sex/"
path_bfile = '/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/'
bim = read.table(paste0(path_bfile, 'chn.bim'))[, c(1, 4, 2)]
bim = bim%>%rename(SNP=V2, CHR=V1, POS=V4)%>%filter(!(CHR==6&POS<=33448354&POS>=28477797))
frq = read.table(paste0(path_bfile, 'chn.frq'), header=1)[, c('SNP', 'A1', 'A2', 'MAF')]%>%rename(FRQ=MAF)

# sample size
n = c('Case_control_16_m', 67773+21905, 'Case_control_18_f', 16077+59721, 
      'Case_control_21_m', 36787+24909,'Case_control_23_f', 6376+7975, 'Case_control_27_f', 5004+33679,
      'QTL_2_m', 85894, 'QTL_4_f', 72390, 'QTL_126_m', 23039, 'QTL_128_f', 7379, 
      'QTL_131_m', 58784, 'QTL_133_f', 13871, 'Case_control_40_f', 605+89731, 'Case_control_102_f', 5552+89731)
n = setNames(data.frame(matrix(n, ncol=2, byrow=T)), c('trait', 'N'))

traits = list.files(path); traits = traits[!grepl('READ', traits)]

for (trait in traits){
  print(trait)
  file_out = paste0(path, '../../clean/ukb_ref/', trait, '.txt')
  if(file.exists(file_out)|file.exists(paste0(file_out, '.gz'))){ next}

  if (trait=='Case_control_40_f'){
    df1 = read.table(paste0(path, trait, '/CeCa.auto.rsq07.mac10.txt.gz'), header=1)
    df2 = read.table(paste0(path, trait, '/CeCa.chrx.rsq07.mac10.txt.gz'), header=1)%>%mutate(CHR = ifelse(CHR=="X", 23, CHR))%>%rename(AF.Cases=AF.Cases.x)
    df = dplyr::bind_rows(df1, df2)
    df = df%>%rename(A1=Allele1, A2=Allele2, P=p.value, FRQ=AF.Cases)
  } else {df = read.table(paste0(path, trait), header=1)}
  df$SNP = NULL
  if ('A1Frq'%in%colnames(df)){ df = df%>%rename(FRQ=A1Frq)} 
  if (trait%in%c('QTL_2_m', 'QTL_4_f')){df = df%>%rename(FRQ=Frq, A1=ALT, A2=REF)}
  if (trait=='Case_control_102_f'){ # A1=Allele2 | README
    df = df%>%rename(A1=Allele2, A2=Allele1, FRQ=AF_Allele2, P=p.value)
  } 
  if (trait=='Case_control_27_f'){
    df = df%>%rename(POS=BP)%>%mutate(BETA=log(df$OR), Z=qnorm((1-df$P)/2), SE=abs(BETA/Z))
    df = df%>%merge(bim, by=c('CHR', 'POS'))%>%merge(frq, by=c('SNP', 'A1', 'A2'))%>%select(-SNP)
    df[is.infinite(df$SE), 'SE'] = 0.01 # arbitary value, bz se may be inf!
  }

  df = df%>%mutate(N=n[n$trait==trait, 'N'])
  print(table(df$CHR))
  df1 = df%>%merge(bim, by=c('CHR', 'POS'))%>%select(CHR, POS, SNP, FRQ, A1, A2, BETA, SE, P, N)%>%na.omit()%>%filter(FRQ>0.01&FRQ<0.99)
  df2 = df1%>%mutate(id=1:nrow(df1))%>%arrange(P)%>%filter(!duplicated(SNP))%>%arrange(id)%>%select(-id)

  write.table(df2, file_out, row.names=F, sep='\t', quote=F)
}

#=====================================================================================
# munge 
#=====================================================================================
### ref with chrx
df1 = read.table('/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/w_hm3.snplist', header=1)
df2 = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt.bim')
df2 = df2[df2[,1]==23, c(2, 5, 6)]
names(df2) = c('SNP', 'A1', 'A2')
df = rbind(df1, df2)
write.table(df, '/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/w_hm3.snplist_withx', row.names=F, quote=F, sep='\t')

### munge
path='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/'
file_ref="/home/yanglab_data/user/zhanghy/soft_slurm//ldsc/w_hm3.snplist_withx"
datas=`ls $path/clean/ukb_ref/ | grep txt | grep -v mtag`
datas=($datas)

len="${#datas[*]}"

for ((i=0; i<$len; i++)) 
  do
data=${datas[$i]}
echo $data
out=`echo $data | sed -e "s/.txt.gz//g"`
file_out=$path/munge/$out

if [ ! -f $file_out.sumstats.gz ] 
then
/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/munge_sumstats.py \
--sumstats $path/clean/ukb_ref/$data \
--out $path/munge/$out \
--merge-alleles $file_ref
fi
done

#=====================================================================================
# nsnp with p<5e-8 
#=====================================================================================
## nsnp with p<5e-8 
cd /home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/clean/
file_out="../../para/nsnp_5e-8.txt"
path_clump="../../clump"

datas=`ls | grep -v mtag | grep txt`
datas=($datas)
len="${#datas[*]}"

rm -rf $file_out; touch $file_out 

for ((i=0; i<$len; i++)) 
do
data=${datas[$i]}
echo $data

if [ `cat $file_out |grep $data | wc -l` == 0 ]
then
nsnp=`awk -F '\t' 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} $(f["P"])<5e-8' <(zcat $data) | wc -l`
echo $data $nsnp >> $file_out
if [ $nsnp -ge 10 -a ! -f ${path_clump}/for_clump/$data ]
then
cat <(head -1 <(zcat $data)) <(awk -F '\t' 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} $(f["P"])<5e-8' <(zcat $data)) \
> ${path_clump}/for_clump/${data/.gz/}
fi; fi; done

#=====================================================================================
# nsnp with p<5e-8 after clump
#=====================================================================================
path="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/"
datas=`awk '$2>=10 {print $1}' $path/para/nsnp_5e-8.txt`
datas=($datas)
len="${#datas[*]}"

for ((i=0; i<$len; i++))
do
data=${datas[$i]}

if [[ $data = *"_m"* ]]
then sex="m"; else sex="f"; fi

file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/chn_${sex}"
file_out=$path/clump/clump_out/${data/.txt.gz/}

if [ ! -f $file_out.clumped ]
then
plink --allow-no-sex \
--clump $path/clump/for_clump/$data \
--bfile $file_bfile \
--clump-kb 1000 --clump-p1 1.0 --clump-p2 1.0 --clump-r2 0.05 \
--out $file_out
fi; done

# collect.r
set = 'bbj/stratify'

res = data.frame()
path_in = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/', set, '/clump/clump_out/')
file_out = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/', set, '/para/nsnp_5e-8_clumped.txt')
file_info = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/', set, '/para/nsnp_5e-8.txt')
info = read.table(file_info)

datas = info[info[,2]>=10, 1]

out = c()
for (data in datas){
  file_in = paste0(path_in, gsub('.txt.gz', '.clumped', data))
  if (! file.exists(file_in)){
    next()
  }
  df = read.table(file_in, header=T, stringsAsFactors = F)
  out = c(out, data, nrow(df))
}

res = matrix(out, ncol=2, byrow=T)
write.table(res, file_out, row.names = F, col.names = F, quote = F, sep = '\t')