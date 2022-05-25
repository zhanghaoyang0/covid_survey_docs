
#=====================================================================================
# download 
#=====================================================================================
### eas
# https://www.linux.com/training-tutorials/wget-and-downloading-entire-remote-directory/
wget -r –level=0 -E –ignore-length -x -k -p -erobots=off -np -N https://ftp.cngb.org/pub/CNSA/data2/CNP0000794/microbiome/


## test
wget -r –level=0 -E –ignore-length -x -k -p -erobots=off -np -N http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/


### eur
info="/home/yanglab_data/user/zhanghy/project/temp/code/gwascatalog_info_20220508.csv"

datas=`cat $info | grep 33462482 | awk -v FS=',' '{print $3}'`
for data in $datas
do 
url="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90011001-GCST90012000/${data}/${data}_buildGRCh37.tsv.gz"
wget -c $url
done

datas=`cat $info | egrep "33208821|33462485" | awk -v FS=',' '{print $3}'`
for data in $datas
do 
url="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90007001-GCST90008000/${data}/${data}_buildGRCh37.tsv.gz"
wget -c $url
done
#=====================================================================================
# unzip 
#=====================================================================================
for i in `ls | grep 7z`
do
7z x $i
done 
#=====================================================================================
# clean 
#=====================================================================================
### prs: select gwas in bfile, map a1, a2, Most PRS software will perform strand-flipping, so not need to filp beta
source('/home/yanglab_data/user/zhanghy/project/temp/code/microbio/source_microbio.r')
pop = 'eas'
map = read.table('/home/yanglab_data/user/zhanghy/db/map/ukb/map_hg38_to_hg19', header = 1, sep = '\t')
bim = read.table(paste0(path_list[[pop]]['bfile'], '.bim'))%>%rename(SNP=V2, A1=V5, A2=V6)
traits = list.files(paste0(path_list[[pop]]['data'], 'raw/CNP0000794/microbiome'))
traits = traits[grepl('7z', traits)]

for (trait in traits){
    print(trait)
    file_out1 = paste0(path_list[[pop]]['data'], 'clean/ukb_ref/', strsplit(trait, '.assoc.')[[1]][1], '.txt')
    file_out2 = paste0(path_list[[pop]]['data'], 'result/prs/clump/for_clump/', strsplit(trait, '.assoc.')[[1]][1], '.txt')
    if (sum(file.exists(file_out1,file_out2))==2 | sum(file.exists(paste0(file_out1, '.gz'),paste0(file_out2, '.gz')))==2){next}
    if (trait %in% c('s_unclassified_Alistipes_sp._HGB5.assoc.raw.7z')){next} # problematic
    df_raw = fread(cmd = paste0('7z e -so ', path_list[[pop]]['data'], 'raw/CNP0000794/microbiome/', trait)) # read 7z
    df = data.frame(df_raw%>%rename(FRQ=freq, BETA=b, SE=se, P=p)%>%mutate(SNP=gsub('chr', '', SNP)))
    cols = str_split_fixed(df$SNP, ":", 2)
    df$CHR_38 = cols[,1]; df$POS_38 = cols[,2]
    df = df%>%select(CHR_38, POS_38, A1, A2, FRQ, BETA, SE, P, N)%>%merge(map, by = c('CHR_38', 'POS_38'))
    df = df%>%filter(!((A1=='T'&A2=='A')|(A1=='A'&A2=='T')|(A1=='C'&A2=='G')|(A1=='G'&A2=='C'))) # new qc
    # remove duplicated snp
    df = df%>%select(CHR, POS, SNP, FRQ, A1, A2, BETA, SE, P, N)%>%na.omit()%>%filter(FRQ>0.01&FRQ<0.99)
    df = df%>%mutate(id=1:nrow(df))%>%arrange(P)%>%filter(!duplicated(SNP))%>%arrange(id)%>%select(-id)
    write.table(df, file_out1, row.names=F, quote=F, sep='\t')
    # df for prs clump
    df2 = prs_filter_gwas(df, bim)
    write.table(df2, file_out2, row.names=F, quote=F, sep='\t')
}

#=====================================================================================
# munge 
#=====================================================================================
path="/home/yanglab_data/user/zhanghy/gwas/summary/microbiome/eas/clean/ukb_ref"
datas=`ls $path`
datas=($datas)
start=0
end="${#datas[*]}"
# end=500
for ((i=$start; i<$end; i++))
do
echo $data
data=${datas[$i]}
file_out=$path/../../munge/${data/.txt.gz/}

if [ ! -f $file_out.sumstats.gz ] 
then
/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/munge_sumstats.py \
--sumstats $path/$data \
--out ${file_out} \
--merge-alleles /home/yanglab_data/user/zhanghy/soft_slurm/ldsc/w_hm3.snplist
fi
done


### check
start=0; end="${#datas[*]}"
for ((i=$start; i<$end; i++))
do
data=${datas[$i]}
file_out=$path/../../munge/${data/.txt.gz/}
t=`zcat $file_out.sumstats.gz  | tail -1  | grep rs | wc -l`
if [ $t != 1 ]
then echo $data
fi
done

#=====================================================================================
# h2 
#=====================================================================================
## ldsc: cal h2
pop="eas"
path="/home/yanglab_data/user/zhanghy/gwas/summary/microbiome/${pop}/"
if [ $pop == "eas" ]
then path_ld="/home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/"
else path_ld="/home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr/"
fi

# run
files=`ls ${path}clean/ukb_ref | grep txt`

for file in $files
do echo $file
file_out=${path}h2/${file/.txt.gz/}
if [ ! -f ${file_out}.log ]
then
/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/ldsc.py \
--h2 ${path}munge/${file/.txt.gz/}.sumstats.gz \
--ref-ld-chr $path_ld --w-ld-chr $path_ld \
--out $file_out
fi
done

### check
for file in $files
do echo $file
file_out=${path}h2/${file/.txt.gz/}
if cat ${file_out}.log  | egrep "Error|error"
then 
rm -rf ${file_out}.log 
fi
done


## collect
pop = 'eas'
path = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/microbiome/', pop, '/h2/')
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
