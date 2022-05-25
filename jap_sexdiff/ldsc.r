#=====================================================================================
# ldsc
#=====================================================================================
##ldsc_pair
path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/'

clump = read.table(paste0(path, 'para/nsnp_5e-8_clumped.txt'))
keep = c(clump[clump[,2]>=10, 1], 'spracklen_nature2020_t2d_m.txt', 'spracklen_nature2020_t2d_f.txt')
keep = gsub('.txt', '', keep)

datas = list.files(paste0(path, '/clean/'))
datas = datas[grepl('txt', datas)]; datas = datas[!grepl('mtag', datas)]
datas = gsub('.txt', '', datas)

datas_m = c(datas[grepl('m', datas)], 'spracklen_nature2020_t2d_m'); datas_m = gsub('.gz', '', datas_m)
datas_f = c(datas[grepl('f', datas)], 'spracklen_nature2020_t2d_f'); datas_f = gsub('.gz', '', datas_f)

pair = rbind(t(combn(datas_m, 2)), t(combn(datas_f, 2)))
pair = pair[pair[,1]%in%keep | pair[,2]%in%keep,]

write.table(pair, paste0(path, 'para/ldsc_pair.txt'), sep='\t', quote=F, row.names=F, col.names=F)

## ldsc
path_bbj='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/'
path_t2d='/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/'
path_ld='/home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/'

data1s=`cut -f 1  $path_bbj/para/ldsc_pair.txt`
data2s=`cut -f 2  $path_bbj/para/ldsc_pair.txt`
data1s=( $data1s ); data2s=( $data2s )

len="${#data1s[*]}"

for ((i=1; i<$len; i++))
do
echo $i
data1=${data1s[$i]}
data2=${data2s[$i]}

if [[ $data1 = *"t2d"* ]]
then file1=$path_t2d/munge/$data1 
else file1=$path_bbj/munge/$data1
fi

if [[ $data2 = *"t2d"* ]]
then file2=$path_t2d/munge/$data2 
else file2=$path_bbj/munge/$data2
fi

file_out=$path_bbj/result/ldsc/$data1\&$data2

if [ ! -f $file_out.log ] 
then
/home/yanglab_data/user/zhanghy/soft_slurm//ldsc/ldsc.py \
--rg  $file1.sumstats.gz,$file2.sumstats.gz \
--ref-ld-chr $path_ld --w-ld-chr $path_ld \
--out $file_out
fi

done

#=====================================================================================
# deal with error file 
#=====================================================================================
## check_ldsc_error.sh
path='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ldsc/'
file_out='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/collect/ldsc_error.log'

rm -rf $file_out

files=`ls $path`
for file in $files
do 
if cat ${path}${file} | egrep "Error|error"
then echo ${path}${file} `cat ${file} | egrep "Error|error"` >> $file_out
fi
done

cat $file_out | awk -F "log" '{print $2}' | sort | uniq -c

## remove file
datas=`ls path`

for data in $datas
do
#echo $data
if cat $data | grep "at least two phenotypes"
then
rm -rf $path$data
fi
done

#=====================================================================================
# lia and constrain
#=====================================================================================
### bmi and t2d
file_script='/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/ldsc.py'
file1='/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/munge/spracklen_nature2020_t2d_f.sumstats.gz'
file2='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/munge/QTL_4_f.sumstats.gz'
path_ld='/home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/'
file_out='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ldsc/spracklen_nature2020_t2d_f&QTL_4_f'

$file_script \
--rg $file1,$file2 \
--ref-ld-chr $path_ld --w-ld-chr $path_ld \
--pop-prev 0.031,nan --samp-prev 0.03802879,nan \
--intercept-h2 1,1 \
--out ${file_out}\_constrain_lia

$file_script \
--rg $file1,$file2 \
--ref-ld-chr $path_ld --w-ld-chr $path_ld \
--pop-prev 0.031,nan --samp-prev 0.03802879,nan \
--out ${file_out}\_lia

#=====================================================================================
# collect_ldsc
#=====================================================================================
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_ldsc.r')

path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ldsc/'
files = list.files(path)
files = files[grepl('log', files)]
res = collect_ldsc(path, files)

write.csv(res, paste0(path, '../collect/ldsc.csv'), row.names = F)