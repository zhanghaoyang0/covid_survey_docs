#=====================================================================================
# ldsc pair
#=====================================================================================
require(data.table)
library(readxl)

path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/'
file_info_gwascata = paste0(path, 'info/gwascatalog_info.xlsx') 

info_gwascata = data.frame(read_excel(file_info_gwascata, sheet=1))
info_gwascata = info_gwascata[is.na(info_gwascata$drop),]
h2_gwascata = read.table(paste0(path, 'para/h2.txt'), sep = '\t', header = 1)
h2 = rbind(h2_gwascata)
keep1 = gsub('.txt', '', h2[h2$h2>0&h2$p<0.05, 'data'])

nsnp_clump_gwascata = read.table(paste0(path, 'para/nsnp_5e-8_clumped.txt'))
nsnp_clump = rbind(nsnp_clump_gwascata)
keep2 = nsnp_clump[nsnp_clump[,2]>=10, 1]
keep2 = sapply(keep2, function(x){gsub('.txt', '', x)})

pair = expand.grid(info_gwascata[info_gwascata$group=='marker', 'label'], info_gwascata[info_gwascata$group=='disease', 'label'])
pair = pair[pair[,1]%in%keep1 & pair[,2]%in%keep1,] # both h2 p < 0.05
pair = pair[pair[,1]%in%keep2 | pair[,2]%in%keep2,] # contain a trait with >=10 iv
names(pair) = c('data1', 'data2')

write.table(pair, paste0(path, 'para/ldsc_pair.txt'), sep='\t', quote=F, row.names=F, col.names=F)
#=====================================================================================
# ldsc 
#=====================================================================================
## ldsc.sh 
#!/usr/bin/bash
batch="$1"
size="$2"

batch_2=`expr $batch + 1`

start=`expr $size \* $batch + 1`
end=`expr $size \* $batch_2 + 1`

start=1
end=99

echo start: $start
echo end: $end

path_in="/bigdat1/user/zhanghy/gwas/summary/gwascatalog/"
path_ld="/bigdat1/user/zhanghy/gwas/plink_file/eur_w_ld_chr/"
path_out="/bigdat1/user/zhanghy/gwas/summary/gwascatalog/result/ldsc/"

for (( i = $start; i <= $end; i++ )) 
  do
echo $i
data1=`awk -v i="$i" 'FNR==i { print $1}' $path_in/para/ldsc_pair.txt`
data2=`awk -v i="$i" 'FNR==i { print $2}' $path_in/para/ldsc_pair.txt`

if [[ $data1 = *"GCST"* ]]
then
file1=/bigdat1/user/zhanghy/gwas/summary/gwascatalog/munge/$data1
elif [[ $data1 = *"fin"* ]]
then
file1=/bigdat1/user/zhanghy/gwas/summary/gwascatalog/munge/$data1
fi

if [[ $data2 = *"GCST"* ]]
then
file2=/bigdat1/user/zhanghy/gwas/summary/gwascatalog/munge/$data2
elif [[ $data2 = *"fin"* ]]
then
file2=/bigdat1/user/zhanghy/gwas/summary/gwascatalog/munge/$data2
fi

file_out=$path_out/$data1\&$data2

if [ ! -f $file_out.log ]
then
/bigdat1/user/zhanghy/soft_slurm/ldsc/ldsc.py \
--rg  $file1.sumstats.gz,$file2.sumstats.gz \
--ref-ld-chr $path_ld \
--w-ld-chr $path_ld \
--out $file_out
fi
done

## ldsc_sbatch.sh
#!/usr/bin/bash
batch="$1"
size="$2"
bash ldsc.sh $batch $size 

## run
for i in {21..22}  # 1..354
do
node=
  sbatch -N1 -n1 -c4 -w compute$node ldsc_sbatch.sh $i 50
done

#=====================================================================================
# deal with error file 
#=====================================================================================
path="/bigdat1/user/zhanghy/gwas/summary/gwascatalog/"

rm -rf ${path}reuslt/collect/ldsc_error.log

files=`ls ${path}result/ldsc/`
for file in $files
do 
if cat ${path}/result/ldsc/$file | egrep "Error|error"
then echo $file `cat ${path}/result/ldsc/$file | egrep "Error|error"` >> ${path}/result/collect/ldsc_error.log 
fi
done

cat ${path}/result/collect/ldsc_error.log | awk -F "log" '{print $2}' | sort | uniq -c


#=====================================================================================
# collect_ldsc
#=====================================================================================
## collect_ldsc.r
library(stringr)

remove_brackets_split = function(x){
  if (length(line) >0) {
    x = str_squish(x) # multiple spaces to single
    x = gsub(')', '', gsub('(', '', x, fixed=T), fixed=T)
    x = strsplit(x, ' ')[[1]]
  }
  return(x)
}

path_input = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/result/ldsc/'
file_out = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/result/collect/ldsc.csv'

temp = list.files(path_input)
files = temp[grepl('log', temp)]

out = c()
file_error = c()
file_warn = c()

col_name = c('data1', 'data2', 'constrain', 'liability', 'nsnp', 'h2_1', 'h2_1_se', 'int_1', 'int_1_se', 'h2_2', 'h2_2_se', 'int_2', 'int_2_se',
             'rg', 'rg_se', 'rg_p', 'int_bi', 'int_bi_se')

write(paste(col_name, collapse = ","), file_out)

for (file in files){
  warn = F
  
  if (sum(grepl("WARN", readLines(paste0(path_input, file))))){
    file_warn = c(file_warn, file)
    next()
  }
  
  data1 = strsplit(file, '&')[[1]][1]
  data2 = gsub('_lia', '', gsub('_constrain', '', gsub('.log', '', strsplit(file, '&')[[1]][2])))
  
  con = file(paste0(path_input, file))
  if (length(readLines(con)) < 62) {
    file_error = c(file_error, file)
    close(con)
    next()
  }
  close(con)
  
  con = file(paste0(path_input, file),open="r")
  n = 1
  while ( TRUE ) {
    line = readLines(con, n = 1)
    # cat(n,line,"\n") # first print and see, then write code
    line = remove_brackets_split(line)
    if ( n > 64 ) {
      break
    }
    if (grepl('_lia', file)){
      liability = 1; shift = 2} else {liability = 0; shift = 0}
    if (grepl('_constrain', file)){
      constrain = 1
      int_1 = 1
      int_1_se = 0
      int_2 = 1
      int_2_se = 0
      int_bi = 1
      int_bi_se = 0
      if (n == 30 + shift){
        nsnp = line[1]
      }
      if (n == 34 + shift){
        h2_1 = line[5]
        h2_1_se = line[6]
      }
      if (n == 41 + shift){
        h2_2 = line[5]
        h2_2_se = line[6]
      }
      if (n == 60 + shift){
        col = line
      }
      if (n == 61 + shift){
        value = line
        rg = value[3]
        rg_se = value[4]
        rg_p = value[6]
      }
    } else {
      constrain = 0
      if (n == 29 + shift){
        nsnp = line[1]
      }
      if (n == 33 + shift){
        h2_1 = line[5]
        h2_1_se = line[6]
      }
      if (n == 36 + shift){
        int_1 = line[2]
        int_1_se = line[3]
      }
      if (n == 41 + shift){
        h2_2 = line[5]
        h2_2_se = line[6]
      }
      if (n == 44 + shift){
        int_2 = line[2]
        int_2_se = line[3]
      }
      if (n == 61 + shift){
        col = line
      }
      if (n == 62 + shift){
        value = line
        rg = value[3]
        rg_se = value[4]
        rg_p = value[6]
        int_bi = value[11]
        int_bi_se = value[12]
      }
    }
    
    n = n+1
  }
  close(con)
  # if (as.numeric(rg_p)>0.05 | length(as.numeric(rg_p))==0){next()}
  if (warn != T){
    print(file)
    append = c(data1, data2, constrain, liability, nsnp, h2_1, h2_1_se, int_1, int_1_se, h2_2, h2_2_se, int_2, int_2_se, 
               rg, rg_se, rg_p, int_bi, int_bi_se)
    write(paste(append, collapse = ","), file_out, append=TRUE)
    #rm(list = col_name)
  }
}

cat('n error file: ', length(file_error), '\n')
cat('n warn file: ', length(file_warn), '\n')

#=====================================================================================
# check and save miss pair
#=====================================================================================
path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/'
pair = read.table(paste0(path, 'para/ldsc_pair.txt'))
pair_out = paste0(path, 'para/ldsc_pair1.txt')
file.remove(pair_out)

for (i in 1:nrow(pair)){
  if (i%%100==0){
    print(i)
  }
  file1 = pair[i, 1]
  file2 = pair[i, 2]
  file = paste0(path, 'result/ldsc/', file1, '&',file2, '.log')
  if (!file.exists(file)){
    write(paste(c(file1, file2), collapse = " "), pair_out, append=TRUE)
  }
}

## remove file not in pair
path = '/bigdat1/user/zhanghy/gwas/summary/gwascatalog/result/ldsc/'
files = list.files(path)
files_in_pair = paste0(pair[,1], '&', pair[,2], '.log')
file.remove(paste0(path, files[!files%in%files_in_pair]))