#=====================================================================================
# ldsc
#=====================================================================================
##ldsc_pair
path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/clean/'
clump = read.table(paste0(path, '../para/nsnp_5e-8_clumped.txt'))
keep = c(gsub('.txt', '', clump[clump[,2]>=10, 1]), 'spracklen_nature2020_t2d', 'spracklen_nature2020_t2d_bmisubset')

datas = list.files(path)
datas = datas[grepl('txt', datas)]
datas = c('spracklen_nature2020_t2d', 'spracklen_nature2020_t2d_bmisubset', gsub('.txt', '', datas))
pair = data.frame(t(combn(datas, 2)))
pair = pair[pair[,1]%in%keep|pair[,2]%in%keep,]

# sort each pair
for (i in 1:nrow(pair)){
  pair[i,] = unlist(sort(pair[i,]))
}

write.table(pair, paste0(path, '../para/ldsc_pair.txt'), row.names = F, quote = F, sep = '\t', col.names = F)

## ldsc
path_bbj='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/'
path_t2d='/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/'
path_ld='/home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/'

data1s=(`cut -f 1  <(cat $path_bbj/para/ldsc_pair.txt| grep -v QTL) `)
data2s=(`cut -f 2  <(cat $path_bbj/para/ldsc_pair.txt| grep -v QTL) `)

len="${#data1s[*]}"

for ((i=1; i<$len; i++))
do
echo $i
data1=${data1s[$i]}
data2=${data2s[$i]}

if [[ $data1 = *"sprack"* ]]
then
file1=$path_t2d/munge/$data1 
else 
file1=$path_bbj/munge/$data1
fi

if [[ $data2 = *"sprack"* ]]
then
file2=$path_t2d/munge/$data2 
else 
file2=$path_bbj/munge/$data2
fi

file_out1=$path_bbj/result/ldsc/${data1}\&${data2}
file_out2=$path_bbj/result/ldsc/${data1}\&${data2}_constrain # for grad use

if [ ! -f $file_out1.log ] 
then
/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/ldsc.py \
--rg  $file1.sumstats.gz,$file2.sumstats.gz \
--ref-ld-chr $path_ld \
--w-ld-chr $path_ld \
--out $file_out1
fi

# if [ ! -f $file_out2.log ] 
# then
/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/ldsc.py \
--rg  $file1.sumstats.gz,$file2.sumstats.gz \
--ref-ld-chr $path_ld \
--w-ld-chr $path_ld \
--intercept-h2 1,1 \
--out $file_out2
# fi

done
#=====================================================================================
# deal with error file 
#=====================================================================================
## check_ldsc_error.sh
path_bbj='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/'

rm -rf ${path_bbj}/result/collect/ldsc_error.log

files=`ls ${path_bbj}/result/ldsc/`
for file in $files
do 
if cat ${path_bbj}/result/ldsc/$file | egrep "Error|error"
then echo $file `cat ${path_bbj}/result/ldsc/$file | egrep "Error|error"` >> ${path_bbj}/result/collect/ldsc_error.log 
fi
done

cat ${path_bbj}/result/collect/ldsc_error.log | awk -F "log" '{print $2}' | sort | uniq -c

## remove file
files=`ls ${path_bbj}result/ldsc |  grep log`
for file in $files
do
#echo $file
if cat $file | grep "at least two phenotypes"
then
rm -rf ${path_bbj}result/ldsc/$file
fi
done

#=====================================================================================
# remove pair with reversed order
#=====================================================================================
pair = read.table('/home/yanglab_data/user/zhanghy/gwas/summary/bbj/para/ldsc_pair.txt')

for (i in 1:nrow(pair)){
  data1 = pair[i, 1]
  data2 = pair[i, 2]
  
  out1 = paste0(data1, '&', data2, '.log')
  out2 = paste0(data1, '&', data2, '_lia.log')
  out3 = paste0(data1, '&', data2, '_constrain_lia.log')
  
  out1_re = paste0(data2, '&', data1, '.log')
  out2_re = paste0(data2, '&', data1, '_lia.log')
  out3_re = paste0(data2, '&', data1, '_constrain_lia.log')
  
  if (file.exists(out1)){
    file.remove(out1_re)
  }
  if (file.exists(out2)){
    file.remove(out2_re)
  }
  if (file.exists(out3)){
    file.remove(out3_re)
  }
}

#=====================================================================================
# collect_ldsc
#=====================================================================================
library(stringr)

remove_brackets_split = function(x){
  if (length(line) >0) {
    x = str_squish(x) # multiple spaces to single
    x = gsub(')', '', gsub('(', '', x, fixed=T), fixed=T)
    x = strsplit(x, ' ')[[1]]
  }
  return(x)
}

path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/result/'

temp = list.files(paste0(path, 'ldsc/'))
files = temp[grepl('log', temp)]
files = files[grepl('14', temp)]
files = files[!grepl('QTL', files)]
files = files[!grepl('spracklen', files)]
files = files[!grepl('lia', files)]

out = c()
file_error = c()
file_warn = c()

col_name = c('data1', 'data2', 'constrain', 'liability', 'nsnp', 'h2_1', 'h2_1_se', 'int_1', 'int_1_se', 'h2_2', 'h2_2_se', 'int_2', 'int_2_se',
             'rg', 'rg_se', 'rg_p', 'int_bi', 'int_bi_se')


for (file in files){
  warn = F; print(file)
  
  if (sum(grepl("WARN", readLines(paste0(path, 'ldsc/', file))))){
    file_warn = c(file_warn, file)
    next()
  }
  data1 = strsplit(file, '&')[[1]][1]
  data2 = gsub('_lia', '', gsub('_constrain', '', gsub('.log', '', strsplit(file, '&')[[1]][2])))
  
  con = file(paste0(path, 'ldsc/', file))
  if (length(readLines(con)) < 62) {
    file_error = c(file_error, file)
    close(con)
    next()
  }
  close(con)
  
  con = file(paste0(path, 'ldsc/', file),open="r")
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
  
  if (warn != T){
    out = c(out, data1, data2, constrain, liability, nsnp, h2_1, h2_1_se, int_1, int_1_se, h2_2, h2_2_se, int_2, int_2_se, 
            rg, rg_se, rg_p, int_bi, int_bi_se)
    rm(list = col_name)
  }
}

res = data.frame(matrix(out, ncol = 18, byrow=T))
colnames(res) = col_name
res = res[res$data1!= '' & res$data2 != '',]

cat('n error file: ', length(file_error), '\n')
cat('n warn file: ', length(file_warn), '\n')

res[,c('h2_1', 'h2_2', 'h2_1_se', 'h2_2_se')] = sapply(res[,c('h2_1', 'h2_2', 'h2_1_se', 'h2_2_se')], as.numeric)
res = res%>%mutate(h2_1_p=2*pnorm(abs(res$h2_1/res$h2_1_se), lower.tail=F), h2_2_p=2*pnorm(abs(res$h2_2/res$h2_2_se), lower.tail=F))

write.csv(res, paste0(path, 'collect/ldsc_t2d.csv'), row.names = F)

#=====================================================================================
# liability ldsc
#=====================================================================================
# t2d lia
data1=Case_control_14_on_QTL_51_ref_ukbchn
data2=Case_control_36

~/soft1/ldsc/ldsc.py \
--rg  $path\munge_out/$data1.sumstats.gz,$path\munge_out/$data2.sumstats.gz \
--ref-ld-chr /home/yanglab_data/pub/gwas_zhanghy/plink_file/1000g/eas_ldscores/ \
--w-ld-chr /home/yanglab_data/pub/gwas_zhanghy/plink_file/1000g/eas_ldscores/ \
--intercept-h2 1,1 \
--out $path\ldsc_out/$data1\&$data2\_constrain_lia \
--samp-prev 0.191,0.1159 \
--pop-prev 0.075,0.0009

~/soft1/ldsc/ldsc.py \
--rg  $path\munge_out/$data1.sumstats.gz,$path\munge_out/$data2.sumstats.gz \
--ref-ld-chr /home/yanglab_data/pub/gwas_zhanghy/plink_file/1000g/eas_ldscores/ \
--w-ld-chr /home/yanglab_data/pub/gwas_zhanghy/plink_file/1000g/eas_ldscores/ \
--out $path\ldsc_out/$data1\&$data2\_lia \
--samp-prev 0.191,0.1159 \
--pop-prev 0.075,0.0009

# t2d lia
path='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/'
data1=Case_control_14
data2=Case_control_36

~/soft1/ldsc/ldsc.py \
--rg  $path\munge/$data1.sumstats.gz,$path\munge/$data2.sumstats.gz \
--ref-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/ \
--w-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/ \
--intercept-h2 1,1 \
--out $path\result/ldsc/$data1\&$data2\_constrain_lia \
--samp-prev 0.191,0.1159 \
--pop-prev 0.075,0.0009

~/soft1/ldsc/ldsc.py \
--rg  $path\munge/$data1.sumstats.gz,$path\munge/$data2.sumstats.gz \
--ref-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/ \
--w-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/ \
--out $path\result/ldsc/$data1\&$data2\_lia \
--samp-prev 0.191,0.1159 \
--pop-prev 0.075,0.0009

# mtcojoed data | without maf filtering
data1s=`ls ~/gwas/data/munge_out/ | grep -E 'Case_control_14*._on_.*ukb.*gz' | sed 's/.sumstats.gz//g'`
data2s=`ls ~/gwas/data/munge_out/ | grep -E 'Case_control_36*._on_.*ukb.*gz' | sed 's/.sumstats.gz//g'`

data1s=( $data1s )
data2s=( $data2s )

# check if equal length and matched, then run
for i in {0..40} 
do
echo $i
data1=${data1s[$i]}
data2=${data2s[$i]}
#if [ `ls $path\ldsc_out | grep $data1\&$data2.log | wc -l` == 0 -o `cat $path\ldsc_out/$data1\&$data2.log | wc -l` == 0 ]
#then
~/soft1/ldsc/ldsc.py \
--rg  $path\munge_out/$data1.sumstats.gz,$path\munge_out/$data2.sumstats.gz \
--ref-ld-chr /home/yanglab_data/pub/gwas_zhanghy/plink_file/1000g/eas_ldscores/ \
--w-ld-chr /home/yanglab_data/pub/gwas_zhanghy/plink_file/1000g/eas_ldscores/ \
--intercept-h2 1,1 \
--out $path\ldsc_out/$data1\&$data2\_constrain_lia \
--samp-prev 0.191,0.1159 \
--pop-prev 0.075,0.0009

~/soft1/ldsc/ldsc.py \
--rg  $path\munge_out/$data1.sumstats.gz,$path\munge_out/$data2.sumstats.gz \
--ref-ld-chr /home/yanglab_data/pub/gwas_zhanghy/plink_file/1000g/eas_ldscores/ \
--w-ld-chr /home/yanglab_data/pub/gwas_zhanghy/plink_file/1000g/eas_ldscores/ \
--out $path\ldsc_out/$data1\&$data2\_lia \
--samp-prev 0.191,0.1159 \
--pop-prev 0.075,0.0009
#fi
done

prev_pop_cataract=0.1879
prev_sam_cataract=0.01396961
prev_pop_t2d=0.1
prev_sam_t2d=0.09538977

# eur
~/soft1/ldsc/ldsc.py \
--rg  /home/yanglab_data/user/zhanghy/gwas/summary/other/t2d-cataract_eur/munge/t2d_eur.sumstats.gz,/home/yanglab_data/user/zhanghy/gwas/summary/other/t2d-cataract_eur/munge/cataract_eur.sumstats.gz \
--ref-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr/ \
--w-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr/ \
--out ~/gwas/data/ldsc/eur/eur_t2d\&eur_cataract_lia \
--samp-prev $prev_sam_t2d,$prev_sam_cataract \
--pop-prev $prev_pop_t2d,$prev_pop_cataract

~/soft1/ldsc/ldsc.py \
--rg  /home/yanglab_data/user/zhanghy/gwas/summary/other/t2d-cataract_eur/munge/t2d_eur.sumstats.gz,/home/yanglab_data/user/zhanghy/gwas/summary/other/t2d-cataract_eur/munge/cataract_eur.sumstats.gz \
--ref-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr/ \
--w-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr/ \
--intercept-h2 1,1 \
--out ~/gwas/data/ldsc/eur/eur_t2d\&eur_cataract_constrain_lia \
--samp-prev $prev_sam_t2d,$prev_sam_cataract \
--pop-prev $prev_pop_t2d,$prev_pop_cataract


# eas and eur counterpart
~/soft1/ldsc/ldsc.py \
--rg  /home/yanglab_data/user/zhanghy/gwas/summary/other/t2d-cataract_eur/munge/cataract_eur.sumstats.gz,/home/yanglab_data/user/zhanghy/gwas/summary/jap/munge/Case_control_36.sumstats.gz \
--ref-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/ \
--w-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/ \
--intercept-h2 1,1 \
--out ~/gwas/data/ldsc/eur/cataract_dou_easref_constrain_lia \
--samp-prev $prev_sam_cataract,0.1159 \
--pop-prev $prev_pop_cataract,0.0009

~/soft1/ldsc/ldsc.py \
--rg  /home/yanglab_data/user/zhanghy/gwas/summary/other/t2d-cataract_eur/munge/cataract_eur.sumstats.gz,/home/yanglab_data/user/zhanghy/gwas/summary/jap/munge/Case_control_36.sumstats.gz \
--ref-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr/ \
--w-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr/ \
--intercept-h2 1,1 \
--out ~/gwas/data/ldsc/eur/cataract_dou_eurref_constrain_lia \
--samp-prev $prev_sam_cataract,0.1159 \
--pop-prev $prev_pop_cataract,0.0009


~/soft1/ldsc/ldsc.py \
--rg  /home/yanglab_data/user/zhanghy/gwas/summary/other/t2d-cataract_eur/munge/t2d_eur.sumstats.gz,/home/yanglab_data/user/zhanghy/gwas/summary/jap/munge/Case_control_14.sumstats.gz \
--ref-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/ \
--w-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/ \
--intercept-h2 1,1 \
--out ~/gwas/data/ldsc/eur/t2d_dou_easref_constrain_lia \
--samp-prev $prev_sam_cataract,0.191 \
--pop-prev $prev_pop_cataract,0.075


~/soft1/ldsc/ldsc.py \
--rg  /home/yanglab_data/user/zhanghy/gwas/summary/other/t2d-cataract_eur/munge/t2d_eur.sumstats.gz,/home/yanglab_data/user/zhanghy/gwas/summary/jap/munge/Case_control_14.sumstats.gz \
--ref-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr/ \
--w-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr/ \
--intercept-h2 1,1 \
--out ~/gwas/data/ldsc/eur/t2d_dou_eurref_constrain_lia \
--samp-prev $prev_sam_t2d,0.191 \
--pop-prev $prev_pop_t2d,0.075

#=====================================================================================
# blood marker
#=====================================================================================
data1=spracklen_nature2020_t2d

data2s=(QTL_7 QTL_9 QTL_11 QTL_13 QTL_15 QTL_17 QTL_21 QTL_23 QTL_25 QTL_27 QTL_29 QTL_31  QTL_39  QTL_43  QTL_45  QTL_47  QTL_49  QTL_53  QTL_57  QTL_59  QTL_61  QTL_71  QTL_75  QTL_77  QTL_79  QTL_81  QTL_83  QTL_85  QTL_87  QTL_89  QTL_93  QTL_97  QTL_99  QTL_105 QTL_107 QTL_109 QTL_111 QTL_113 QTL_115 QTL_117 QTL_119 QTL_121)

for data2 in ${data2s[@]}
do
echo $data2

python ~/soft1/ldsc/ldsc.py \
--rg /home/yanglab_data/user/zhanghy/gwas/summary/other/eas/munge/$data1.sumstats.gz,/home/yanglab_data/user/zhanghy/gwas/summary/bbj/munge/$data2\.sumstats.gz \
--ref-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/ \
--w-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/ \
--pop-prev 0.075,nan \
--samp-prev 0.195,nan \
--intercept-h2 1,1 \
--out /home/yanglab_data/user/zhanghy/gwas/summary/bbj/result/ldsc/$data1\&$data2\_constrain_lia

python ~/soft1/ldsc/ldsc.py \
--rg /home/yanglab_data/user/zhanghy/gwas/summary/other/eas/munge/$data1.sumstats.gz,/home/yanglab_data/user/zhanghy/gwas/summary/bbj/munge/$data2\.sumstats.gz \
--ref-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/ \
--w-ld-chr /home/yanglab_data/user/zhanghy/gwas/plink_file/eas_ldscores/ \
--pop-prev 0.075,nan \
--samp-prev 0.195,nan \
--out /home/yanglab_data/user/zhanghy/gwas/summary/bbj/result/ldsc/$data1\&$data2\_lia

done