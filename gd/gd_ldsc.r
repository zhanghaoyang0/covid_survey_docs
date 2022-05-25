# gd
#=====================================================================================
## liability ldsc
#=====================================================================================
path_in_fin="/home/yanglab_data/user/zhanghy/gwas/summary/finngen/"
path_in_giant="/home/yanglab_data/user/zhanghy/gwas/summary/giant/"
patn_in_ukb="/home/yanglab_data/user/zhanghy/gwas/summary/ukb/eur/"
path_ld="/home/yanglab_data/user/zhanghy/gwas/plink_file/eur_w_ld_chr/f/"
path_out="/home/yanglab_data/user/zhanghy/gwas/summary/giant/result/ldsc/"

pop_prev_t1d=0.003
pop_prev_t2d=0.031
pop_prev_gd=0.041

sam_prev_t1d=0.006
sam_prev_t2d=0.042
sam_prev_gd=0.046

data1s=(finr5_263 E11_f E10_f)
data2s=(giant_bmi_f_2015 giant_hipadjbmi_f_2015 giant_hip_f_2015 giant_wcadjbmi_f_2015 giant_wc_f_2015 giant_whr_f_2015 giant_whradjbmi_f_2015 giant_height_f_2013 giant_weight_f_2013 sitting_height_f)

for data1 in ${data1s[@]}
do 
for data2 in ${data2s[@]}
do

if [[ $data1 = *"fin"* ]]
then file1=$path_in_fin/munge/$data1
elif [[ $data1 = *"giant"* ]]
then file1=$path_in_giant/munge/$data1
else file1=$patn_in_ukb/munge/$data1
fi

if [[ $data2 = *"fin"* ]]
then file2=$path_in_fin/munge/$data2
elif [[ $data2 = *"giant"* ]]
then file2=$path_in_giant/munge/$data2
else file2=$patn_in_ukb/munge/$data2
fi

if [ $data1 = 'finr5_263' ]
then pop_prev=$pop_prev_gd; sam_prev=$sam_prev_gd
elif [ $data1 = 'E11_f' ]
then pop_prev=$pop_prev_t2d; sam_prev=$sam_prev_t2d
else pop_prev=$pop_prev_t1d; sam_prev=$sam_prev_t1d
fi

file_out=$path_out/$data1\&$data2\_lia
echo $file_out

# if [ ! -f $file_out.log ]
# then
/home/yanglab_data/user/zhanghy/soft_slurm/ldsc/ldsc.py \
--rg  $file1.sumstats.gz,$file2.sumstats.gz \
--ref-ld-chr $path_ld \
--w-ld-chr $path_ld \
--pop-prev $pop_prev,nan \
--samp-prev $sam_prev,nan \
--out $file_out
# fi

done
done

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

path_input = '/home/yanglab_data/user/zhanghy/gwas/summary/giant/'
file_out = '/home/yanglab_data/user/zhanghy/gwas/summary/giant/result/collect/ldsc.csv'

temp = list.files(path_input)
files = temp[grepl('log', temp)]

out = c()
file_error = c()
file_warn = c()

col_name = c('data1', 'data2', 'constrain', 'liability', 'nsnp', 'h2_1', 'h2_1_se', 'int_1', 'int_1_se', 'h2_2', 'h2_2_se', 'int_2', 'int_2_se',
             'rg', 'rg_se', 'rg_p', 'int_bi', 'int_bi_se')

write(paste(col_name, collapse = ","), file_out)

data1s = c('finr5_263', 'E10_f', 'E11_f')
data2s = c('giant_bmi_f_2015', 'giant_hipadjbmi_f_2015', 'giant_hip_f_2015', 'giant_wcadjbmi_f_2015', 
           'giant_wc_f_2015', 'giant_whr_f_2015', 'giant_whradjbmi_f_2015', 'giant_height_f_2013', 'giant_weight_f_2013',
           'sitting_height_f')

for (data1 in data1s){
  for (data2 in data2s){
    file = paste0(path_input, 'result/ldsc/', data1, '&', data2, '_lia.log')
    warn = F
    
    if (sum(grepl("WARN", readLines(file)))){
      file_warn = c(file_warn, file)
      next()
    }
    
    con = file(file)
    if (length(readLines(con)) < 62) {
      file_error = c(file_error, file)
      close(con)
      next()
    }
    close(con)
    
    con = file(file,open="r")
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
    # if (warn != T){next}
    append = c(data1, data2, constrain, liability, nsnp, h2_1, h2_1_se, int_1, int_1_se, h2_2, h2_2_se, int_2, int_2_se, 
               rg, rg_se, rg_p, int_bi, int_bi_se)
    write(paste(append, collapse = ","), file_out, append=TRUE)
  }
}


cat('n error file: ', length(file_error), '\n')
cat('n warn file: ', length(file_warn), '\n')

# cal h2 pval
df = read.csv(file_out)
df$h1_p = 2*pnorm(abs(df$h2_1/df$h2_1_se), lower.tail=FALSE)
df$h2_p = 2*pnorm(abs(df$h2_2/df$h2_2_se), lower.tail=FALSE)
write.csv(df, file_out, row.names=F)
