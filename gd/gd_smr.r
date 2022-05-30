#=====================================================================================
# convert to cojo format
#=====================================================================================
source('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/source/source_mr.r')
get_gd_info(F)

datas = c('E11_f', 'finr5_263', 'giant_wc_f_2015', 'giant_whr_f_2015', 'giant_bmi_f_2015', 'giant_hip_f_2015',
  'giant_whradjbmi_f_2015', 'giant_wcadjbmi_f_2015.txt', 'giant_hipadjbmi_f_2015.txt.gz')

for (data in c('giant_whradjbmi_f_2015', 'giant_wcadjbmi_f_2015', 'giant_hipadjbmi_f_2015')){
  print(data)
  path_in = path_list[sapply(names(path_list), function(x){grepl(x, data)})][[1]]
  file_out = paste0(path_list[['giant']], 'result/gcta/for_cojo/', data)
  if (file.exists(file_out)){next}
  df = read.table(paste0(path_in, 'clean/ukb_ref/', data, '.txt.gz'), header=1)
  df = df%>%select(SNP, A1, A2, FRQ, BETA, SE, P, N)%>%rename(freq=FRQ, b=BETA, se=SE, p=P)
  write.table(df, file_out, row.names = F, sep = '\t', quote = F)
}


# eas bmi | for grad
data='QTL_4_f'
file_in = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/clean/', data, '.txt.gz')
file_out = paste0(path_list[['giant']], 'result/gcta/for_cojo/', data)
df = read.table(file_in, header=1)
df = df%>%select(SNP, A1, A2, FRQ, BETA, SE, P, N)%>%rename(freq=FRQ, b=BETA, se=SE, p=P)
write.table(df, file_out, row.names = F, sep = '\t', quote = F)




#=====================================================================================
# pre- and post-meta smr
#=====================================================================================
path_in="/home/yanglab_data/user/zhanghy/gwas/summary/giant/result/gcta/for_cojo/"
path_out="/home/yanglab_data/user/zhanghy/gwas/summary/giant/result/gcta/smr/"
path_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/chn_f"

for data in `ls $path_in`
do
/home/yanglab_data/user/zhanghy/soft_slurm/gcta_1.93.2beta/smr_Linux \
--bfile $path_bfile \
--gwas-summary $path_in$data \
--beqtl-summary /home/yanglab_data/user/zhanghy/gwas/smr_file/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense \
--diff-freq-prop 0.5 \
--out $path_out\eqtl/$data
done

# compare
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/gwas/summary/giant/result/gcta/smr/'
map = read.csv('/home/yanglab_data/user/zhanghy/db/source/gencode_summ_hg19.csv')
gene_unique = names(table(map$gene))[table(map$gene)==1]
map = map[map$gene %in% gene_unique,]

res = data.frame()
datas = list.files(paste0(path, 'eqtl')); datas = datas[grepl('smr', datas)]; datas = gsub('.smr', '', datas)

for (qtl in c('eqtl')){
  print(qtl)
  for (data in datas){
    print(data)
    df = read.csv(paste0(path, qtl, '/', data, '.smr'), sep = '\t')
    print(nrow(df))
    if (nrow(df) == 0){next}
    df = cbind(qtl, data, df)
    res = rbind(res, df)
  }
}

res = res%>%na.omit()%>%filter(!grepl('mtag', data))

res1 = res[res$p_HEIDI>=0.05 & res$nsnp_HEIDI>=10 & !is.na(res$p_SMR),]

res2 = data.frame()
for (data_ in datas){
  sub = res1%>%filter(data==data_)
  sub = sub%>%mutate(ntest=nrow(sub))%>%filter(p_SMR<0.05/nrow(sub))
  res2 = rbind(res2, sub)
}

res2$gene = NA
for (i in 1:nrow(res2)){
  print(i)
  chr = res2[i, 'ProbeChr']
  pos = res2[i, 'Probe_bp']
  gene = map[map$chr == chr & map$start<=pos & map$end>=pos, 'gene']
  res2[i, 'gene'] = paste(gene, collapse=';')
}


write.csv(res2, paste0(path, '../../collect/smr_ttt.csv'), row.names = F)
