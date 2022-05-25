#=====================================================================================
# bbj bsed 
#=====================================================================================
### merge eqtl
zcat chr1_cis_eqtl_mapping_nofilt_nomulti_with_alleles.txt.gz | head -1 > eQTL_6_5e-8
for i in {1..22}
do
echo $i
zcat "chr${i}_cis_eqtl_mapping_nofilt_nomulti_with_alleles.txt.gz" | awk '$8<5e-8' >> eQTL_6_5e-8 
done


### map with ukb chn bim
library(dplyr)
complement = function(string){
    string = strsplit(string, '')[[1]]
    if (sum(string%in%c('A', 'T', 'C', 'G')==F)>0){return(NA); break}
    comp = sapply(string, function(x){switch(x, "A" = "T", "C" = "G", "T" = "A", "G" = "C", return(NA))})
    comp = paste(comp, collapse='')
    return(comp)
}

path_bbj = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/'
pthres = '1e-5'

file_bfile = paste0('/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/chn')
bim = read.table(paste0(file_bfile, '.bim'))%>%rename(SNP=V2, CHR=V1, POS=V4, A1=V5, A2=V6)%>%select(-V3)
frq = read.table(paste0(file_bfile, '.frq'), header=1)

df_raw = read.table(paste0(path_bbj, 'raw/eQTL_6/Peripheral_blood/eQTL_6_', pthres), header=1)
df = df_raw%>%mutate(FDR=p.adjust(p.value, 'fdr'), 
    CHR = sapply(SNP, function(x){t = strsplit(x, ':')[[1]][1]; t = gsub('chr', '', t); return(t)}))%>%
    rename(A1=ALT, A2=REF)%>%select(-SNP)

info = df%>%merge(bim, by=c('CHR', 'POS'))
info = info%>%mutate(A1.z=sapply(A1.x, complement), A2.z=sapply(A2.x, complement))

match = info%>%filter(A1.x == A1.y & A2.x== A2.y)%>%pull(SNP)
rmatch = info%>%filter(A1.x == A2.y & A2.x== A1.y)%>%pull(SNP)
cmatch = info%>%filter(A1.z == A1.y & A2.z== A2.y)%>%pull(SNP)

info = info%>%mutate(A1.x=ifelse(SNP%in%cmatch, A1.z, A1.x), A2.x=ifelse(SNP%in%cmatch, A2.z, A2.x),
                    t.stat=ifelse(SNP%in%rmatch, -t.stat, t.stat))
info = info%>%rename(A1=A1.x, A2=A2.y)%>%select(CHR, POS, A1, A2, gene, beta, t.stat, p.value, SNP)

df_out = info%>%select(SNP, beta, t.stat, p.value, FDR)%>%rename('t-stat'=t.stat, 'p-value'=p.value)

write.table(df_out, paste0(path_bbj, 'clean/eqtl/bbj_blood_', pthres, '.txt'), sep='\t', quote=F, row.names=F)

## add in info in esi and epi
snp_info = bim%>%merge(frq%>%select(SNP, MAF), by='SNP')

esi = read.table(paste0(path_bbj, 'clean/eqtl/bbj_blood_', pthres, '.esi'))%>%rename(SNP=V2)
esi = esi%>%merge(snp_info, all.x=T)
esi = esi%>%select(CHR, SNP, V3, A1, A2, MAF)
write.table(esi, paste0(path_bbj, 'clean/eqtl/bbj_blood_', pthres, '.esi'), sep='\t', row.names=F, col.names=F, quote=F)


epi = read.table(paste0(path_bbj, 'clean/eqtl/bbj_blood_', pthres, '.epi'))

info%>%select(CHR, SNP, POS, A1, A2, )

### make bsed
path="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/clean/eqtl/"
pthres="1e-5"
smr --eqtl-summary mateQTL.txt --matrix-eqtl-format --make-besd --out mybesd 

/home/yanglab_data/user/zhanghy/soft_slurm/gcta_1.93.2beta/smr_Linux \
--eqtl-summary ${path}bbj_blood_${pthres}.txt \
--matrix-eqtl-format --make-besd --out ${path}bbj_blood_${pthres} 


## check smrsnp frq in two population
frq_eur = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt.frq', header=1)
frq1 = frq%>%filter(SNP%in%smr$topSNP)
frq2 = frq_eur%>%filter(SNP%in%smr$topSNP)

#=====================================================================================
# convert to cojo format
#=====================================================================================
## sprack
library(dplyr)
path_list = list('Case'='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/clean/',
    'spracklen'='/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/clean/',
    'out'='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/gcta/for_cojo/')

for (data in c('Case_control_37_f', 'Case_control_37_m', 'spracklen_nature2020_t2d_f', 'spracklen_nature2020_t2d_m')){
  print(data)
  path_in = path_list[sapply(names(path_list), function(x){grepl(x, data)})][[1]]
  df = read.table(paste0(path_in, data, '.txt.gz'), header=1)
  df = df%>%select(SNP, A1, A2, FRQ, BETA, SE, P, N)%>%rename(freq=FRQ, b=BETA, se=SE, p=P)
  write.table(df, paste0(path_list[['out']], data), row.names = F, sep = '\t', quote = F)
}

#=====================================================================================
# pre- and post-meta smr
#=====================================================================================
path_in="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/gcta/for_cojo/"
path_out="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/gcta/smr/"

sex='f'
path_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/chn_${sex}"
path_eqtl="/home/yanglab_data/user/zhanghy/gwas/smr_file/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense"
path_eqtl="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/clean/eqtl/bbj_blood_5e-8"

for data in spracklen_nature2020_t2d_${sex} Case_control_37_${sex} spracklen_nature2020_t2d_${sex}_and_Case_control_37_${sex}_mtag spracklen_nature2020_t2d_${sex}_and_Case_control_37_${sex}_mtag
do
/home/yanglab_data/user/zhanghy/soft_slurm/gcta_1.93.2beta/smr_Linux \
--bfile $path_bfile \
--gwas-summary $path_in$data \
--beqtl-summary $path_eqtl \
--diff-freq-prop 1 \
--peqtl-smr 1e-5\
--out $path_out\eqtl/chn/$data 
done

# compare
path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/gcta/smr/'
map = read.csv('/home/yanglab_data/user/zhanghy/db/source/gencode_summ_hg19.csv')
gene_unique = names(table(map$gene))[table(map$gene)==1]
map = map[map$gene %in% gene_unique,]

res = data.frame()
datas = c('Case_control_37_m', 'Case_control_37_f', 'spracklen_nature2020_t2d_m', 'spracklen_nature2020_t2d_f',
          'spracklen_nature2020_t2d_m_and_Case_control_37_m_mtag', 'spracklen_nature2020_t2d_f_and_Case_control_37_f_mtag')

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

res = na.omit(res)
res1 = res[res$p_HEIDI>=0.05 & res$nsnp_HEIDI>=10 & !is.na(res$p_SMR),]
res1 = res1[res1$p_SMR < 0.05/19250,]

res1$gene = NA
for (i in 1:nrow(res1)){
  print(i)
  chr = res1[i, 'ProbeChr']
  pos = res1[i, 'Probe_bp']
  gene = map[map$chr == chr & map$start<=pos & map$end>=pos, 'gene']
  res1[i, 'gene'] = paste(gene, collapse=';')
}

write.csv(res1, paste0(path, '../../collect/t2d_cat/smr.csv'), row.names = F)


