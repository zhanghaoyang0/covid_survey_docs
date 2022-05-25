#=====================================================================================
# coloc
# s= the proportion of samples are cases, when type=”cc”. It also needs the minor allele frequency.
#=====================================================================================
# library(remotes)
# install_github("chr1swallace/coloc",build_vignettes=TRUE)
library(coloc)
library(dplyr)
library(readxl)
library(locuscomparer)

file_eqtl = '/home/yanglab_data/user/zhanghy/gwas/summary/eqtl/clean/GTEx_v8_FDXR.txt'
file_info = '/home/yanglab_data/user/zhanghy/gwas/summary/gwascatalog/info/gwascatalog_info.xlsx'
info = data.frame(read_excel(file_info, sheet=1))
eqtl = read.table(file_eqtl, header=T)
files =  list.files('/home/yanglab_data/user/zhanghy/gwas/summary/gwascatalog/clean/')

out = c()
for (file in files){
  print(file)
  file_gwas = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/gwascatalog/clean/', file)
  trait = gsub('.txt', '', file)
  
  # gwas = read.table(file_gwas, header=T)
  # input = merge(eqtl, gwas, by="SNP", all=FALSE, suffixes=c("_eqtl","_gwas"))
  # result = coloc.abf(dataset1=list(pvalues=input$P_gwas, type="quant", N=input$N_gwas),
  #                    dataset2=list(pvalues=input$P_eqtl, type="quant", N=input$N_eqtl), MAF=input$FRQ_eqtl)
  # result = coloc.abf(dataset1=list(pvalues=input$P_gwas, type="cc", s=0.33, N=input$N_gwas),
  #                    dataset2=list(pvalues=input$P_eqtl, type="quant", N=input$N_eqtl), MAF=input$FRQ_eqtl)
  # result$data = trait
  # save(result, file=file_out)

  file_out = paste0('~/', trait, '.rdata')
  load(paste0('~/', trait, '.rdata'))
  
  out = c(out, trait, result$summary)
  
  # need_result = result$results %>% filter(SNP.PP.H4 > 0.75)
  # print(need_result)
}

res = data.frame(matrix(out, ncol=7, byrow=T))
names(res) = c('trait', names(result$summary))

res$trait_mean = res$trait
for (i in c(unique(c(res$trait)))){
  res[,c('trait_mean')][res[,c('trait')]==i] = info[info$label==i, 'trait']
}

write.csv(res, paste0('~/fdxr_coloc_collect.csv'), row.names=F)



info = data.frame(read_excel(file_info, sheet=1))

## focal and plot

file = 'GCST008055.txt'
file_gwas = paste0('/home/yanglab_data/user/zhanghy/gwas/summary/gwascatalog/clean/', file)
trait = 'GCST008055'

gwas = read.table(file_gwas, header=T)
input = merge(eqtl, gwas, by="SNP", all=FALSE, suffixes=c("_eqtl","_gwas"))
result = coloc.abf(dataset1=list(pvalues=input$P_gwas, type="quant", N=input$N_gwas),
                   dataset2=list(pvalues=input$P_eqtl, type="quant", N=input$N_eqtl), MAF=input$FRQ_eqtl)
need_result = result$results %>% filter(SNP.PP.H4 > 0.75)
need_result
'rs2121646'

gwas1 = gwas[gwas$CHR==17&
               gwas$POS>eqtl[1,'POS']-50*1000&
               gwas$POS<eqtl[nrow(eqtl),'POS']+50*1000,]

write.table(gwas1, '~/GCST008055_FDXR.txt', row.names=F, quote=F, sep='\t')


gwas_fn="~/GCST008055_FDXR.txt"
eqtl_fn="/home/yanglab_data/user/zhanghy/gwas/summary/eqtl/clean/GTEx_v8_FDXR.txt"
marker_col="SNP"
pval_col="P"


png('~/t.png')

locuscompare(in_fn1=gwas_fn, in_fn2=eqtl_fn, title1="GWAS", title2="eQTL", snp='rs2121646',
             marker_col1= marker_col, pval_col1=pval_col, marker_col2=marker_col, pval_col2=pval_col)

dev.off()


png('~/t1.png')

locuscompare(in_fn1=gwas_fn, in_fn2=eqtl_fn, title1="GWAS", title2="eQTL", snp='rs690514',
             marker_col1= marker_col, pval_col1=pval_col, marker_col2=marker_col, pval_col2=pval_col)

dev.off()
