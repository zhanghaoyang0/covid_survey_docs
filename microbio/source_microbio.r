library(dplyr)
library(stringr)
library(data.table)
complement = function(string){
    string = strsplit(string, '')[[1]]
    if (sum(string%in%c('A', 'T', 'C', 'G')==F)>0){return(NA); break}
    comp = sapply(string, function(x){switch(x, "A" = "T", "C" = "G", "T" = "A", "G" = "C", return(NA))})
    comp = paste(comp, collapse='')
    return(comp)
}

prs_filter_gwas = function(df, bim){
    info = df[,c('SNP', 'A1', 'A2')]%>%merge(bim[,c('SNP', 'A1', 'A2')], 'SNP')
    info = info%>%mutate(A1.z=sapply(A1.x, complement), A2.z=sapply(A2.x, complement))
    match = info%>%filter(A1.x == A1.y & A2.x== A2.y)%>%pull(SNP)
    rmatch = info%>%filter(A1.x == A2.y & A2.x== A1.y)%>%pull(SNP)
    cmatch = info%>%filter(A1.z == A1.y & A2.z== A2.y)%>%pull(SNP)
    out = df%>%filter(SNP%in%c(match, rmatch, cmatch))%>%select(SNP, A1, BETA, P)
    return(out)
}

h2l_R2 <- function(k, r2, p) {
    x = qnorm(1-k); z = dnorm(x); i = z/k
    C = k*(1-k)*k*(1-k)/(z^2*p*(1-p))
    theta = i*((p-k)/(1-k))*(i*((p-k)/(1-k))-x)
    h2l_R2 = C*r2 / (1 + C*theta*r2)
}

eas = list(); eur = list(); path_list = list()
eas[['data']] = '/home/yanglab_data/user/zhanghy/gwas/summary/microbiome/eas/'
eas[['bfile']] = '/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/chn'
eur[['data']] = '/home/yanglab_data/user/zhanghy/gwas/summary/microbiome/eur/'
eur[['bfile']] = '/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt'
path_list[['eas']] = eas; path_list[['eur']] = eur