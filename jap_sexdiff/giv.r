
#=====================================================================================
# select iv 
#=====================================================================================
library(dplyr)
source('/home/yanglab_data/user/zhanghy/project/temp/code/source/source_mr.r')
get_t2dcat_diff_info('all')
datas = unique(as.vector(pair))
datas = datas[grepl('0.8|sp|37', datas)]

for (pthres in c('1e-1', '1e-2', '1e-3', '1e-4', '1e-5')){
    print(pthres)
    for (data in datas){
        print(data)
        file_out = paste0(path_list[['bbj']], 'result/giv/clump/for_clump/', pthres, '/', data, '.txt')
        if (file.exists(file_out)){next}
        path_in = path_list[sapply(names(path_list), function(x){grepl(x, data)})][[1]]
        df_raw = read.table(paste0(path_in, 'clean/', data, '.txt.gz'), header=1, sep = '\t')
        df = df_raw%>%select(SNP, A1, A2, P)%>%filter(P<as.numeric(pthres))
        write.table(df, file_out, row.names=F, quote=F, sep='\t')
    }
}

### clump
pthres="1e-1"

for data in spracklen_nature2020_t2d Case_control_37 E11_0.8 cataract-adjE11_0.8
do for sex in m f 
do
file_gwas="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/giv/clump/for_clump/${pthres}/${data}_${sex}.txt"
file_out="/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/giv/clump/clump_out/${pthres}/${data}_${sex}"

if [[ $data = *"E11"*  ||  $data = *"cat"* ]]
then file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/brt_${sex}"
else file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/chn_${sex}"
fi
# clump
if [ ! -f $file_out.clumped ]
then 
plink \
    --bfile $file_bfile \
    --clump $file_gwas \
    --clump-p1 1 --clump-r2 0.05 --clump-kb 1000 \
    --clump-snp-field SNP --clump-field P \
    --out $file_out
fi 
awk 'NR!=1{print $3}' ${file_out}.clumped >  "${file_out}_clumped-snp"
# make geno
plink \
--recode A \
--bfile $file_bfile \
--extract ${file_out}_clumped-snp \
--out $file_out
done;done

#=====================================================================================
# make geno file 
#=====================================================================================
collect = data.frame()
for (i in 1:22){
    print(i)
    df = read.table(paste0(i, '.hg19_multianno.txt'), fill=T, header=1)
    collect = rbind(collect, df)
}


## prepare genotype
# recodeA make genotype to 0,1,2 (number of MA), recodeD make genotype to 0,1 (het or not), recodeAD make genotype to both col
pop="chn"
file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/${pop}"
file_snp="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/prune_1000kb/${pop}.prune.in"
file_out="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/prune_1000kb/${pop}"






source('/home/yanglab_data/user/zhanghy/project/temp/code/sex_diff/giv_source.r')

pheno_raw = read.table(paste0(path_ukb, 'data/icd10_body_self-report_blood_clean.txt'), header=1)%>%
    mutate(cat=ifelse(grepl('H25|H26|H28', icd), 1, 0), t2d=ifelse(grepl('E11', icd), 1, 0))%>%
    select(FID, sex, age, ethnicity, t2d, cat, center, batch)
pop = 'eas'; sex = 'm'

# geno first col
file_geno = paste0(path_ukb, 'bfile/prune_1000kb/', ifelse(pop=='eas', 'chn', 'brt'), '.raw')
geno = read.table(pipe(paste0("awk '{print $1}' ", file_geno)), header=1)

# filter pheno to geno fid
set.seed(1)
pheno = pheno_raw%>%filter(ethnicity==ifelse(pop=='eas', 5, 1001), sex==ifelse(sex=='m', 1, 0), FID%in%geno$FID)
pheno = pheno%>%mutate(group=sample_3set(pheno, 'cat'))%>%arrange(group)%>%select(FID, age, t2d, cat, group)
table(pheno$group, pheno$cat)

# map geno to pheno fid
map_row = c()
for (i in 1:nrow(pheno)){map_row = c(map_row, which(geno$FID==pheno[i, 'FID']))}
geno = data.frame(FID=geno[map_row, 'FID'])
table(geno$FID%in%pheno$FID)

# geno by block
ncol = as.numeric(system(paste0("head -1 ", file_geno, "| awk '{print NF}'"), intern = TRUE)) # from col 7 to ncol
# blocks = c(seq(7, ncol, 5), ncol)
# command = paste("awk '{for(i=", blocks[i], ';i<=', blocks[i+1], ';i++) printf $i" "; print ""}', "'", file_geno)

command = paste("awk '{for(i=1;i<=", ncol, ';i++) printf $i" "; print ""}', "'", file_geno)
markers = as.matrix(read.table(pipe(command), header=1)[map_row,])
markers[is.na(markers)] = 0 # impute by zhanghy, no sure if is right

allocate_mat(pheno) 
trait1 = 't2d'; trait2 = 'cat'
T = pheno[,trait1]; y = pheno[,trait2]; age = pheno$age

res <- foreach(marker=iter(markers,by='col'), .combine=rbind) %dopar% {
    # Split the sample in two: #
    X1 = cbind(ones1,marker[1:break1]) # Regressor matrix for first half
    X2 = cbind(ones2,marker[(break1+1):break2]) # Regressor matrix for second
    ## Normal GWAS ##
    beta1_y = fols(y[1:break1],X1)  # OLS on half of the sample
    beta2_y = fols(y[(break1+1):break2],X2)  # Other half
    beta_y_full = fols(y[1:break2],rbind(X1,X2)) # OLS on the full sample
    beta_T = fols(T[(break1+1):break2],X2) # GWAS for T on the second sample
    ## Conditional GWAS ##
    beta1_y_cond =  fols(y[1:break1],cbind(X1,T[1:break1], age[1:break1]))  # OLS conditioning on T
    beta2_y_cond =  fols(y[(break1+1):break2],cbind(X2,T[(break1+1):break2], age[(break1+1):break2]))
    beta_y_cond_full = fols(y[1:break2],cbind(rbind(X1,X2),T[1:break2], age[1:break2]))
    # only return the coefficients for the markers, not the constant
    return(c(beta1_y[2],beta2_y[2],beta_T[2],beta1_y_cond[2], beta2_y_cond[2], beta_y_full[2], beta_y_cond_full[2]))
}

save(res, file=paste0(path_bbj, 'result/giv/', pop, '_', sex, '_', trait1, '&', trait2, '_res.rdata'))