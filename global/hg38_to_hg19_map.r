### hg38_to_hg19_map
#=====================================================================================
# 1000g: make hg 38 to hg 19 dict
#=====================================================================================
## hg19.bed
library(dplyr)
file_bim1 = '/home/yanglab_data/user/zhanghy/gwas/bfile/1000g/no_mhc/eur/eur_nomhc.bim'
file_bim2 = '/home/yanglab_data/user/zhanghy/gwas/bfile/1000g/no_mhc/eas/eas_nomhc.bim'
file_bed = '/home/yanglab_data/user/zhanghy/db/map/1000g/hg19.bed'

df1 = read.table(file_bim1)[,c(1, 4)]
df2 = read.table(file_bim2)[,c(1, 4)]
df = rbind(df1, df2)%>%distinct()
df[,1] = paste0('chr', df[,1])
df[,3] = as.integer(df[,2]+1) # non-integer will cause error in crossmap
df[,4] = 1:nrow(df)

write.table(df, file_bed, quote=F, row.names=F, col.names=F)

## crossmap: hg19 to hg38
/home/yanglab_data/user/zhanghy/soft_slurm/liftover/liftOver \
/home/yanglab_data/user/zhanghy/db/map/1000g/hg19.bed \
/home/yanglab_data/user/zhanghy/soft_slurm/crossmap/hg19ToHg38.over.chain.gz \
/home/yanglab_data/user/zhanghy/db/map/1000g/hg38_map.bed \
/home/yanglab_data/user/zhanghy/db/map/1000g/hg38_unmap.bed

## make dict for map
library(dplyr)
bed_hg38 = read.table('/home/yanglab_data/user/zhanghy/db/map/1000g/hg38_map.bed')[,c(1, 2, 4)]
bed_hg19 = read.table('/home/yanglab_data/user/zhanghy/db/map/1000g/hg19.bed')[,c(1, 2, 4)]

bim1 = read.table('/home/yanglab_data/user/zhanghy/gwas/bfile/1000g/no_mhc/eur/eur_nomhc.bim')[,c(2,1,4)]
bim2 = read.table('/home/yanglab_data/user/zhanghy/gwas/bfile/1000g/no_mhc/eas/eas_nomhc.bim')[,c(2,1,4)]
bim = rbind(bim1, bim2)%>%distinct()
colnames(bim) = c('SNP', 'CHR', 'POS')

colnames(bed_hg19) = c('CHR', 'POS', 'id')
colnames(bed_hg38) = c('CHR_38', 'POS_38', 'id')

map = bed_hg38 %>% merge(bed_hg19, by = 'id')
map$id = NULL

map$CHR = gsub('chr', '', map$CHR)
map$CHR_38 = gsub('chr', '', map$CHR_38)

map = map %>% merge(bim, by = c('CHR', 'POS'))

write.table(map, '/home/yanglab_data/user/zhanghy/db/map/1000g/map_hg38_to_hg19', row.names = F, quote = F, sep = '\t')
#=====================================================================================
# ukb: make hg 38 to hg 19 dict
#=====================================================================================
## hg19.bed
library(dplyr)
file_bim1 = '/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt.bim'
file_bim2 = '/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/chn.bim'
file_bed = '/home/yanglab_data/user/zhanghy/db/map/ukb/hg19.bed'

bim1 = read.table(file_bim1); bim2 = read.table(file_bim2)
bim = rbind(bim1, bim2)%>%mutate(t=paste0(V1, ':',V4))%>%distinct(t,.keep_all = TRUE)%>%select(-t)
bim = bim%>%rename(CHR=V1, POS=V4)%>%filter(!(CHR==6&POS<=33448354&POS>=28477797))

bed = bim[,c(1,4)]%>%mutate(CHR=paste0('chr', CHR), STOP=as.integer(POS+1), ID=1:nrow(bim))%>%mutate(CHR=ifelse(CHR=='chrX', 'chr23', CHR))
write.table(bed, file_bed, quote=F, row.names=F, col.names=F)

## crossmap: hg19 to hg38
path="/home/yanglab_data/user/zhanghy/db/map/ukb/"
/home/yanglab_data/user/zhanghy/soft_slurm/liftover/liftOver \
${path}hg19.bed \
/home/yanglab_data/user/zhanghy/soft_slurm/crossmap/hg19ToHg38.over.chain.gz \
${path}hg38_map.bed \
${path}hg38_unmap.bed

## make dict for map, replace chrx with chr23
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/db/map/ukb/'
bed_hg38 = read.table(paste0(path, 'hg38_map.bed'))[,c(1, 2, 4)]%>%rename(CHR_38=V1, POS_38=V2, id=V4)%>%mutate(CHR_38=ifelse(grepl('X', CHR_38), 'chr23', CHR_38))
bed_hg19 = read.table(paste0(path, 'hg19.bed'))[,c(1, 2, 4)]%>%rename(CHR=V1, POS=V2, id=V4)%>%mutate(CHR=ifelse(grepl('X', CHR), 'chr23', CHR))

bim_formap = bim[,c(2,1,4)]%>%rename(SNP=V2)

map = bed_hg38%>%merge(bed_hg19, by='id')%>%select(-id)
map = map%>%mutate(CHR=gsub('chr', '', CHR), CHR_38=gsub('chr', '', CHR_38))%>%merge(bim_formap, by = c('CHR', 'POS'))
map = map%>%mutate(CHR=ifelse(grepl('X', CHR), 23, CHR), CHR_38=ifelse(grepl('X', CHR), 23, CHR_38))

write.table(map, paste0(path, 'map_hg38_to_hg19'), row.names = F, quote = F, sep = '\t')