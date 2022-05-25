#=====================================================================================
# no-kinship_list
# 22021 kinship
# 22009 pc
#=====================================================================================
## get kinship id and kinship field
all_field = as.vector(unlist(read.table('/bigdat1/pub/DownloadUKB/Converted/Main1.tab', nrows=1)))
# kinship
keep_field = c('f.eid2', 'f.22021.0.0')
n_field = sapply(keep_field, function(x){which(all_field==x)})
n_field

## extract
cat /bigdat1/pub/DownloadUKB/Converted/Main1.tab | cut -f 1,10042,10000-10009 \
> /bigdat1/user/zhanghy/gwas/ukb/para/temp.tab

## extract list
path = '/bigdat1/user/zhanghy/gwas/ukb/para/'
df = read.table(paste0(path, 'temp.tab'), header=1)
# kinship
table(df[,'f.22021.0.0'])
out = df[!is.na(df[,'f.22021.0.0'])&df[,'f.22021.0.0']==0, 'f.eid2']
out = data.frame(cbind(out, out))
names(out) = c('FID', 'IID')
write.table(out, paste0(path, 'no-kinship_list'), row.names=F, quote=F)
#=====================================================================================
# pop list
#=====================================================================================
file_in="/bigdat1/user/zhanghy/gwas/ukb/data/icd10_clean.txt"
file_out_chn="/bigdat1/user/zhanghy/gwas/ukb/para/chn_list"
file_out_brt="/bigdat1/user/zhanghy/gwas/ukb/para/brt_list"

awk '$5==5 || $5=="ethnicity" {print $1, $2}' $file_in > $file_out_chn
awk '$5==1001 || $5=="ethnicity" {print $1, $2}' $file_in > $file_out_brt

#=====================================================================================
# snp list with > 0.8 info
#=====================================================================================
cd /bigdat1/user/zhanghy/gwas/ukb/info/snp_info/
  
touch temp

for i in {1..22} X XY
do
echo $i
cat $path\ukb_mfi_chr$i\_v3.txt | awk '$8>0.8 {print $2}' >> temp
done

sort temp | uniq -c > ../../para/snp_info0.8_list

rm -rf temp
#=====================================================================================
# qc 
# dup name, dup pos list
# 1 rs151120166 715142  A
# 1 rs151120166 715142  T
# 1 1:13289_CCT_C  13289 C
# 1 rs568318295  13289 T
#=====================================================================================
path_in = '/home/yanglab_data/user/xiuxh/bfile/ori/'
path_out = '/home/yanglab_data/user/zhanghy/gwas/ukb/para/dupsnp/'

for (i in c(1:22, 'X')){
  print(i)
  file_out = paste0(path_out, 'chr', i)
  if (file.exists(file_out)){
    next()
  }
  df = read.table(paste0(path_in, 'chr', i, '.bim'))[,c(2,4)]
  snp_dup = unique(df[duplicated(df[,1])|duplicated(df[,2]), 1])
  write.table(snp_dup, file_out, row.names=F, col.names=F, quote=F)
  print('done')
}

#=====================================================================================
# extract brt bfile, qc
# https://www.cnblogs.com/esctrionsit/p/13415050.html
#=====================================================================================
## geno, mind, maf, hwe
# qc_1.sh
#!/usr/bin/bash
pop="$1"
chr="$2"

file_bfile="/home/yanglab_data/user/xiuxh/bfile/ori/chr${chr}"
file_keep_snp="/home/yanglab_data/user/zhanghy/gwas/ukb/para/snp_info0.8_list"
file_exclude_snp="/home/yanglab_data/user/zhanghy/gwas/ukb/para/dupsnp/chr${chr}"
file_keep_sample="/home/yanglab_data/user/zhanghy/gwas/ukb/para/${pop}_list"
file_out="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/${pop}/qc1_bychr/chr${chr}_1"

/home/xiuxh/tool/plink2/plink2 \
--bfile $file_bfile \
--extract $file_keep_snp \
--keep $file_keep_sample \
--exclude $file_exclude_snp \
--mind 0.05 \
--geno 0.05 \
--maf 0.01 \
--hwe 1e-6 \
--make-bed \
--out $file_out

# run
pop="chn"
chr="3"
compute=`frq | sort -k1 | awk '{print $2}' | head -1`
sbatch -N1 -n1 -c10 -o chr${chr}_1.log -w compute7 qc_1.sh $pop $chr

#======================================
## rename id in bim with longname
## PLINK 1.9 does not scale well to length-80+ variant IDs
# move
pop="chn"
new_path="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/${pop}_qc/qc_1_bychr/ori_bim/"
mkdir $new_path

mv ${new_path}../*bim ${new_path}

# rename
for i in $(seq 1 22) X;do
echo $i
cat ${new_path}chr${i}_1.bim \
| awk -F '\t' -v OFS='\t' '{if($2 !~ /^rs/) {print $1,$1"_"$4,$3,$4,$5,$6} else {print $0}}' \
> ${new_path}../chr${i}_1.bim
done

#======================================
## collect sample to keep
pop = "chn"
path = paste0('/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/', pop, '_qc/qc_1_bychr/')

id = c()
for (i in c(1:22, 'X')){
  print(i)
  id = c(id, read.table(paste0(path, 'chr', i, '_1.fam'))[,1])
}

frq = table(id)
remove = names(frq)[frq!=23]
out = setNames(data.frame(cbind(remove, remove)), c('FID', 'IID'))

write.table(out, paste0(path, '../', pop, '_qcbychr_outlier.txt'), row.names=F, quote=F)

#======================================
## merge
## qc_1_merge.sh
#!/usr/bin/bash
pop="$1"

path="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/${pop}_qc/"
file_bfile="${path}/qc_1_bychr/chr1_1"
file_merge_list="${path}${pop}_qcbychr_bfilelist"
file_remove_sample="${path}${pop}_qcbychr_outlier.txt"
file_out="${path}${pop}_1"

/home/xiuxh/tool/plink/plink \
--bfile $file_bfile \
--make-bed \
--remove $file_remove_sample \
--merge-list $file_merge_list \
--out $file_out

# run
pop="chn"
sbatch -N1 -n1 -c10 -o qc_1_merge.log -w compute7 qc_1_merge.sh $pop 

cat qc_1_merge.log

#======================================
### checksex
## qc_2_checksex.sh
#!/usr/bin/bash
pop="$1"
file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/${pop}_qc/${pop}_1"

plink \
--bfile $file_bfile \
--check-sex  

# run
sbatch -N1 -n1 -c20 -o qc_2_checksex.log -w compute5 qc_2_checksex.sh chn

cat qc_2_checksex.log

# rename
pop="chn"
path="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/${pop}_qc/"

rm -rf ${path}plink.hh
mv ${path}plink.sexcheck ${path}${pop}_1.checksex

#======================================
### check het
## qc_2_prune.sh
#!/usr/bin/bash
pop="$1"
file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/${pop}_qc/${pop}_1"
file_exclude_snp="/home/yanglab_data/user/zhanghy/gwas/tutorial/GWA_tutorial-master/1_QC_GWAS/inversion.txt"
file_out="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/${pop}_qc/${pop}_1_indepSNP"

plink \
--bfile $file_bfile \
--exclude $file_exclude_snp \
--range \
--indep-pairwise 50 5 0.2 \
--out $file_out

pop="chn"
sbatch -N1 -n1 -c10 -o qc_2_prune.log -w compute8 qc_2_prune.sh $pop

cat qc_2_prune.log

## qc_2_checkhet.sh
#!/usr/bin/bash
pop="$1"
path="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/${pop}_qc/"
file_bfile=${path}${pop}_1
file_snp_extract=${path}${pop}_1_indepSNP.prune.in
file_out=${path}${pop}_1

plink \
--bfile $file_bfile \
--extract $file_snp_extract \
--het \
--out $file_out

# run
pop="chn"
sbatch -N1 -n1 -c10 -o qc_2_checkhet.log -w compute10 qc_2_checkhet.sh $pop

cat qc_2_checkhet.log
#======================================
### remove sex and het outlier

## sex, het outlier 
# sex
pop="chn"
path="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/${pop}_qc/"
file_sex_outlier="${path}${pop}_1_sexoutlier.txt"

grep "PROBLEM" ${path}${pop}_1.checksex| awk '{print$1,$2}'> $file_sex_outlier

# het
file_script="/home/yanglab_data/user/zhanghy/gwas/tutorial/GWA_tutorial-master/1_QC_GWAS/heterozygosity_outliers_list.R"
file_in="${path}${pop}_1.het"
file_het_outlier="${path}${pop}_1_hetoutlier.txt"

Rscript --no-save $file_script $file_in $file_het_outlier

# merge
file_outlier="${path}${pop}_1_outlier.txt"

cat <(echo FID IID) <(cat <(cat $file_het_outlier|awk '{print $1, $2}'|grep -v FID)\
                      <(cat $file_sex_outlier)|sort|uniq)\
> $file_outlier

wc -l *outlier*

## remove
## qc_2.sh
#!/usr/bin/bash
pop="$1"
path="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/${pop}_qc/"
file_bfile=${path}${pop}_1
file_outlier="${path}${pop}_1_outlier.txt"
file_out=${path}../${pop}1

plink \
--bfile $file_bfile \
--remove $file_outlier \
--make-bed \
--out $file_out 

pop="chn"
sbatch -N1 -n1 -c10 -o ${pop}_2.log -w compute10 qc_2.sh $pop

cat ${pop}_2.log
#=====================================================================================
# frq
#=====================================================================================
## frq.sh
#!/usr/bin/bash
pop="$1"
path="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/"
file_bfile=${path}${pop}
file_out=${path}${pop}

plink \
--bfile $file_bfile \
--freq \
--out $file_out 

pop="brt"
sbatch -N1 -n1 -c10 -o ${pop}_frq.log -w compute7 frq.sh $pop
cat ${pop}_frq.log
#=====================================================================================
# prune (for pca, bz using all snp would out of memory) 
# reason: https://speciationgenomics.github.io/pca/ 
# One of the major assumptions of PCA is that the data we use is indpendent. i.e. there are no spurious correlations among the measured variables. 
# This is obviously not the case for most genomic data as allele frequencies are linkage. So we need to prune.
# indep-pairwise, window=50 Kb. window step=5, r2 threshold=0.2(i.e. the threshold of linkage we are willing to tolerate)
#=====================================================================================
## prune.sh
#!/usr/bin/bash
pop="$1"
path="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/"
file_bfile=${path}${pop}
file_out=${path}${pop}

plink \
--bfile $file_bfile \
--indep-pairwise 50 5 0.2 \
--out $file_out 

pop="brt"
sbatch -N1 -n1 -c10 -o ${pop}_prune.log -w compute7 prune.sh $pop
cat ${pop}_prune.log


#=====================================================================================
# pca (with pruned snp) 
#=====================================================================================
## pca.sh
#!/usr/bin/bash
pop="$1"
file_bfile="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/${pop}"
file_pruned="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/${pop}.prune.in"
file_out="/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/${pop}"

/home/xiuxh/tool/plink2/plink2 \
--bfile $file_bfile \
--extract $file_pruned \
--pca approx 10 \
--out $file_out

pop="brt"
sbatch -N1 -n1 -c20 -o ${pop}_pca.log -w compute7 pca.sh $pop
cat ${pop}_pca.log


## unrelated bfile (if need)
file_bfile="/bigdat1/user/zhanghy/gwas/ukb/bfile/chn"
file_keep_sample="/bigdat1/user/zhanghy/gwas/ukb/para/no-kinship_list"
file_out="/bigdat1/user/zhanghy/gwas/ukb/bfile/chn_unrelated"

plink \
--bfile $file_bfile \
--keep $file_keep_sample \
--make-bed \
--out $file_out