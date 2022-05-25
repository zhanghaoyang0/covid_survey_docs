info="/home/yanglab_data/user/zhanghy/project/temp/code/gwascatalog_info_20220508.csv"

datas=`cat $info | grep 33462482 | awk -v FS=',' '{print $3}'`
for data in $datas
do 
url="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90011001-GCST90012000/${data}/${data}_buildGRCh37.tsv.gz"
wget -c $url
done

datas=`cat $info | egrep "33208821|33462485" | awk -v FS=',' '{print $3}'`
for data in $datas
do 
url="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90007001-GCST90008000/${data}/${data}_buildGRCh37.tsv.gz"
wget -c $url
done



