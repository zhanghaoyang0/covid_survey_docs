


## compare
## ctrl is normal, not t1d, t2d, gd
library(dplyr)
path = '/home/yanglab_data/user/zhanghy/gwas/ukb/phewas/'
load(paste0(path, 'data/for_phewas.rdata'))


edu = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/data/edu_etc.tab', header=1, nrows=5)

df = pheno%>%mutate(T1D=grepl('E10', traits), T2D=grepl('E11', traits), GD=grepl('nc_1221_nc', traits), batch=as.factor(batch), center=as.factor(center), whr=wc/hc)%>%filter(sex==0)
df = df%>%select(T1D, T2D, GD, age, height, sitting_height, weight, bmi, hc, wc, whr)
write.csv(df, '/home/yanglab_data/user/zhanghy/gwas/summary/gd_pheno.csv', row.names=F)

out = c()
for (i in c('T1D', 'T2D', 'GD')){
  for (j in c('height', 'sitting_height', 'weight', 'bmi', 'hc', 'wc', 'whr')){
    sub = df%>%select(all_of(i), all_of(j), age)%>%na.omit()
    case = sub[sub[,i]==T,]
    ctrl = sub[sub[,i]==F&sub$T1D==F&sub$T2D==F&sub$GD==F,]

    

    formula = as.formula(paste0(i,'~',j,'+age'))
    mod = glm(formula, rbind(case, ctrl), family='binomial')
    out = c(out, i, j, mean(case[,j]), sd(case[,j]), mean(case[,j]), sd(case[,j]), summary(mod)$coefficient[2, c(1, 2, 4)], nrow(sub))
  }
}

res = data.frame(matrix(out, ncol=6, byrow=T))%>%mutate(X3=as.numeric(X3), X4=as.numeric(X4), X5=as.numeric(X5))
colnames(res) = c('trait', 'var', 'b', 'se', 'p', 'n')
res%>%filter(p<0.05/nrow(res))




