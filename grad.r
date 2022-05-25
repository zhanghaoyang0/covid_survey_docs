
#=====================================================================================
# dis_info
#=====================================================================================
library(dplyr)
drop_data = c('Case_control_2', 'Case_control_8', 'Case_control_49', 'Case_control_15', 'Case_control_55', 
  'Case_control_95', 'Case_control_90','Case_control_99', 'Case_control_105', 'Case_control_84')

path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/'
ldsc = read.csv(paste0(path, 'result/collect/ldsc_t2d.csv'))%>%filter(!(data1%in%drop_data|data2%in%drop_data))

h2 = read.csv(paste0(path, 'para/h2.txt'), sep='\t')
map = read.csv(paste0(path, 'para/para.csv'))%>%select(id, case, ctrl, Main.group, Disease.type)%>%rename(data=id)
map[map=='Colorectal Cancer_1'] = 'Colorectal Cancer'

datas = ldsc%>%filter(data1=='Case_control_14'|data2=='Case_control_14')%>%select(data1, data2)
datas = unique(unlist(datas))
diseases = datas[!grepl('QTL', datas)]
h2 = h2%>%filter(data%in%paste0(diseases))
nsnp1 = ldsc%>%filter(data1=='Case_control_14')%>%select(data2, nsnp)%>%rename(data=data2)%>%distinct()
nsnp2 = ldsc%>%filter(data2=='Case_control_14')%>%select(data1, nsnp)%>%rename(data=data1)%>%distinct()
nsnp = rbind(nsnp1, nsnp2)

# dis_info
info = h2%>%filter(!data%in%drop_data)%>%merge(nsnp, 'data')%>%merge(map, 'data')%>%mutate(prop=case/(case+ctrl))
write.csv(info, paste0(path, 'ttt.csv'))

# # t2d h2
# cat Case_control_36\&spracklen_nature2020_t2d.log 
# 0.1546 (0.0135)  2.30175e-30

# cat Case_control_14\&Case_control_36.log
# 0.0973 (0.0088) 2.031777e-28

 (0.0035)
0.0073 (0.0025)


2*pnorm(abs(0.0338/0.0028), lower.tail=F)
2*pnorm(abs(0.0073/0.0025), lower.tail=F)
#=====================================================================================
# disease_ldsc plot
#=====================================================================================
library(ggplot2)
library(viridis)
library(scales)
library(dplyr)

# make t2d as data1
df1 = ldsc%>%filter(data1=='Case_control_14')%>%mutate(data=data2)
df2 = ldsc%>%filter(data2=='Case_control_14')%>%mutate(data=data1, h2_2=h2_1, h2_2_se=h2_1_se, h2_2_p=h2_1_p)
df = rbind(df1, df2)%>%select(-h2_1, -h2_1_se, -h2_1_p, -data1)
df = df%>%merge(map, by='data')%>%rename(category=Disease.type)%>%rename(label=Main.group)
write.csv(df, paste0(path, 'ttt.csv'))

df1 = df%>%filter(constrain==0)
df2 = df%>%filter(constrain==1)

sub1 = cbind('ldsc', df1$label, df1$category, df1$rg, df1$rg-1.96*df1$rg_se, df1$rg+1.96*df1$rg_se, df1$rg_p)
sub2 = cbind('ldsc_constrain', df2$label, df2$category, df2$rg, df2$rg-1.96*df2$rg_se, df2$rg+1.96*df2$rg_se, df2$rg_p)
sub3 = cbind('h2_2', df1$label, df1$category, df1$h2_2, df1$h2_2-1.96*df1$h2_2_se, df1$h2_2+1.96*df1$h2_2_se, df1$h2_2_p)
sub4 = cbind('h2_2_constrain', df2$label, df2$category, df2$h2_2, df2$h2_2-1.96*df2$h2_2_se, df2$h2_2+1.96*df2$h2_2_se, df2$h2_2_p)

ldsc_p = setNames(data.frame(rbind(sub1, sub2)), c('method', 'trait', 'Category', 'rg', 'rg_l', 'rg_u', 'p'))

ldsc_p[, c('rg', 'rg_l', 'rg_u', 'p')] = sapply(ldsc_p[, c('rg', 'rg_l', 'rg_u', 'p')], as.numeric)
ldsc_p$trait = factor(ldsc_p$trait, levels = unique(df1[order(df1$rg, decreasing = T), 'label']))

ldsc_p$method = factor(ldsc_p$method, levels = c('ldsc', 'ldsc_constrain'))

mr_trait = unique(unlist(ldsc%>%filter(constrain==0&rg_p<0.025)%>%select(data1, data2)))
mr_trait = mr_trait[mr_trait!='Case_control_14']

## ldsc pair
pair = rbind(cbind('Case_control_14', mr_trait),cbind(mr_trait, 'Case_control_14'))
pair = setNames(data.frame(pair), c('data1', 'data2'))
write.table(pair, paste0(path, 'para/mr_pair.txt'), row.names = F, sep = '\t', quote = F)

#=====================================================================================
# biomarker ldsc plot
#=====================================================================================
{
  library(ggplot2)
  library(viridis)
  library(scales)
  library(dplyr)
  source('/home/yanglab_data/user/zhanghy/project/temp/code/source_mr.r')
  path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/'

  map = get_bloodmarker('')
  markers = get_bloodmarker('marker')

  df = read.csv(paste0(path, 'result/collect/ldsc_blood.csv'))
  df = df[(df$data1=='spracklen_nature2020_t2d' & df$data2%in%markers) | (df$data2=='spracklen_nature2020_t2d' & df$data1%in%markers),]
  df = df %>% merge(map, by.x = 'data2', by.y = 'marker')
  df[,c('rg', 'rg_se', 'rg_p')] = sapply(df[,c('rg', 'rg_se', 'rg_p')], as.numeric)

  # get marker label
  for (i in map$marker){
    df[df==i] = map[map$marker==i, 'label']
  }

  df1 = df[df$constrain==0&df$liability==1,]
  df2 = df[df$constrain==1&df$liability==1,]

  sub1 = cbind('ldsc', df1$label, df1$category, df1$rg, df1$rg-1.96*df1$rg_se, df1$rg+1.96*df1$rg_se, df1$rg_p)
  sub2 = cbind('ldsc_constrain', df2$label, df2$category, df2$rg, df2$rg-1.96*df2$rg_se, df2$rg+1.96*df2$rg_se, df2$rg_p)

  ldsc = setNames(data.frame(rbind(sub1, sub2)), c('method', 'trait', 'Category', 'rg', 'rg_l', 'rg_u', 'p'))

  ldsc[, c('rg', 'rg_l', 'rg_u', 'p')] = sapply(ldsc[, c('rg', 'rg_l', 'rg_u', 'p')], as.numeric)
  ldsc$trait = factor(ldsc$trait, levels = unique(df1[order(df1$rg, decreasing = T), 'label']))

  ldsc$method = factor(ldsc$method, levels = c('ldsc', 'ldsc_constrain'))
  mr_trait = as.character(ldsc[ldsc$method=='ldsc'&ldsc$p<0.025, 'trait'])

  label = c()
  pos = c()
  category = c()

  measures = c('ldsc', 'ldsc_constrain')
  levels(ldsc$method) = measures

  for (method in measures){
    temp = ldsc[ldsc$method==method,]
    pos = c(pos, sapply(levels(ldsc$trait), function(x){return(temp[temp$trait == x, 'rg_u'])}))
    category = c(category, sapply(levels(ldsc$trait), function(x){return(temp[temp$trait == x, 'Category'])}))
    
    label = c(label, sapply(levels(ldsc$trait), function(x){if(temp[temp$trait == x, 'p'] <= 0.001/2){return('***')}
      else if(temp[temp$trait == x, 'p'] <= 0.01/2){return('**')}
      else if(temp[temp$trait == x, 'p'] <= 0.05/2){return('*')}
      else{return('')}}))
  }

  label_star = data.frame(method=as.vector(sapply(measures, function(x){rep(x, nrow(temp))})), 
                          trait = c(rep(levels(ldsc$trait), 4)), label = label, Category = category, rg = pos + 0.05)
  label1 = label_star%>%filter(label=='*')
  label2 = label_star%>%filter(label=='**')
  label3 = label_star%>%filter(label=='***')

  pdf(paste0(path, 'ttt.pdf'), width=10, height=7)

  scaleFUN <- function(x) sprintf("%.1f", x) # digit place
  ggplot(ldsc, aes(x=trait, y=rg, color = Category)) +  
    geom_pointrange(aes(ymin=rg_l, ymax=rg_u), size = 0.4) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.spacing = unit(2, "lines"),
          strip.background = element_blank(), strip.text.y = element_blank()) + # legend.position = "none"
    labs(x ="Biomarker") +  ylab('')  + 
    facet_grid(vars(method), scales = "free") +
    scale_y_continuous(labels=scaleFUN) + 
    geom_text(data = label1,label = '*') +
    geom_text(data = label2,label = '**') +
    geom_text(data = label3,label = '***') +
    theme(axis.text.x=element_text(size=12))

  dev.off()
}

#=====================================================================================
# disease mr
#=====================================================================================
col_name = c('method', 'data1', 'data2', 'nsnp', 'b', 'b_l', 'b_u', 'p')
path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/'

# lcv
df = read.csv(paste0(path, 'result/collect/lcv.csv'))
df1 = df[df[,2]%in%mr_trait & df[,1] == 'Case_control_14',]
df2 = df%>%filter(data2=='Case_control_14'&(data1=='Case_control_1'|data1=='Case_control_11'|data1=='Case_control_103'))
df2 = df2%>%mutate(gcp.pm=-gcp.pm, t=data1, data1=data2, data2=t)%>%select(-t)
lcv = rbind(df1, df2)%>%mutate(method='LCV',nsnp=lcv_snp_n, b=gcp.pm, b_l=gcp.pm-1.96*gcp.pse, b_u=gcp.pm+1.96*gcp.pse, p=lcv_p)%>%
  select(method, data1, data2, nsnp, b, b_l, b_u, p)%>%distinct()

# gsmr
df = read.csv(paste0(path, 'result/collect/gsmr_t2d.csv'))
gsmr = df%>%mutate(method='GSMR',b=bxy, b_l=bxy-1.96*bxy_se, b_u=bxy+1.96*bxy_se, p=bxy_pval)%>%
  select(method, data1, data2, nsnp, b, b_l, b_u, p)%>%distinct()

# mr
df = read.csv(paste0(path, 'result/collect/mr_t2d.csv'))
mr = data.frame()
for (method in c('MR.Egger', 'Inverse.variance.weighted', 'Weighted.median', 'Weighted.mode')){
  print(method)
  nsnp = paste0(method, '_nsnp'); b = paste0(method, '_b'); se = paste0(method, '_se'); p = paste0(method, '_pval')
  sub = setNames(cbind(method, df$data1, df$data2, df[nsnp], df[b], df[b]-1.96*df[se], df[b]+1.96*df[se], df[p]), col_name)
  mr = rbind(mr, sub)
}

# cause
df = read.csv(paste0(path, 'result/collect/cause_t2d.csv'))
coef = matrix(as.numeric(unlist(strsplit(gsub('\\(|,|\\)', '', df$gamma_causal), ' '))), ncol=3, byrow=T)
z = coef[,1]/((coef[,3]-coef[,2])/(2*1.96))
p = pnorm(abs(z), lower.tail=FALSE)*2
cause = setNames(data.frame(cbind('CAUSE', df$data1, df$data2, df$nsnp, coef, p)), col_name)

res = data.frame(rbind(lcv, cause, gsmr, mr))
res[, c('b', 'b_l', 'b_u', 'p')] = sapply(res[, c('b', 'b_l', 'b_u', 'p')], as.numeric)


res[,c('or', 'or_l', 'or_u')] = exp(res[,c('b', 'b_l', 'b_u')])
res[res$method == 'LCV', c('or', 'or_l', 'or_u')] = res[res$method == 'LCV', c('b', 'b_l', 'b_u')] # b is gcp

res[res == 'MR.Egger'] = 'MR-Egger'
res[res == 'Inverse.variance.weighted'] = 'IVW'
res[res == 'Weighted.median'] = 'Weighted median'
res[res == 'Weighted.mode'] = 'Weighted mode'


write.csv(res, paste0(path, 'result/collect/collect_t2d.csv'), row.names=F)
#=====================================================================================
# disease mr heatmap
#=====================================================================================
require(ggplot2) # ggplot2
require(reshape2) # melt data
library(forestplot)
path = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/'

convert_p_notation = function(x){
  n1 = round(as.numeric(strsplit(x, 'e-')[[1]][1]), 2)
  n2 = round(as.numeric(strsplit(x, 'e-')[[1]][2]), 0)
  if (is.na(n2)){
    out = n1
  } else {out = as.expression(bquote(.(n1)~x~10^-.(n2)))} # forest plot only accpet expression
  return(out)
}

df = read.csv(paste0(path, 'result/collect/collect_t2d.csv'), stringsAsFactors = F)
df = df%>%mutate(data=ifelse(data1=='Case_control_14', data2, data1))%>%merge(map)%>%rename(label=Main.group, category=Disease.type)

## heatmap | mr
methods = c("CAUSE", 'GSMR', 'IVW', 'MR-Egger', "Weighted median", "Weighted mode")
labels = unique(df$label)
mat = matrix(nrow = length(methods)*2, ncol = length(labels))
colnames(mat) = labels
methods = c(methods, paste0(methods, '_re'))
rownames(mat) = methods 
mat_p = mat

for (i in methods){
  for (j in labels){
    if (!grepl('_re', i)){
      b = df%>%filter(data1=='Case_control_14'&label==j&method==i)%>%pull(or)
      p = df%>%filter(data1=='Case_control_14'&label==j&method==i)%>%pull(p)
    } else {
      b = df%>%filter(data2=='Case_control_14'&label==j&df$method==gsub('_re', '', i))%>%pull(or)
      p = df%>%filter(data2=='Case_control_14'&label==j&df$method==gsub('_re', '', i))%>%pull(p)
    }
    mat[rownames(mat)==i, colnames(mat)==j] = b
    mat_p[rownames(mat_p)==i, colnames(mat_p)==j] = p
  }
}

# sort
temp = df[,c('label','category')]%>%distinct%>%arrange(desc(category))
# idx = order(mat_p[row.names(mat_p)=='IVW',], decreasing=T)
mat = mat[, temp$label]
mat_p = mat_p[, temp$label]
lcv_label = colnames(mat) = colnames(mat_p) = paste0( temp$label, ' (', temp$category, ')')

text = melt(mat)
text_p = melt(mat_p)
text_label = sapply(1:nrow(text_p), function(x){if(text_p[x, 3] <= 0.001/13){return('***')}
  else if(text_p[x, 3] <= 0.01/13){return('**')}
  else if(text_p[x, 3] <= 0.05/13){return('*')}
  else{return('')}})
text.label=round(text$value,2)

pdf(paste0(path, 'ttt.pdf'), width=10, height=8)
ggplot(text, aes(x=Var1, y=Var2, fill=value,label=text_label))+ 
  geom_tile() + 
  geom_text() +
  labs(x="",y="",fill = 'OR') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size=12)) +
  scale_fill_distiller(palette="Spectral", guide = "legend") +
  scale_x_discrete(breaks = c(row.names(mat)), labels = rep(row.names(mat)[1:6], 2))
dev.off()
        

## heatmap | lcv
methods = c("LCV", '1', '2', '3', "4", "5")

mat = matrix(nrow = length(methods)*2, ncol = length(labels))
colnames(mat) = lcv_label
methods = c(methods, paste0(methods, '_re'))
rownames(mat) = methods 
mat_p = mat

for (i in methods){
  for (j in lcv_label){
    if (i == 'LCV'){
      b = df%>%filter(data1=='Case_control_14'&paste0(label, ' (', category, ')')==j&method==i)%>%pull(b)
      p = df%>%filter(data1=='Case_control_14'&paste0(label, ' (', category, ')')==j&method==i)%>%pull(p)
    } else {
      or = 0
      p = 1
    }
    mat[rownames(mat)==i, colnames(mat)==j] = b
    mat_p[rownames(mat_p)==i, colnames(mat_p)==j] = p
  }
}

text = melt(mat)
text_p = melt(mat_p)
text_label = sapply(1:nrow(text_p), function(x){if(text_p[x, 3] <= 0.001/13){return('***')}
  else if(text_p[x, 3] <= 0.01/13){return('**')}
  else if(text_p[x, 3] <= 0.05/13){return('*')}
  else{return('')}})
text.label=round(text$value,2)

# pdf: 8*6
pdf(paste0(path, 'ttt1.pdf'), width=10, height=8)
ggplot(text, aes(x=Var1, y=Var2, fill=value,label=text_label))+ 
  geom_tile() + 
  geom_text() +
  labs(x="",y="",fill = 'GCP') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_distiller(palette="Spectral", guide = "legend") +
  scale_x_discrete(breaks = c(row.names(mat)), labels = rep(row.names(mat)[1:6], 2))
dev.off()



## forest
for (trait in c('Cl', 'K')){
  sub = df[(df$data1==trait | df$data2==trait) & df$method != 'LCV',]
  
  sub = sub[order(factor(sub$method, levels=c("CAUSE", 'GSMR', 'IVW', 'MR-Egger', "Weighted median", "Weighted mode"))),]
  sub = sub[order(sub$data1, decreasing = T),]

  tabletext=list(list('Method', 'CAUSE', 'GSMR', 'IVW', 'MR-Egger', 'Weighted median', 'Weighted mode'))
  
  legend = c(paste0('T2D to ', trait), paste0(trait, ' to T2D'))
  
  # x_min = ceiling(min(sub$b_l)*10)/10 - 0.1
  # x_max = floor(max(sub$b_u)*10)/10 + 0.1
  # step = 0.2
  if (trait=='Cl'){
    x_min = -0.6
    x_max = 0.6
    step = 0.2
  } else{
    x_min = -0.2
    x_max = 0.3
    step = 0.1
  }
  
  png(paste0('D:/nutstore/project/gwas/paper/t2d_blood/plot/forest/', trait, '_forest.png'), width=900, height=600, res=120)
  
  forestplot(tabletext,
             mean  = cbind(c(NA, sub$b[1:6]), c(c(NA, sub$b[7:12]))),
             lower = cbind(c(NA, sub$b_l[1:6]), c(c(NA, sub$b_l[7:12]))),
             upper = cbind(c(NA, sub$b_u[1:6]), c(c(NA, sub$b_u[7:12]))),
             new_page = T, boxsize=0.1, line.margin=0.25, is.summary=c(TRUE, rep(FALSE,6)), 
             cex=15.5, vertices = TRUE, xlog=FALSE, zero=0, clip=c(x_min, x_max), 
             legend = legend,
             # title = "B.                                                                                                                                                                                                                                                                                                                       ",
             xlab="Causal Effect (Liability β)", col=fpColors(box=c("blue", "red"), line="darkblue", summary="royalblue"),
             legend_args = fpLegend(pos = list(x=0.8, y=0.9), gp=gpar(col="#CCCCCC", fill="#F9F9F9")),
             xticks = seq(from = x_min, to = x_max, by = step),
             txt_gp = fpTxtGp(ticks = gpar(fontfamily = "", cex=1),
                              xlab  = gpar(fontfamily = "", cex = 1),
                              title = gpar(fontfamily = "", cex=1, align='l'),
                              legend = gpar(fontfamily = "", cex = 1)))
  
  dev.off()
}

#=====================================================================================
# eas female t2d bmi
#=====================================================================================
mr = data.frame()
load('/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/mr/spracklen_nature2020_t2d_f&QTL_4_f.rdata')
load('/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/mr/spracklen_nature2020_t2d_f&QTL_4_f.rdata')

df = read.csv(paste0(path, 'result/collect/t2d_cat/mr.csv'))
mr = data.frame()
for (method in c('MR.Egger', 'Inverse.variance.weighted', 'Weighted.median', 'Weighted.mode')){
  sub = df[,c(names(df)[grepl(method, names(df))], 'nsnp_x', 'data1', 'data2', 'snp_r2.exposure', 'snp_r2.outcome', 'pthres')]%>%mutate(t=method)
  names(sub) = c('nsnp', 'b_raw', 'se', 'p', 'nsnp_x', 'data1', 'data2', 'r2_data1', 'r2_data2', 'pthres', 'method')
  mr = rbind(mr, sub)
}
