### hist with line, used in rbmo
library(ggplot2)
library(dplyr)
setwd('/home/yanglab_data/user/zhanghy/project/slurm_gwas_code/')

df = read.csv('plot_demo/data/gest.csv')
outcome = 'gest'; var = 'p3_SO2_24h'
start = round(min(df[, var]), 1); end = round(max(df[, var]), 1)
newx = seq(start, end, by=(end-start)/20)
df_plot = as.data.frame(table(cut(df[,var], breaks=newx, include.lowest=T), df[,outcome]))
colnames(df_plot) = c('range', 'type', 'count')
mat = matrix(df_plot$count, ncol=2)
rate = (mat[,1]/(mat[,1]+mat[,2]))/rowSums(mat)

fill_level = c('1', '0')
colour = c("lightblue", "#E69F00")
df_line = as.data.frame(cbind(df_plot[df_plot$type==1, 1], rate))
colnames(df_line) = c('range', 'rate')

png('ttt.png')    
ggplot(df_plot) + 
    geom_col(aes(x=range, y=count, fill=factor(type, levels=fill_level)))  +
    geom_smooth(data=df_line,aes(x = range, y = rate*1000),  method = "lm", se= T, color = "black", size = 1.5) +
    scale_fill_manual(values=colour) +
    labs(x = '', y='') +
    theme(legend.title=element_blank(), axis.text.x = element_text(angle=90, hjust=1), legend.position = "none")+ 
    scale_y_continuous(sec.axis = sec_axis(~ . / 1000))
dev.off()