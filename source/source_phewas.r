library(dplyr)
library(PheWAS)

info = setNames(data.frame(matrix(c('TSC2', 16, 2097986, 2139492, 'HUWE1', 23, 53559057, 53713664, 'USP7', 16, 8985954, 9057763, 'TRIP12', 2, 230628553, 230787902,'FMR1', 23, 146993437, 147032645, 
'HACE1', 6, 105175969, 105307794, 'BRCA1', 17, 41196312, 41277381,'MKRN2', 3, 12598586, 12625212, 'CHD8', 14, 21853358, 21924282, 'FDXR', 17, 72858619, 72869119), ncol=4, byrow=T)), c('gene', 'chr', 'start', 'end'))
info = info%>%mutate(chr=as.numeric(chr), start=as.numeric(start), end=as.numeric(end))