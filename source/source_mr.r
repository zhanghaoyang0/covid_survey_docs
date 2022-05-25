# source.r
library(dplyr)
library(TwoSampleMR)
options(stringsAsFactors = F)

get_gsmr_para = function(x, pop='eas'){
  if (x == 'gsmr'){
    if (pop == 'eas'){
      n_ref <<- 481    # Sample size of the 1000g eas (nrow of fam)
    } else if (pop == 'eur') {n_ref <<- 489} 
    else {print('need pop'); 
      break()}
    gwas_thresh <<- 5e-8    # GWAS threshold to select SNPs as the instruments for the GSMR analysis
    single_snp_heidi_thresh <<- 0.01    # p-value threshold for single-SNP-based HEIDI-outlier analysis | default is 0.01
    multi_snp_heidi_thresh <<- 0.01    # p-value threshold for multi-SNP-based HEIDI-outlier analysis | default is 0.01
    nsnps_thresh <<- 1   # the minimum number of instruments required for the GSMR analysis | default is 10
    heidi_outlier_flag <<- T    # flag for HEIDI-outlier analysis
    ld_r2_thresh <<- 0.05    # LD r2 threshold to remove SNPs in high LD 
    ld_fdr_thresh <<- 0.05   # FDR threshold to remove the chance correlations between the SNP instruments
    gsmr2_beta <<- 1     # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development  | default is 0
  }
}

get_pair = function(type = 'eas'){
  if (type == 'eas'){
    pair = read.table('~/gwas/para/pair_mtcojo.txt')[,1]
    file1 = paste0('Case_control_14_on_', pair, '_ref_ukbchn')
    file2 = paste0('Case_control_36_on_', pair, '_ref_ukbchn')
    pair = cbind(file1, file2)
    colnames(pair) = c('data1', 'data2')
    
    lcv = read.csv('~/gwas/result/lcv.csv')
    lcv = lcv[lcv$lcv_p<0.05/13, 1:2] # 6 bi-dir + lcv
    lcv = lcv[!grepl('_on_',lcv[,1]) & !grepl('_on_', lcv[,2]),]
    
    # t2d and all qtl
    ldsc = read.csv('~/gwas/result/ldsc.csv')
    ldsc = ldsc[(ldsc$data1=='Case_control_14'|ldsc$data2=='Case_control_14') & (grepl('QTL', ldsc$data1)|grepl('QTL', ldsc$data2)) & ldsc$p<0.05, 1:2]
    
    add = rbind(lcv, ldsc)
    add = rbind(add[,2:1], setNames(add[,2:1], names(add)[1:2]))
    
    pair = unique(rbind(pair, add))
  } else if (type == 'eur'){
    pair = t(combn(c('t2d', 'alt', 'urpota'), 2))
    pair = data.frame(rbind(pair, pair[,2:1]))
    pair = pair[pair[,1]=='t2d' | pair[,2]=='t2d',]
    colnames(pair) = c('data1', 'data2')
  }
  return(pair)
}


get_bloodmarker = function(type){
  markers = c('QTL_7',	'QTL_9',	'QTL_11',	'QTL_13',	'QTL_15',	'QTL_17',	'QTL_19',	'QTL_21',	'QTL_23',	'QTL_25',	'QTL_27',	'QTL_29',	'QTL_31',	'QTL_39',	'QTL_43',	'QTL_45',
              'QTL_47',	'QTL_49',	'QTL_51',	'QTL_53',	'QTL_57',	'QTL_59',	'QTL_61',	'QTL_71',	'QTL_75',	'QTL_77',	'QTL_79',	'QTL_81',	'QTL_83',	'QTL_85',	'QTL_87',	'QTL_89',
              'QTL_93',	'QTL_97',	'QTL_99',	'QTL_105',	'QTL_107',	'QTL_109',	'QTL_111',	'QTL_113',	'QTL_115',	'QTL_117',	'QTL_119',	'QTL_121')
  
  labels = c('AG',	'ALP',	'ALT',	'APTT',	'AST',	'Alb',	'BS',	'BUN',	'Baso',	'CK',	'CRP',	'Ca',	'Cl',	'Eosino',	'Fbg',	'GGT',	'HDL-C',	'Hb',	'HbA1c',	
             'Ht',	'K',	'LDH',	'LDL-C',	'Lym',	'MCH',	'MCHC',	'MCV',	'Mono',	'NAP',	'Na',	'Neutro',	'P',	'PT',	'Plt',	'RBC',	'TBil',	'TC',	'TG',	
             'TP',	'UA',	'WBC',	'ZTT',	'eGFR',	'sCr')
  
  categorys = c('Protein',	'Liver-related',	'Liver-related',	'Other biochemical',	'Liver-related',	'Protein',	'Metabolic',	'Kidney-related',	'Hematological',	'Other biochemical',	
                'Other biochemical',	'Electrolyte',	'Electrolyte',	'Hematological',	'Other biochemical',	'Liver-related',	'Metabolic',	'Hematological',	'Metabolic',	'Hematological',	'Electrolyte',
                'Other biochemical',	'Metabolic',	'Hematological',	'Hematological',	'Hematological',	'Hematological',	'Hematological',	'Protein',	'Electrolyte',	'Hematological',	'Electrolyte',	
                'Other biochemical',	'Hematological',	'Hematological',	'Liver-related',	'Metabolic',	'Metabolic',	'Protein',	'Kidney-related',	'Hematological',	'Liver-related',	'Kidney-related',	'Kidney-related')
  
  if (type == 'marker'){
    return(markers)
  } else if (type == 'label'){
    return(labels)
  } else {
    return(setNames(data.frame(cbind(markers, labels, categorys)), c('marker', 'label', 'category')))
  }
}

# Total protein (TP) | both
# Albumin/globulin ratio (AG) | liver 
# Albumin (Alb) | liver 
# Non-albumin protein (NAP) | kidney
get_lkmarker = function(type){
  markers = c('QTL_7', 	'QTL_9', 	'QTL_11', 	'QTL_15', 	'QTL_17', 	'QTL_21', 	'QTL_45', 	'QTL_83', 	'QTL_105', 	'QTL_111', 	'QTL_113', 	'QTL_117', 	'QTL_119', 	'QTL_121' )
  labels = c('AG', 	'ALP', 	'ALT', 	'AST', 	'Alb', 	'BUN', 	'GGT', 	'NAP', 	'TBil', 	'TP', 	'UA', 	'ZTT', 	'eGFR', 	'sCr')
  categorys = c('Liver', 	'Liver', 	'Liver', 	'Liver', 	'Liver', 	'Kidney', 	'Liver', 	'Kidney', 	'Liver', 	'Both', 	'Kidney', 	'Liver', 	'Kidney', 	'Kidney')
  if (type == 'marker'){
    return(markers)
  } else if (type == 'label'){
    return(labels)
  } else {
    return(setNames(data.frame(cbind(markers, labels, categorys)), c('marker', 'label', 'category')))
  }
}


get_bbj_info = function(){
  path_list = list('QTL|Case'='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/clean/',
    't2d'='/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/clean/',
    'bbj'='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/',
    'cause_ref'='/home/yanglab_data/user/zhanghy/project/mr_server/db/cause_ref/eas/')
  path_list <<- path_list
  datas_proxy <<- c('QTL_43', 'QTL_89', 'QTL_117', 'Case_control_32', 'Case_control_36', 'Case_control_38','Case_control_45',
    'Case_control_51', 'Case_control_61', 'Case_control_71', 'Case_control_75', 'Case_control_79', 'Case_control_88')
}

get_bbj_stratify_info = function(){
  path_list = list('QTL|Case'='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/',
    't2d'='/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/',
    'bbj'='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/',
    'ref'='/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/cause_ref/')
  path_list <<- path_list
  datas_proxy <<- c('')
}

get_gd_info = function(get_bim = T){
  if (get_bim==T){bim = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt.bim'); bim <<- bim%>%filter(V1==23)}
  expos = c('giant_bmi_f_2015', 'giant_hipadjbmi_f_2015', 'giant_hip_f_2015', 'giant_wcadjbmi_f_2015', 
          'giant_wc_f_2015', 'giant_whradjbmi_f_2015', 'giant_whr_f_2015', 'giant_height_f_2013', 'giant_weight_f_2013',
          'sitting_height_f')
  outcomes = c('finr5_263', 'E10_f', 'E11_f')
  pair = data.frame()
  for (outcome in outcomes){pair = rbind(pair, cbind(expos, outcome))}
  colnames(pair) = c('data1', 'data2')
  pair = rbind(pair, setNames(pair[,2:1], c('data1', 'data2')))
  pair <<- pair
  datas_proxy <<- c('E10_f', 'finr5_263') 

  path_list = list('E1|sit'='/home/yanglab_data/user/zhanghy/gwas/summary/ukb/eur/',
    'fin'='/home/yanglab_data/user/zhanghy/gwas/summary/finngen/',
    'giant'='/home/yanglab_data/user/zhanghy/gwas/summary/giant/',
    'ref'='/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/cause_ref/')
  path_list <<- path_list

  ld = read.table('/home/yanglab_data/user/zhanghy/gwas/summary/giant/result/ld/ld_out/gd_0.6.ld', header=1)
  ld <<- ld

  get_r2_para <<- function(file){
    if(!file%in%c('E10_f', 'E11_f', 'finr5_263')){
      print('file error, can not get r2 para'); q()
    }
    if (file=='E10_f'){ncase = 1357; nctrl = 228837; prev = 0.003}
    if (file=='E11_f'){ncase = 8754; nctrl = 221440; prev = 0.031}
    if (file=='finr5_263'){ncase = 5687; nctrl = 117892; prev = 0.041}
    out = list()
    out['ncase'] = ncase; out['nctrl'] = nctrl; out['prev'] = prev
    return(out)
  }
}


get_t2dcat_diff_info = function(pop, get_bim = T){
  if (pop=='eas'){
    path_list = list('QTL|Case'='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/',
      'spracklen'='/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/',
      'ref'='/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/cause_ref/',
      'bbj'='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/')
    path_list <<- path_list
    ld = list()
    ld[['m']] = read.table('/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ld/ld_out/t2d-cat_m_eas_0.6.ld', header=1)
    ld[['f']] = read.table('/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ld/ld_out/t2d-cat_f_eas_0.6.ld', header=1)
    ld <<- ld
    if (get_bim==T){bim = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/chn.bim'); bim <<- bim%>%filter(V1==23)}
    pair = matrix(c('spracklen_nature2020_t2d_f', 'Case_control_37_f', 'spracklen_nature2020_t2d_m', 'Case_control_37_m', 
                'Case_control_96_f', 'Case_control_37_f', 'Case_control_96_m', 'Case_control_37_m'), ncol=2, byrow = T)
    pair = rbind(pair, pair[, c(2, 1)])
    colnames(pair) = c('data1', 'data2')
    pair <<- pair
    datas_proxy <<- c('Case_control_37_f', 'Case_control_37_m') 
    get_r2_para <<- function(file){
      if(!grepl('Case|sprack', file)){print('file error, can not get r2 para'); q()}
      if (file=='spracklen_nature2020_t2d_m'){ncase = 28027; nctrl = 89312; prev = 0.086}; if (file=='spracklen_nature2020_t2d_f'){ncase = 27370; nctrl = 135055; prev = 0.068}
      if (file=='spracklen_nature2020_bmiadjt2d_m'){ncase = 28027; nctrl = 89312; prev = 0.086}; if (file=='spracklen_nature2020_bmiadjt2d_f'){ncase = 27370; nctrl = 135055; prev = 0.068}
      if (file=='Case_control_96_m'){ncase = 25705; nctrl = 82774; prev = 0.059}; if (file=='Case_control_96_f'){ncase = 14545; nctrl = 87841; prev = 0.07}
      if (file=='Case_control_37_m'){ncase = 11641; nctrl = 97706; prev = 0.406}; if (file=='Case_control_37_f'){ncase = 12981; nctrl = 90125; prev = 0.4233} 
      out = list()
      out['ncase'] = ncase; out['nctrl'] = nctrl; out['prev'] = prev
      return(out)
    }
  } 
  if (pop=='eur'){
    path_list = list('bbj'='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/',
      'ukb'='/home/yanglab_data/user/zhanghy/gwas/summary/ukb/eur/',
      'ref'='/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/mr_ref/cause_ref/')
    path_list <<- path_list
    ld = list()
    ld[['m']] = read.table('/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ld/ld_out/t2d-cat_m_eur_0.6.ld', header=1)
    ld[['f']] = read.table('/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ld/ld_out/t2d-cat_f_eur_0.6.ld', header=1)
    ld <<- ld
    if (get_bim==T){bim = read.table('/home/yanglab_data/user/zhanghy/gwas/ukb/bfile/brt.bim'); bim <<- bim%>%filter(V1==23)}
    path = '/home/yanglab_data/user/zhanghy/gwas/summary/ukb/eur/'
    path <<- path
    pair = matrix(c('E11_m', 'cataract_m', 'E11_f', 'cataract_f', 
      'E11_0.8_f', 'cataract-adjE11_0.8_f', 'E11_0.8_m','cataract-adjE11_0.8_m'), ncol=2, byrow = T)
    pair = rbind(pair, pair[, c(2, 1)])
    colnames(pair) = c('data1', 'data2')
    pair <<- pair
    datas_proxy <<- c('Case_control_37_f', 'Case_control_37_m') 
    get_r2_para <<- function(file){
      if(!grepl('E11|cat', file)){print('file error, can not get r2 para')}
      if (file=='E11_m'){ncase = 14059; nctrl = 181145; prev = 0.038}; if (file=='E11_f'){ncase = 8754; nctrl = 221440; prev = 0.031}
      if (file=='cataract_m'){ncase = 12137; nctrl = 183067; prev = 0.3002}; if (file=='cataract_f'){ncase = 15713; nctrl = 214481; prev = 0.4069}
      if (file=='E11_0.8_m'){ncase = 11223; nctrl = 144940; prev = 0.038}; if (file=='E11_0.8_f'){ncase = 7004; nctrl = 177151; prev = 0.031}
      if (file=='cataract-adjE11_0.8_m'){ncase = 9698; nctrl = 146465; prev = 0.3002}; if (file=='cataract-adjE11_0.8_f'){ncase = 12611; nctrl = 171544; prev = 0.4069}
      out = list()
      out['ncase'] = ncase; out['nctrl'] = nctrl; out['prev'] = prev
      return(out)
    }
  }
  if (pop=='all'){
    path <<- '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/'
    pair = matrix(c('spracklen_nature2020_t2d_f', 'Case_control_37_f', 'spracklen_nature2020_t2d_m', 'Case_control_37_m', 
    'Case_control_96_f', 'Case_control_37_f', 'Case_control_96_m', 'Case_control_37_m','E11_m', 'cataract_m', 'E11_f', 'cataract_f', 
    'E11_0.8_f', 'cataract-adjE11_0.8_f', 'E11_0.8_m','cataract-adjE11_0.8_m', 
    'cataract-adjE11_0.8_f', 'E11_0.8_f', 'cataract-adjE11_0.8_m', 'E11_0.8_m'), ncol=2, byrow = T)
    pair = rbind(pair, pair[, c(2, 1)])
    colnames(pair) = c('data1', 'data2')
    pair <<- pair
    path_list = list('QTL|Case'='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/',
      'spracklen'='/home/yanglab_data/user/zhanghy/gwas/summary/other/eas/',
      'E11|cataract'='/home/yanglab_data/user/zhanghy/gwas/summary/ukb/eur/',
      'bbj'='/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/')
    path_list <<- path_list
    datas_proxy <<- c('Case_control_37_f', 'Case_control_37_m') 
    ld_eas = list()
    ld_eas[['m']] = read.table('/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ld/ld_out/t2d-cat_m_eas_0.6.ld', header=1)
    ld_eas[['f']] = read.table('/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ld/ld_out/t2d-cat_f_eas_0.6.ld', header=1)
    ld_eur = list()
    ld_eur[['m']] = read.table('/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ld/ld_out/t2d-cat_m_eur_0.6.ld', header=1)
    ld_eur[['f']] = read.table('/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/result/ld/ld_out/t2d-cat_f_eur_0.6.ld', header=1)
    ld = list(); ld[['eas']] = ld_eas; ld[['eur']] = ld_eur
    ld <<- ld
    prev_t2d_m_eas<<-0.086; prev_t2d_f_eas<<-0.068; prev_cataract_m_eas<<-0.406; prev_cataract_f_eas<<-0.4233
    prev_t2d_m_eur<<-0.038; prev_t2d_f_eur<<-0.031; prev_cataract_m_eur<<-0.3002; prev_cataract_f_eur<<-0.4069
  }
}

get_lia_factor = function(prev_expo, prev_out){
  z_expo = dnorm(qnorm(prev_expo))
  z_out = dnorm(qnorm(prev_out))
  factor = z_expo*prev_out*(1-prev_out)/z_out/prev_expo/(1-prev_expo)
  return(factor)
}

get_harmo_res = function(df1_raw, df2_raw, pthres){
  df1 = df1_raw%>%filter(P<as.numeric(pthres))
  df2 = df2_raw%>%filter(SNP%in%df1$SNP)

  df1 = df1%>%rename(pval.exposure=P, effect_allele.exposure=A1, other_allele.exposure=A2, samplesize.exposure=N, beta.exposure=BETA, se.exposure=SE, eaf.exposure=FRQ)%>%
    mutate(id.exposure=data1, exposure=data1)
  df2 = df2%>%rename(pval.outcome=P, effect_allele.outcome=A1, other_allele.outcome=A2, samplesize.outcome=N, beta.outcome=BETA, se.outcome=SE, eaf.outcome=FRQ)%>%
    mutate(id.outcome=data2, outcome=data2)

  dat1 = harmonise_data(df1, df2)
  dat2 = dat1%>%filter(mr_keep==T)%>%
    select(CHR.x, POS.x, SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, samplesize.exposure)%>%
    rename(CHR=CHR.x, POS=POS.x, A1=effect_allele.exposure, A2=other_allele.exposure, FRQ=eaf.exposure, BETA=beta.exposure, SE=se.exposure, P=pval.exposure, N=samplesize.exposure)
  res = list(dat1, dat2)
  return(res)
}

get_mr_res = function(df, snp){
  df = df%>%filter(SNP%in%snp)
  out = mr(df)
  pleio = mr_pleiotropy_test(df)[,5:7]
  colnames(pleio) = paste0('pleio_', colnames(pleio))
  
  res = c()
  for (method in out$method){
    sub = out[out$method==method, c('nsnp', 'b', 'se', 'pval')]
    colnames(sub) = sapply(colnames(sub), function(x){paste0(method, '_', x)})
    res = c(res, sub)
  }
  
  res = c(out$id.exposure[1], out$id.outcome[1], unlist(res), pleio, pthres)
  names(res)[1:2] = c("data1", "data2")
  names(res)[length(res)] = 'pthres'
  return(res)
}

get_gsmr_res = function(df, snp, mat){
  df = df%>%filter(SNP%in%snp)
  gsmr_data = df%>%rename('a1'='effect_allele.exposure', 'a2'='other_allele.exposure', 'bzx_pval'='pval.exposure', 'bzy_pval'='pval.outcome', 'bzx'='beta.exposure',
                        'bzy'='beta.outcome', 'bzx_se'='se.exposure', 'bzy_se'='se.outcome', 'bzx_n'='samplesize.exposure', 'bzy_n'='samplesize.outcome', 'a1_freq'='eaf.exposure')
  
  ldrho = mat[rownames(mat) %in% gsmr_data$SNP, colnames(mat)%in%gsmr_data$SNP]
  snp_coeff_id = rownames(ldrho)
  
  gsmr_results = try(gsmr(gsmr_data$bzx, gsmr_data$bzx_se, gsmr_data$bzx_pval, gsmr_data$bzy, gsmr_data$bzy_se, gsmr_data$bzy_pval,
                          ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh = gwas_thresh, single_snp_heidi_thresh, 
                          multi_snp_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta))
  if (grepl('non-conformable arrays', gsmr_results[[1]])){
    gsmr_results = try(gsmr(gsmr_data$bzx, gsmr_data$bzx_se, gsmr_data$bzx_pval, gsmr_data$bzy, gsmr_data$bzy_se, gsmr_data$bzy_pval,
                            ldrho, snp_coeff_id, n_ref, F, gwas_thresh = gwas_thresh, single_snp_heidi_thresh, multi_snp_heidi_thresh, 
                            nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta))
    heidi = 0
  } else {heidi = 1}
  res = c(data1, data2, length(gsmr_results$used_index), gsmr_results$bxy, gsmr_results$bxy_se, gsmr_results$bxy_pval, heidi, pthres)
  names(res) = c('data1', 'data2', 'gsmr_nsnp', 'gsmr_b', 'gsmr_se', 'gsmr_pval', 'heidi', 'pthres')
  return(res)
}

get_cause_res = function(df, pruned, params){
  res = try(cause(X=df, variants = pruned, param_ests = params))
  if ('try-error' %in% class(res)){
    res = cause(X=df, variants = pruned, param_ests = params, force=TRUE)
    res$force = 1
  }else{res$force = 0}
  # stat
  t = summary(res)[4]$tab
  t = c(t[1,3:4], t[2, 2:4])
  names(t) = c('cause_share_eta', 'cause_share_q', 'cause_causal_gamma', 'cause_causal_eta', 'cause_causal_q')
  t = strsplit(gsub('\\(|,|\\)', '', t), ' ')
  t = as.data.frame(sapply(t, as.numeric))
  t = sapply(colnames(t), function(x){b = t[1, x];se=(t[3, x]-t[2, x])/2/1.96;p=pnorm(abs(b/se), lower.tail=FALSE)*2;return(c(b, se, p))})
  newname = as.vector(sapply(colnames(t), function(x){paste0(x, c('_b', '_se', '_pval'))}))
  stat = as.vector(t); names(stat) = newname
  stat[['cause_nsnp']] = length(res$data$snp)
  stat[['cause_elpd_pval']] = summary(res)$p
  res$stat = stat
  return(res)
}

get_snp_pleio = function(data2, data1, path, type, ld){
  pthres_re = ifelse(data2 %in% datas_proxy, '1e-5', '5e-8')
  sex = ifelse(grepl('_m', data2), 'm', 'f')
  if (type=='mr'){
    file_clump_re = paste0(path, 'result/clump/clump_out/', pthres_re, '/', data2, '&', data1, '.clumped')
    snp_re = read.table(file_clump_re, header=1)[,3]
  }
  if (type=='cause'){
    file_prune_re = paste0(path, 'result/cause/for_cause/prune/t2d_cat/', data2, '&', data1, '.rdata')
    load(file_prune_re)
    snp_re = res$pruned
  }
  snp_pleio = unique(c(snp_re, unlist(ld%>%filter(SNP_A%in%snp_re|SNP_B%in%snp_re)%>%select(SNP_A, SNP_B))))
  return(snp_pleio)
}

get_snp_pleio_gd = function(data2, path, type){
  pthres_re = ifelse(data2 %in% datas_proxy, '1e-5', '5e-8')
  if (type=='mr'){
    file_clump_re = paste0(path, 'result/clump/clump_out/', pthres_re, '/', data2, '&', data1, '.clumped')
    snp_re = read.table(file_clump_re, header=1)[,3]
  }
  if (type=='cause'){
    file_prune_re = paste0(path, 'result/cause/for_cause/prune/', data2, '&', data1, '.rdata')
    load(file_prune_re)
    snp_re = res$pruned
  }
  snp_pleio = unique(c(snp_re, unlist(ld%>%filter(SNP_A%in%snp_re|SNP_B%in%snp_re)%>%select(SNP_A, SNP_B))))
  return(snp_pleio)
}

get_mvharmo_res = function(df1_raw, df2_raw, df3_raw, data1, data2, data3, pthres){
  snp1 = df1_raw%>%filter(P<as.numeric(pthres))%>%pull(SNP)
  snp2 = df2_raw%>%filter(P<as.numeric(pthres))%>%pull(SNP)
  snp_keep = c(snp1, snp2)
  df1 = df1_raw%>%filter(SNP%in%snp_keep)
  df2 = df2_raw%>%filter(SNP%in%snp_keep)
  df3 = df3_raw%>%filter(SNP%in%snp_keep)

  df1 = df1%>%rename(pval.exposure=P, effect_allele.exposure=A1, other_allele.exposure=A2, samplesize.exposure=N, beta.exposure=BETA, se.exposure=SE, eaf.exposure=FRQ)%>%
    mutate(id.exposure=data1, exposure=data1)
  df2 = df2%>%rename(pval.exposure=P, effect_allele.exposure=A1, other_allele.exposure=A2, samplesize.exposure=N, beta.exposure=BETA, se.exposure=SE, eaf.exposure=FRQ)%>%
    mutate(id.exposure=data2, exposure=data2)
  df3 = df3%>%rename(pval.outcome=P, effect_allele.outcome=A1, other_allele.outcome=A2, samplesize.outcome=N, beta.outcome=BETA, se.outcome=SE, eaf.outcome=FRQ)%>%
    mutate(id.outcome=data3, outcome=data3)

  mvdat = mv_harmonise_data(rbind(df1, df2), df3)
  snp = rownames(mvdat$exposure_beta)
  res = list(mvdat, snp)
  return(res)
}
