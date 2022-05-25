library(stringr)
remove_brackets_split = function(x){
  if (length(line) >0) {
    x = str_squish(x) # multiple spaces to single
    x = gsub(')', '', gsub('(', '', x, fixed=T), fixed=T)
    x = strsplit(x, ' ')[[1]]
  }
  return(x)
}

collect_ldsc = function(path, files){
    out = c()
    file_error = c()
    file_warn = c()
    col_name = c('data1', 'data2', 'constrain', 'liability', 'nsnp', 'h2_1', 'h2_1_se', 'int_1', 'int_1_se', 'h2_2', 'h2_2_se', 'int_2', 'int_2_se',
                'rg', 'rg_se', 'rg_p', 'int_bi', 'int_bi_se')

    for (file in files){
    warn = F
    print(file)
    if (sum(grepl("WARN", readLines(paste0(path, file))))){
        file_warn = c(file_warn, file); next()
    }
    data1 = strsplit(file, '&')[[1]][1]
    data2 = gsub('_lia', '', gsub('_constrain', '', gsub('.log', '', strsplit(file, '&')[[1]][2])))

    con = file(paste0(path, file))
    if (length(readLines(con)) < 62) {
        file_error = c(file_error, file)
        close(con)
        next()
    }
    close(con)
    con = file(paste0(path, file),open="r")
    n = 1
    while ( TRUE ) {
        line = readLines(con, n = 1)
        # cat(n,line,"\n") # first print and see, then write code
        line = remove_brackets_split(line)
        if ( n > 64 ) {
        break
        }
        if (grepl('_lia', file)){
        liability = 1; shift = 2} else {liability = 0; shift = 0}
        if (grepl('_constrain', file)){
        constrain = 1
        int_1 = 1
        int_1_se = 0
        int_2 = 1
        int_2_se = 0
        int_bi = 1
        int_bi_se = 0
        if (n == 30 + shift){
            nsnp = line[1]
        }
        if (n == 34 + shift){
            h2_1 = line[5]
            h2_1_se = line[6]
        }
        if (n == 41 + shift){
            h2_2 = line[5]
            h2_2_se = line[6]
        }
        if (n == 60 + shift){
            col = line
        }
        if (n == 61 + shift){
            value = line
            rg = value[3]; rg_se = value[4]; rg_p = value[6]
        }
        } else {
        constrain = 0
        if (n == 29 + shift){
            nsnp = line[1]
        }
        if (n == 33 + shift){
            h2_1 = line[5]; h2_1_se = line[6]
        }
        if (n == 36 + shift){
            int_1 = line[2]; int_1_se = line[3]
        }
        if (n == 41 + shift){
            h2_2 = line[5]; h2_2_se = line[6]
        }
        if (n == 44 + shift){
            int_2 = line[2]; int_2_se = line[3]
        }
        if (n == 61 + shift){
            col = line
        }
        if (n == 62 + shift){
            value = line
            rg = value[3]; rg_se = value[4]; rg_p = value[6]; int_bi = value[11]; int_bi_se = value[12]
        }
        }
        n = n+1
    }
    close(con)
    if (warn != T){
        out = c(out, data1, data2, constrain, liability, nsnp, h2_1, h2_1_se, int_1, int_1_se, h2_2, h2_2_se, int_2, int_2_se, 
                rg, rg_se, rg_p, int_bi, int_bi_se)
        rm(list = col_name)
    }
    }

    res = data.frame(matrix(out, ncol = 18, byrow=T))
    colnames(res) = col_name
    res = res[res$data1!= '' & res$data2 != '',]

    cat('n error file: ', length(file_error), '\n')
    cat('n warn file: ', length(file_warn), '\n')

    res[,c('h2_1', 'h2_2', 'h2_1_se', 'h2_2_se')] = sapply(res[,c('h2_1', 'h2_2', 'h2_1_se', 'h2_2_se')], as.numeric)
    res$h2_1_p = 2*pnorm(abs(res$h2_1/res$h2_1_se), lower.tail=F)
    res$h2_2_p = 2*pnorm(abs(res$h2_2/res$h2_2_se), lower.tail=F)

    return(res)
}


