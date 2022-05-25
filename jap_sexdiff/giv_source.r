library(MASS)
library(foreach)
library(doParallel)
library(dplyr)

path_ukb = '/home/yanglab_data/user/zhanghy/gwas/ukb/'
path_bbj = '/home/yanglab_data/user/zhanghy/gwas/summary/bbj/stratify/'

# parallel
no_cores <- detectCores() - 1
cat('\nno cores available: ',detectCores(), '\n')
cat('no cores used: ',no_cores, '\n')
cl <- makeCluster(no_cores, type="FORK") # FORK does not work on windows
registerDoParallel(cl) # Register paralell pool so GWAS utilise multiple cores


sample_3set = function(df, var){
    sample = 1:nrow(df)
    for (i in unique(df[,var])){
        t = sample[df[,var]==i]
        t12 = sample(t, 0.8*length(t))
        t1 = sample(t12, 0.5*length(t12))
        t2 = t12[!t12%in%t1]
        t3 = t[!t%in%t12]
        sample[sample%in%t1] = 'a'; sample[sample%in%t2] = 'b'; sample[sample%in%t3] = 'c'
    }
    return(sample)
}

allocate_mat = function(df){
    break1 <<- table(df$group)[1]; break2 <<- table(df$group)[1]+table(df$group)[2]
    rep = nrow(df)-break2
    ones1 <<- matrix(1L,nrow=break1,ncol=1)  # First GWAS subsample
    ones2 <<- matrix(1L,nrow=break2-break1,ncol=1)  # Second GWAS subsample
    onesr <<- matrix(1L,nrow=rep,ncol=1)  # Replication sample
    GIV_U <<- matrix(0L,nrow=rep,ncol=3)
    GIV_C <<- matrix(0L,nrow=rep,ncol=3)
    OLS_T <<- matrix(0L,nrow=rep,ncol=2)
    OLS_T_Sy <<-  matrix(0L,nrow=rep,ncol=3)
    OLS_T_Sy_cond <<- matrix(0L,nrow=rep,ncol=3)
    MR <<- matrix(0L,nrow=rep,ncol=2)
    EMR <<- matrix(0L,nrow=rep,ncol=3)
    EMR2 <<- matrix(0L,nrow=rep,ncol=3)
    cov_GIV_U <<- array(0L, c(3, 3, rep))
    cov_GIV_C <<- array(0L, c(3, 3, rep))
    cov_OLS_T <<- array(0L, c(2, 2, rep))
    cov_OLS_T_Sy <<- array(0L, c(3, 3, rep))
    cov_OLS_T_Sy_cond <<- array(0L, c(3, 3, rep))
    cov_MR <<- array(0L, c(2, 2, rep))
    cov_EMR <<- array(0L, c(3, 3, rep))
    cov_EMR2 <<- array(0L, c(3, 3, rep))
}





# Fast OLS
fols = function (y, x) {
     XtX = crossprod(x)
     Xty = crossprod(x, y)
     solve(XtX, Xty)
}

# calculate variance-covariance matrix of OLS
olscov = function(y,X,bhat,n) {
    k = dim(bhat)[1]
    e = y - X %*% bhat
    s2 = crossprod(e)/(n-k)
    return(solve(crossprod(X)) * s2[1,1]) # first element of 1 by 1 matrix
}

# calculate variance-covariance matrix of 2SLS
ivcov = function(y,X,Xpred,bhat,n) {
    e = y - X %*% bhat
    dim(X)
    dim(bhat)
    s2 = crossprod(e)/n
    return(solve(crossprod(Xpred)) * s2[1,1]) # first element of 1 by 1 matrix
}