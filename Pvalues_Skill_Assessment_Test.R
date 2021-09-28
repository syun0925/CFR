library(fields)
library(fda)

setwd("/Users/sooinyun/Dropbox/Research/CFR/Revision/Comparison Random Field")

#####################################################################################
## 1. DATA load
#####################################################################################

bcc <- read.table('target_bcc.dat', sep='') 
ccsm <- read.table('target_ccsm.dat', sep='')
giss <- read.table('target_giss.dat', sep='')
ipsl <- read.table('target_ipsl.dat', sep='')
mpi <- read.table('target_mpi.dat', sep='')

# reconstruction
recon4 <- c('cca', 'ridge','ttls', 'ttlh')

cca.bcc <- read.table('recon_cca_bcc80.dat', sep='')
ridge.bcc <- read.table('recon_ridge_bcc80.dat', sep='') 
ttls.bcc <- read.table('recon_ttls_bcc80.dat', sep='') 
ttlh.bcc <- read.table('recon_ttlshyb_bcc80.dat', sep='') 

cca.ccsm <- read.table('recon_cca_ccsm80.dat', sep='')
ridge.ccsm <- read.table('recon_ridge_ccsm80.dat', sep='') 
ttls.ccsm <- read.table('recon_ttls_ccsm80.dat', sep='') 
ttlh.ccsm <- read.table('recon_ttlshyb_ccsm80.dat', sep='') 

cca.giss <- read.table('recon_cca_giss80.dat', sep='')
ridge.giss <- read.table('recon_ridge_giss80.dat', sep='') 
ttls.giss <- read.table('recon_ttls_giss80.dat', sep='') 
ttlh.giss <- read.table('recon_ttlshyb_giss80.dat', sep='') 

cca.ipsl <- read.table('recon_cca_ipsl80.dat', sep='')
ridge.ipsl <- read.table('recon_ridge_ipsl80.dat', sep='') 
ttls.ipsl <- read.table('recon_ttls_ipsl80.dat', sep='') 
ttlh.ipsl <- read.table('recon_ttlshyb_ipsl80.dat', sep='') 

cca.mpi <- read.table('recon_cca_mpi80.dat', sep='')
ridge.mpi <- read.table('recon_ridge_mpi80.dat', sep='') 
ttls.mpi <- read.table('recon_ttls_mpi80.dat', sep='') 
ttlh.mpi <- read.table('recon_ttlshyb_mpi80.dat', sep='') 

gridloc <- read.table('locs_pnas.dat', sep='') 

######################################################
## 2. Smooth data: B-Splines
######################################################

sp1 <- create.bspline.basis(rangeval=c(-177.5,177.5), nbasis=20, norder=4) ## What is coefficient?
sp2 <- create.bspline.basis(rangeval=c(-77.5,82.5), nbasis=9, norder=4)

eval.sp1 <- eval.basis(gridloc[,1], sp1)
eval.sp2 <- eval.basis(gridloc[,2], sp2)

## creat design matrix for least square esitmates of coeff ##

eval.sp <- matrix(NA, 1732, ncol(eval.sp1)*ncol(eval.sp2))
for (i in 1:1732){
  eval.sp[i,] <- kronecker(eval.sp1[i,], eval.sp2[i,])   
}

## Smoothing Estimate
Smooth = function(data) { 
  
  data = eval.sp%*%solve(t(eval.sp)%*%eval.sp, t(eval.sp)%*%t(data)) 
  return(data)
}

## Smooth all data
cca.bcc <- t(Smooth(cca.bcc))
ridge.bcc <- t(Smooth(ridge.bcc))
ttls.bcc <- t(Smooth(ttls.bcc))
ttlh.bcc <- t(Smooth(ttlh.bcc))

cca.ccsm <- t(Smooth(cca.ccsm))
ridge.ccsm <- t(Smooth(ridge.ccsm))
ttls.ccsm <- t(Smooth(ttls.ccsm))
ttlh.ccsm <- t(Smooth(ttlh.ccsm))

cca.giss <- t(Smooth(cca.giss))
ridge.giss <- t(Smooth(ridge.giss))
ttls.giss <- t(Smooth(ttls.giss))
ttlh.giss <- t(Smooth(ttlh.giss))

cca.ipsl <-  t(Smooth(cca.ipsl))
ridge.ipsl <-  t(Smooth(ridge.ipsl))
ttls.ipsl <-  t(Smooth(ttls.ipsl))
ttlh.ipsl <-  t(Smooth(ttlh.ipsl))

cca.mpi <- t(Smooth(cca.mpi))
ridge.mpi <- t(Smooth(ridge.mpi))
ttls.mpi <- t(Smooth(ttls.mpi))
ttlh.mpi <-  t(Smooth(ttlh.mpi))

bcc  <- t(Smooth(bcc))  #1732*1146
ccsm <- t(Smooth(ccsm))
giss <- t(Smooth(giss))
ipsl <- t(Smooth(ipsl))
mpi  <- t(Smooth(mpi))

models <- c('bcc', 'ccsm', 'giss', 'ipsl', 'mpi')
cca = paste0("cca.",models)
ridge = paste0("ridge.",models)
ttls =paste0("ttls.",models)
ttlh = paste0("ttlh.",models)
model = c('bcc', 'ccsm', 'giss', 'ipsl', 'mpi', cca, ridge,ttls,ttlh)

################################################################################################
## 3. Skill Assessment Reconstruction vs. All (850-1850 vs. 1850-1995) of MEAN
################################################################################################
#######################################################
## 3.1. TEST Mean based on each PC
#######################################################
################################################################################################
## Mean comparison test of Climate simulation models vs. Climate Field Reconstructions
## Test Mean of each PC 
################################################################################################

N=1000  # number of time points
D=1732  # number of spatial locations
d=5     # number of PCs to test
m=0

model.ratio <- matrix(NA, 20, d)
model.mean.ts = model.cov.ts  = matrix(NA, 20, d)  ## d is the user chosen number of basis function
s = Sys.time()
for(i in 1:length(models)){ ## 5 Climate simulation models
  for (j in 1:4){  ## 4 different CFRs
    
    xx <- get(model[i])[1:N,]
    yy <- get(model[5*j+i])[1:N,]
    
    upward <- 0.5*(rowMeans(xx) + rowMeans(yy))
    xx <- t(xx-upward)
    yy <- t(yy-upward)
    
    demean.x <- xx - rowMeans(xx) ## remove time mean at each spatial location
    demean.y <- yy - rowMeans(yy)
    
    demean.x <- demean.x*sqrt(abs(cos(gridloc[,2]*pi/180)))
    demean.y <- demean.y*sqrt(abs(cos(gridloc[,2]*pi/180)))
    
    covar   = (demean.x %*% t(demean.x))/N
    
    re.all= eigen(covar)
    
    ve.all <- (re.all$vector[,1:d])*sqrt(D)
    values.all <- re.all$values[1:d]/D
    m = m + 1
    
    ve.all <- (re.all$vector[,1:d])*sqrt(D)
    values.all <- re.all$values[1:d]/D
    
    ####### Mean ##########
    xx <- xx*sqrt(abs(cos(gridloc[,2]*pi/180)))
    yy <- yy*sqrt(abs(cos(gridloc[,2]*pi/180)))
    scores.un.x <- t(xx) %*% ve.all/D  ## The test statistics is codmputed using the original data so it keeps the original spatial varying mean
    scores.un.y <- t(yy) %*% ve.all/D
    
    for (k in 1:d){
      test <- try(test.sn.mean.each(k))
      if (inherits(test, "error")) next
      if (!inherits(test, "error")) model.mean.ts[m,k] <- test	
    }
    
    print(c(i,j))
  }
}
e = Sys.time(); e-s
##############################################
## Derive the p-values of the mean comparison 
##############################################
setwd('/Users/sooinyun/Dropbox/Research/CFR/Revision/rdata')
load('crit.RData')    ## Load empirical distribution of W_q from Zhang and Shao (2015)

### model.cov.ts and model.ts are charactor matrix, transform them into numeric
model1.mean.ts <- matrix(NA,20, d)
for(i in 1:20){
  for (j in 1:length(which(!is.na(model.mean.ts[i,])))){
    model1.mean.ts[i,j] <- as.numeric(model.mean.ts[i,j])
  }
}
model.mean.ts <- model1.mean.ts
model.mean.pv <- model.mean.ts   ## mean surface
for (i in 1:20){
  for (j in 1:length(which(!is.na(model.mean.ts[i,])))){
    if (!is.na(model.mean.ts[i,j])) model.mean.pv[i,j] <- mean(critiq[j,]>model.mean.ts[i,j])		
  }	
}

## adjusted p-value using FDR for mean
model.mean.pv1 <- c(model.mean.pv)
tmp <- !is.na(model.mean.pv1)
tmp1 <- p.adjust(model.mean.pv1[tmp], method='fdr')
model.mean.pv1[tmp] <- tmp1
model.mean.pval <- matrix(model.mean.pv1, nrow=20)

#######################################################
## 3.2. TEST Mean based on cumulative PCs
#######################################################
################################################################################################
## Mean comparison test of Climate simulation models vs. Climate Field Reconstructions
## Test Mean of cumulative PCs
################################################################################################

d=50    ## User choice: number of PCs to test
m=0

model.ratio <- matrix(NA, 20, d); model.thresh <- matrix(NA, 20, 3)
model.mean.ts =  matrix(NA, 20, d)  ## 58 is the maxim numbe of 80% variation

s = Sys.time()
for(i in 1:length(models)){
  for (j in 1:4){
    
    xx <- get(model[i])[1:N,]
    yy <- get(model[5*j+i])[1:N,]
    
    upward <- 0.5*(rowMeans(xx) + rowMeans(yy))
    xx <- t(xx-upward)
    yy <- t(yy-upward)
    
    demean.x <- xx - rowMeans(xx) ## remove time mean at each spatial location
    demean.y <- yy - rowMeans(yy)
    
    demean.x <- demean.x*sqrt(abs(cos(gridloc[,2]*pi/180)))
    demean.y <- demean.y*sqrt(abs(cos(gridloc[,2]*pi/180)))
    
    covar   = (demean.x %*% t(demean.x))/(ncol(demean.x))
    
    re.all= eigen(covar)
    
    ve.all <- (re.all$vector[,1:d])*sqrt(D)
    values.all <- re.all$values[1:d]/D
    ratio.all <- re.all$values/sum(re.all$values)
    m = m + 1
    
    model.ratio[m,] <- cumsum(ratio.all[1:d])
    tmp <- model.ratio[m,]
    model.thresh[m,] <- d - c(sum(tmp>0.85), sum(tmp>0.9), sum(tmp>0.95)) + 1
    
    d1 <- model.thresh[m,1]    ### Subject to change depending on the cumulative percentage
    
    ve.all <- (re.all$vector[,1:d1])*sqrt(D)
    values.all <- re.all$values[1:d1]/D
    
    ####### Mean ##########
    xx <- xx*sqrt(abs(cos(gridloc[,2]*pi/180)))
    yy <- yy*sqrt(abs(cos(gridloc[,2]*pi/180)))
    scores.un.x <- t(xx) %*% ve.all/D  ## The test statistics is codmputed using the original data so it keeps the original spatial varying mean
    scores.un.y <- t(yy) %*% ve.all/D
    
    for (k in 1:d1){
      test <- try(test.sn.mean(k))
      if (inherits(test, "error")) next
      if (!inherits(test, "error")) model.mean.ts[m,k] <- test	
    }
    print(c(i,j))
  }
}
e = Sys.time(); e-s

##############################################
## Derive the p-values of the mean comparison 
##############################################

### model.cov.ts and model.ts are charactor matrix, transform them into numeric
model1.mean.ts <- model1.cov.ts <- matrix(NA,20, 50)
for(i in 1:20){
  for (j in 1:length(which(!is.na(model.mean.ts[i,])))){
    model1.mean.ts[i,j] <- as.numeric(model.mean.ts[i,j])
  }
  model.mean.ts <- model1.mean.ts
  model.mean.pv <- model.mean.ts   ## mean surface
  for (i in 1:20){
    for (j in 1:length(which(!is.na(model.mean.ts[i,])))){
      if (!is.na(model.mean.ts[i,j])) model.mean.pv[i,j] <- mean(critiq[j,]>model.mean.ts[i,j])		
    }	
  }
  
  ## adjusted p-value using FDR for mean
  model.mean.pv1 <- c(model.mean.pv)
  tmp <- !is.na(model.mean.pv1)
  tmp1 <- p.adjust(model.mean.pv1[tmp], method='fdr')
  model.mean.pv1[tmp] <- tmp1
  model.cummean.pval <- matrix(model.mean.pv1, nrow=20)

  ################################################################################################
  ## 4. Skill Assessment Reconstruction vs. All (850-1850 vs. 1850-1995) of COVARIANCE
  ################################################################################################

  ## Select ENSO dominant region
  # Nino3 Region
  id0 = which(gridloc[,1]<=-90 & gridloc[,1]>=-150 & abs(gridloc[,2]) <=5)
  near = which((gridloc[,1]<= -35 | gridloc[,1]>=150) & abs(gridloc[,2])<=20)
  
  # Correlation test function
  Corr.test = function(y) { 
    tmp = apply(y[,id0],1, mean)
    res = apply(y,2, function(x) cor.test(x, tmp)$p.value)	
    res[id0] = NA
    return(res)
  }
  
  bcc.cor.p = Corr.test(bcc[1:N,])
  ccsm.cor.p = Corr.test(ccsm[1:N,])
  giss.cor.p = Corr.test(giss[1:N,])
  ipsl.cor.p = Corr.test(ipsl[1:N,])
  mpi.cor.p = Corr.test(mpi[1:N,])

  ## Choose region with significant correlation
  bcc.sig = (which(bcc.cor.p <0.1))
  ccsm.sig = (which(ccsm.cor.p <0.1))
  giss.sig = (which(giss.cor.p <0.1))
  ipsl.sig = (which(ipsl.cor.p <0.1))
  mpi.sig = (which(mpi.cor.p <0.1))
  
  tel.sig = intersect(intersect(intersect(intersect(bcc.sig, ccsm.sig), giss.sig), ipsl.sig),mpi.sig)
  tel.sig = tel.sig[!(tel.sig %in% near)]    ### ENSO dominant region
  
  ## (1) INPUT 0 in covariance matrix
  # We will exclude the region that show insignificant correlation by replacing the values to zero
  Cov.Zero = function (dat) { 
    dat[1:length(id0), 1:length(id0)]  = 0
    dat[25:nrow(dat), 25:nrow(dat)]=0
    dat
  }
 
  ###################################################################
  ## 4.1. TEST Covariance based on each PC (parallel computation)
  ###################################################################
  # parameters
  id.tel = c(id0, tel.sig)   ### Select ENSO dominant region
  D=length(id.tel); d=5      ### Number of PCs to test
    
  # Choose models to compare
  library(doParallel)
  registerDoParallel(cores=12)
  detectCores()
  i=5   ## Test i-th climate simulation model vs. its CFRs
  DatSet = lapply(model[c(i,5*(1:4)+i)], function(x) get(x))
  
  s = Sys.time()
  ptime <- system.time({
    r <- foreach(j=1:4,.combine='rbind',.export=c('test.sn.cov.each')
                 , .packages = c("fields", "fda")) %dopar% {
                   
                   xx <- DatSet[[1]][1:1000,id.tel]
                   yy <- DatSet[[j+1]][1:1000,id.tel]
                   
                   upward <- 0.5*(rowMeans(xx) + rowMeans(yy))
                   xx <- t(xx-upward)
                   yy <- t(yy-upward)
                   
                   demean.x <- xx - rowMeans(xx) ## remove time mean at each spatial location
                   demean.y <- yy - rowMeans(yy)
                   
                   demean.x <- demean.x*sqrt(abs(cos(gridloc[id.tel,2]*pi/180)))
                   demean.y <- demean.y*sqrt(abs(cos(gridloc[id.tel,2]*pi/180)))
                   
                   covar   = (demean.x %*% t(demean.x))/N
                   covar   = covar[1:24,-c(1:24)]
                   
                   re.svd = svd(covar)
                   u.svd = re.svd$u; v.svd = re.svd$v; 
                   model.cov.ts = c()
                   
                   for (l in 1:d){
                     
                     test <- try(test.sn.cov.each(l))
                     if (inherits(test, "error")) next
                     if (!inherits(test, "error")) model.cov.ts[l] <- test ### Test stat of the covariance test
                   }
                   model.cov.ts
                 }
  })
  ptime
  
  #########################################################
  ## Derive the p-values of the covariance comparison 
  #########################################################
  model1.cov.ts <- matrix(NA,4, 5)
  for ( i in 1:4){
    for (j in 1:5){
      model1.cov.ts[i,j] <- as.numeric(r[i,j])
    }
  }
  
  model.cov.ts <- model1.cov.ts
  
  model.cov.pv <- model.cov.ts   ## covaraince structure
  for (i in 1:4){
    for (j in 1:5){		
      if (!is.na(model.cov.ts[i,j])) model.cov.pval[i,j] <- mean(critiq[1,]>model.cov.ts[i,j])		### P-values of testing each PC
    }	
  }
  
  #######################################################
  ## 4.2. TEST Covariance based on cumulative PCs
  #######################################################
  # parameters
  D=length(id.tel)
  d=24
  s = Sys.time()
  # Choose models to compare
  model.cov.ts  = matrix(NA, 20, d)  ### d is the user chosen number 
  m=0; model.thresh = matrix(nrow=20, ncol=3); 
  model.ratio = matrix(nrow=20, ncol=d)
  
  for(i in 1:length(models)){
    for (j in 1:4){
      
      xx <- get(model[i])[1:1000,id.tel]
      yy <- get(model[5*j+i])[1:1000,id.tel]
      
      upward <- 0.5*(rowMeans(xx) + rowMeans(yy))
      xx <- t(xx-upward)
      yy <- t(yy-upward)
      
      demean.x <- xx - rowMeans(xx) ## remove time mean at each spatial location
      demean.y <- yy - rowMeans(yy)
      
      demean.x <- demean.x*sqrt(abs(cos(gridloc[id.tel,2]*pi/180)))
      demean.y <- demean.y*sqrt(abs(cos(gridloc[id.tel,2]*pi/180)))
      
      covar   = (demean.x %*% t(demean.x))/N
      covar   = covar[1:24,-c(1:24)]
      
      re.svd = svd(covar)
      u.svd = re.svd$u; v.svd = re.svd$v
      m = m + 1
      
      model.ratio[m,] <- cumsum(re.svd$d/sum(re.svd$d))
      tmp <- model.ratio[m,]
      model.thresh[m,] <- d - c(sum(tmp>0.85), sum(tmp>0.9), sum(tmp>0.95)) + 1
      
      d1 <- model.thresh[m,3]    ### Subject to change depending on the cumulative percentage
      
      model.cov.tmp = c()
      for (l in 1:d1){
        
        test <- try(test.sn.cov(l))
        if (inherits(test, "error")) next
        if (!inherits(test, "error")) model.cov.ts[m,l] <- test 
      }
      print(c(i,j))
    }
  }
  e = Sys.time(); e-s
  
  ### Save the test statistics of five climate model comparisons: 
  ### Store test statistics of covariance test of ith climate model at: 'CovComp_Cumi.Rdata'
  #setwd('C:/Users/ysoo8/Dropbox/Research/Climate Model_data/Re 2019/Revision/reanalysis')
  # save(model.cov.ts, file=paste('CovComp_Cum', i,'.Rdata', sep=''))
  model.cov.ts.bind = c()   # store test statistics of all five models
  model.ratio.bind = c()    # store the ratio of the variability of cumulative PCs
  for ( i in 1:5) { 
    datname = paste('CovComp_Cum', i,'.Rdata', sep='')
    load(datname)
    model.cov.ts.bind = rbind(model.cov.ts.bind, model.cov.ts)
    model.ratio.bind = rbind(model.ratio.bind, model.ratio)
  }
  
  #########################################################
  ## Derive the p-values of the covariance comparison 
  #########################################################
  model1.cov.ts <- matrix(NA,20,5)
  for(i in 1:20){
    for (j in 1:5){
      model1.cov.ts[i,j] <- as.numeric(model.cov.ts.bind[i,j])
    }
  }

  model.cov.pv <- model1.cov.ts   ## covaraince structure
  for (i in 1:20){
    for (j in 1:5){		
      if (!is.na(model.cov.ts.bind[i,j])) model.cov.pv[i,j] <- mean(critiq[1,]>model.cov.ts.bind[i,j])		
    }	
  }
  
  ## adjusted p-value using FDR for covariance
  model.cov.pv1 <- c(model.cov.pv)
  tmp <- !is.na(model.cov.pv1)
  tmp1 <- p.adjust(model.cov.pv1[tmp], method='fdr')
  model.cov.pv1[tmp] <- tmp1
  model.covcum.pval <- matrix(model.cov.pv1, nrow=20)   ### P-values of testing cumulative PCs
  
  