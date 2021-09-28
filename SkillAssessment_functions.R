######################################################################
## SKill assessment comparison test functions 
######################################################################

#######################################################
## Mean test of each PC
#######################################################
##d2 is the d2-th basis functions (or PCs) being considered;

test.sn.mean.each<-function(d2)  {
  ProX.SN<-as.matrix(scores.un.x[,d2])
  ProY.SN<-as.matrix(scores.un.y[,d2])
  
  ###SN test
  
  N0<-2*N  ##When the sample sizes are unequal, replace N0 by N1+N2.
  index.X<-floor((1:N0)*(N/N0))
  index.Y<-floor((1:N0)*(N/N0))
  aX<-sum(index.X==0)
  aY<-sum(index.Y==0)
  VX<-(diag(c(1:N)^(-1))%*%apply(ProX.SN,2,cumsum))[index.X,]
  VY<-(diag(c(1:N)^(-1))%*%apply(ProY.SN,2,cumsum))[index.Y,]
  V<-rbind(matrix(0,aX,1),as.matrix(VX))-rbind(matrix(0,aY,1),as.matrix(VY))
  
  de.V<-t(V)-V[N0,]
  D.self<-de.V%*%diag(c(1:N0)^2)%*%t(de.V)
  G.self<-N0^3*t(V[N0,])%*%solve(D.self, tol = 1e-30)%*%t(t(V[N0,]))
  return(G.self)
}

#######################################################
## Mean test of cumulative PCs
#######################################################
##d1 is the number of basis functions (or PCs) being considered; it can be chosen using the cumulative variation rule.
test.sn.mean<-function(d1)  {
  ProX.SN<-as.matrix(scores.un.x[,1:d1])
  ProY.SN<-as.matrix(scores.un.y[,1:d1])
  
  ###SN test
  
  N0<-2*N  ##When the sample sizes are unequal, replace N0 by N1+N2.
  index.X<-floor((1:N0)*(N/N0))
  index.Y<-floor((1:N0)*(N/N0))
  aX<-sum(index.X==0)
  aY<-sum(index.Y==0)
  VX<-(diag(c(1:N)^(-1))%*%apply(ProX.SN,2,cumsum))[index.X,]
  VY<-(diag(c(1:N)^(-1))%*%apply(ProY.SN,2,cumsum))[index.Y,]
  V<-rbind(matrix(0,aX,d1),as.matrix(VX))-rbind(matrix(0,aY,d1),as.matrix(VY))
  
  de.V<-t(V)-V[N0,]
  D.self<-de.V%*%diag(c(1:N0)^2)%*%t(de.V)
  G.self<-N0^3*t(V[N0,])%*%solve(D.self, tol = 1e-30)%*%t(t(V[N0,]))
  return(G.self)
}

#######################################################
## Covariance test of each PC
#######################################################
test.sn.cov.each<-function(d)  ##d is the number of basis functions (or PCs) being considered; it can be chosen using the cumulative variation rule.
{
  p<-1
  # Covariance
  N0<-2*N
  new.covxy  = matrix(nrow = 1000, ncol = p)
  
  for ( j in 1:1000){
    
    new.tmpxy = c()
    tmp.x = (demean.x[,1:j] %*% t(demean.x[,1:j]))/j	# Cx.1,..., Cx.1000
    tmp.y = (demean.y[,1:j] %*% t(demean.y[,1:j]))/j	# Cy.1,..., Cy.1000
    tmp   = tmp.x - tmp.y
    
    # Multiply the eigenfunction (dimension reduction)
    tmp = tmp[1:24, -c(1:24)]
    covar.xy  = t(u.svd[,d]) %*% tmp %*% v.svd[,d]		# force the not interested pairs to be 0
    new.covxy[j,] = covar.xy		
    
  }
  #######Self-normalization
  N0<-2*N  ##When the sample sizes are unequal, replace N0 by N1+N2.
  index.X<-floor((1:N0)*(N/N0))
  aX<-sum(index.X==0)
  
  re.xy<- new.covxy[index.X,]
  
  hat.a<-rbind(matrix(0,aX,p),as.matrix(re.xy))
  de.hat.a<-t(hat.a[,])-hat.a[N0,]
  D.self<-de.hat.a%*%diag(c(1:N0)^2)%*%t(de.hat.a)
  G.self<-N0^3*t(hat.a[N0,])%*%solve(D.self, tol = 1e-30)%*%t(t(hat.a[N0,]))
  return(G.self)
  
}
#######################################################
## Covariance test of cumulative PCs
#######################################################
test.sn.cov<-function(d1)  ##d is the number of basis functions (or PCs) being considered; it can be chosen using the cumulative variation rule.
{
  p<-d1*(d1+1)/2
  
  # Covariance
  N0<-2*N
  new.covxy  = matrix(nrow=N, ncol=p)
  
  for ( j in 1:1000){
    new.tmpxy = c()
    tmp.x = (demean.x[,1:j] %*% t(demean.x[,1:j]))/j	# Cx.1,..., Cx.1000
    tmp.y = (demean.y[,1:j] %*% t(demean.y[,1:j]))/j	# Cy.1,..., Cy.1000
    tmp   = tmp.x - tmp.y
    
    # Multiply the eigenfunction (dimension reduction)
    tmp = tmp[1:24, -c(1:24)]
    covar.xy  = t(u.svd[,1:d1]) %*% tmp %*% v.svd[,1:d1]		# force the not interested pairs to be 0
    new.covxy[j,] = covar.xy[upper.tri(covar.xy, diag=TRUE)]	
  }
  #######Self-normalization
  N0<-2*N  ##When the sample sizes are unequal, replace N0 by N1+N2.
  index.X<-floor((1:N0)*(N/N0))
  aX<-sum(index.X==0)
  re.xy<- new.covxy[index.X,]
  
  hat.a<-rbind(matrix(0,aX,p),as.matrix(re.xy))
  de.hat.a<-t(hat.a[,])-hat.a[N0,]
  D.self<-de.hat.a%*%diag(c(1:N0)^2)%*%t(de.hat.a)
  G.self<-N0^3*t(hat.a[N0,])%*%solve(D.self, tol = 1e-30)%*%t(t(hat.a[N0,]))
  return(G.self)
}
