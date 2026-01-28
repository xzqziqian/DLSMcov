
# 
# Extension of Sewell & Chen (2015) to include dynamic non-network covariate evolution. Based on the dynamic latent space model code from https://publish.illinois.edu/yuguo/software/ and https://sewell.lab.uiowa.edu/software
# 
# Reference: Sewell, D. K., & Chen, Y. (2015). Latent space models for dynamic networks. Journal of the american statistical association, 110(512), 1646-1657.


#######################
###     Startup     ###
#######################

n <- dim(Y)[1]
TT <- dim(Y)[3]
ncov = dim(C)[2]

###
###Missing Data
###
set.seed(1)
# if(MissData){
#   Yaggreg <- array(0,dim(Y[,,1]))
#   for(tt in 1:TT) Yaggreg <- Yaggreg + Y[,,tt]
#   Yaggreg[which(Yaggreg>1,arr.ind=TRUE)] <- 1
#   SPaggreg <- shortest.paths(graph=graph.adjacency(Yaggreg),mode="all")
#   SPaggreg[which(SPaggreg==Inf,arr.ind=TRUE)] <- 5
# 
#   outDeg <- inDeg <- denom <- array(0, dim = (n))
#   for(tt in 1:TT){
#     denom[c(1:n)[Missing[[tt]]]] <- denom[c(1:n)[Missing[[tt]]]] + 1
#     if(length(Missing[[tt]])==0) denom <- denom + 1
#     inDeg <- inDeg + colSums(Y[,,tt])/(n-1-length(Missing[[tt]])+as.numeric(1:n%in%Missing[[tt]]))*(n-1)
#     outDeg <- outDeg + rowSums(Y[,,tt])
#   }
#   inDeg= inDeg/TT
#   outDeg <- round(outDeg/denom)
# 
#   for(tt in 1:TT){
#     for(i in Missing[[tt]]){
#       Probs <- inDeg[-i]/SPaggreg[i,-i]
#       ind <- sample(size=outDeg[i],x=c(1:n)[-i],prob=Probs,replace=T)
#       Y[ind,i,tt] <- Y[i,ind,tt]  <- 1
#     }
#   }
# }

Bout <- Bin <- numeric(N)
w <- matrix(0,n,N)
s2 <- t2 <- numeric(N)
# AccRate <- matrix(0,N,3+n*TT); AccRate[1,] <- rep(NA,3+n*TT)
# colnames(AccRate) <- c("Bin","Bout","weights",
#                        paste("X",rep(1:n,TT),rep(1:TT,each=n),sep=","))
AccRate <- numeric(3+ncov+n*TT)
names(AccRate) <- c("Bin","Bout","Bin","weights",
                    paste("X",rep(1:n,TT),rep(1:TT,each=n),sep=","))


### covariates
B1 <- array(dim = c(ncov, N))
B1[, 1] <- rep(0, ncov)
sigmaZ <- array(dim = c(ncov, N))
sigmaZ[,1] <- rep(0.001, ncov)
phi <- array(dim = c(ncov, N))
phi[, 1] <- rep(0, ncov)

###
###Weights
###
for(tt in 1:TT){
  w[,1] <- w[,1] + apply(Y[,,tt],1,sum) + 
    apply(Y[,,tt],2,sum)
}
w[,1] <- w[,1]/sum(Y)/2
if(sum(w==0)>0){
  w[,1] <- w[,1]+1e-5
  w[,1] <- w[,1]/sum(w[,1])
}
w[,1] <- w[,1]/sum(w[,1])

###
###Initial Latent Positions (GMDS, Sarkar and Moore, 2005)
###
dissim <- array(0,dim=dim(Y))
for(tt in 1:TT) dissim[,,tt] <- shortest.paths(graph=graph.adjacency(Y[,,tt]),mode="all")
dissim[which(dissim==Inf,arr.ind=TRUE)] <- max(c(dissim[which(dissim!=Inf,arr.ind=TRUE)]))
X <- list()
X[[1]] <- array(0,c(p,TT,n))
X[[1]][,1,] <- t(cmdscale(d=dissim[,,1],k=p))
temp.lambda <- 10
H <- matrix(-1/n,n,n)+diag(n)

for(tt in 2:TT){
  temp <- 1/(1+temp.lambda)*H%*%(-1/2*dissim[,,tt]^2)%*%H +
    temp.lambda/(1+temp.lambda)*t(X[[1]][,tt-1,])%*%X[[1]][,tt-1,]
  temp <- eigen(temp)
  X[[1]][,tt,] <- t(temp$vectors[,1:p]%*%diag(temp$values[1:p])^(1/2))
  X[[1]][,tt,] <- t(vegan::procrustes(X=t(X[[1]][,tt-1,]),Y=t(X[[1]][,tt,]),scale=FALSE)$Yrot)
}

###
###Initial \beta_{IN} and \beta_{OUT}; scale latent positions
###
initialize.wrap <- function(x){
  -c.initialize1(X[[1]],c(n,p,TT),Y,XSCALE=1/n,
                 BETAIN=x[1],BETAOUT=x[2],w[,1])
}
initialize.grad.wrap <- function(x){
  -c.initialize1.grad(X[[1]],c(n,p,TT),Y,XSCALE=1/n,
                      BETAIN=x[1],BETAOUT=x[2],w[,1])[2:3]
}
Optim <- optim(par=c(1,1),fn=initialize.wrap,
               gr=initialize.grad.wrap,method="BFGS")
X[[1]] <- X[[1]]/n
Bin[1] <- max(Optim$par[1],1e-4)
Bout[1] <- max(Optim$par[2],1e-4)

Xit0 <- t(X[[1]][,1,])
for(tt in 2:TT)Xit0 <- rbind(Xit0,t(X[[1]][,tt,]))
Xit0 <- Xit0 -
  kronecker(rep(1,n*TT),matrix(apply(Xit0,2,mean),1,p))

###
###Priors and further initializations:
###
nuIn <- Bin[1]
nuOut <- Bout[1]
xiIn <- xiOut <- 10
xiB1 <- rep(10, ncov)
nuB1 <- rep(0,ncov)
thetaZ <- 2
phiZ <- 10
nuPhi <- 0
xiPhi <- 10

#Tau^2
t2[1] <- sum(X[[1]][,1,]*X[[1]][,1,])/(n*p)
shapeT2 <- 2.05
scaleT2 <- (shapeT2-1)*t2[1]

#Sigma^2
s2[1] <- 0
for (ti in 1:(TT-1)){
  s2[1] <- sum((X[[1]][,ti+1,]-X[[1]][,ti,])^2)+s2[1]
}
s2[1] <- s2[1]/(TT-1)/n/p
# shapeS2 <- 9
# scaleS2 <- 1.5
shapeS2 <- 2.05
scaleS2 <- (shapeS2-1)*s2[1]


nuA = 0
xiA = 10
alpha = numeric(N)
alpha[1] <- 0

xi0 = 10
nu0 = 0
phi0 = 10
theta0 = 2
sigma0 = array(0.0001, dim = c(ncov, N))
mu0 = array(0, dim = c(ncov, N))

###
###log-likelihood approximation subsampling
###
if(llApprox){
  DEGREE <- array(0,c(n,2,TT))
  for(tt in 1:TT){
    DEGREE[,1,tt] <- colSums(Y[,,tt])
    DEGREE[,2,tt] <- rowSums(Y[,,tt])
  }
  dinmax <- max(DEGREE[,1,])
  doutmax <- max(DEGREE[,2,])
  ELIN <- array(0,c(n,dinmax,TT))
  ELOUT <- array(0,c(n,doutmax,TT))
  
  for(tt in 1:TT){
    for(i in 1:n){
      ind <- which(Y[,i,tt]==1)
      if(length(ind)>0){
        ELIN[i,1:length(ind),tt] <- ind
      }
      ind <- which(Y[i,,tt]==1)
      if(length(ind)>0){
        ELOUT[i,1:length(ind),tt] <- ind
      }
    }
  }
  
  n0<-max(n0,dinmax+10,doutmax+10)
  
  edgeList <- list()
  for(i in 1:n){
    edgeList[[i]] <- which(Y[i,,]==1|Y[,i,]==1,arr.ind=TRUE)[,1]
    edgeList[[i]] <- unique(edgeList[[i]])
  }
  
  SUBSEQ = matrix(0,n,n0)
  for(i in 1:n){
    #   SUBSEQ[i,] <- sample(c(1:n)[-i],size=n0,replace=FALSE) #Simple random sample
    nOnes <- round(length(edgeList[[i]])/n*n0) #stratified sampling
    if(length(edgeList[[i]])>0){ nOnes <- max(nOnes,1) }
    SUBSEQ[i,1:nOnes] <- sample(edgeList[[i]],size=nOnes,replace=FALSE)
    SUBSEQ[i,(nOnes+1):n0] <- sample(c(1:n)[-c(i,edgeList[[i]])],size=n0-nOnes,replace=FALSE)
  }
}

