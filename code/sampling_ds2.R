
### data loading


setwd("~/Library/CloudStorage/GoogleDrive-zxu9@nd.edu/My Drive/projects/DLSMcov/shandong")
load("shandong.RData")
Y.pre = shandong[[1]][1:3]
n = 71
TT = 3
Y = array(as.numeric(unlist(Y.pre)), dim=c(n, n, TT))

# Y[is.na(Y)] <-0
ncov = 2
C = array(NA, dim = c(n, ncov, TT))
for (t in 1:TT) {
  c1 <- shandong[[2]][[t]][, 'lone']
  c2 <- shandong[[2]][[t]][, 'extroversion']
  for (i in 1:n){
    C[i, 1, t] = c1[i]
    C[i, 2, t] = c2[i]
  }
}

# C[,1,] <- scale(C[,1,], scale = F, center = T)
# C[,2,] <- scale(C[,2,], scale = F, center = T)
# C[is.na(C)] <- 0
C[,1,] <- C[,1,] - mean(C[,1,], na.rm = T)
C[,2,] <- C[,2,] - mean(C[,2,], na.rm = T)

save(Y, C, file = "ds2min.RData")

### setup

#Number of MCMC iterations
N = 1000000
#Dimension of the Euclidean latent space
p=2
#Use log likelihood approximation (BOOLEAN)?
#If TRUE, how large a subsample n0?
llApprox = FALSE
if(llApprox) n0 = 100
#Are there missing edges?
MissData = T
#If TRUE, construct Missing: Missing[[tt]]= c(1,2,3) => at time tt we have no row data on actors 1,2&3 
if(MissData){
  Missing <- list() #Enter as list manually, or if NAs are in data run the following:
  mindex = which(is.na(Y),arr.ind=TRUE)
  for(tt in 1:dim(Y)[3]){
    Missing[[tt]] = unique(mindex[which(mindex[,3]==tt),1])
  }
  Y[mindex] = 0
}
cmindex = which(is.na(C),arr.ind=TRUE)
C[is.na(C)] <- 0
cr <- range(C)


#MCMC tuning parameters
tuneX <-   0.0075
tuneBIO <- 0.1
tuneB1 <- c(0.5, 0.5, 0.5)
Kappa <-   175000
burnin = round(N/10)

require("igraph")
require("MCMCpack")
require("inline")
require("RcppArmadillo")
require("vegan")

source("functions_noa.R")
# source("sim.R")




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
nchain = 1


###Run MCMC


# pb <- txtProgressBar(min=2,max=500,style=3)
system.time({
  for (chain in 1:nchain){
    if (chain == 2){
      soure("initialize2.R")
    }else{
      source("initialize.R")
    }
    
    load("sdres_newpr_centerts.RData")
    nt <- 500000
    Bin[1:nt] <- chaindat$Bin
    Bout[1:nt] <- chaindat$Bout
    B1[, 1:nt] <- chaindat$B1
    s2[1:nt] <- chaindat$s2
    t2[1:nt] <- chaindat$t2
    mu0[,1:nt] <- chaindat$mu0
    sigma0[,1:nt] <- chaindat$sigma0
    phi[,1:nt] <- chaindat$phi
    sigmaZ[,1:nt] <- chaindat$sigmaZ
    AccRate <- chaindat$AccRate
    X <- chaindat$X
    w[,1:nt] <- chaindat$w
    C = chaindat$C
    Y = chaindat$Y
    
    set.seed(1*chain)
    pb <- txtProgressBar(min=2,max=N,style=3)
    for(it in 825690:N){
      RN <- rnorm(n*TT*p)
      RNBIO <- rnorm(2)
      RNB1 <- rnorm(ncov)
      RNA <- 0.02*rnorm(1)
      
      if(llApprox){
        if(it%%100==0){
          SUBSEQ = matrix(0,n,n0)
          for(i in 1:n){
            nOnes <- round(length(edgeList[[i]])/n*n0) #stratified sampling
            if(length(edgeList[[i]])>0){ nOnes <- max(nOnes,1) }
            SUBSEQ[i,1:nOnes] <- sample(edgeList[[i]],size=nOnes,replace=FALSE)
            SUBSEQ[i,(nOnes+1):n0] <- sample(c(1:n)[-c(i,edgeList[[i]])],size=n0-nOnes,replace=FALSE)
          }
        }
      }
      if(llApprox){
        Draws <- c.update2(X[[it-1]],c(n,p,TT,dinmax,doutmax),tuneX,Y,
                           Bin[it-1],Bout[it-1],tuneBIO,w[,it-1],
                           t2[it-1],s2[it-1],xiIn,xiOut,nuIn,
                           nuOut,CAUCHY=0,RN,RNBIO,ELOUT,ELIN,SUBSEQ,DEGREE)
      }else{
        Draws <- c.update1(X[[it-1]],c(n,p,TT,1),tuneX,Y,
                           Bin[it-1],Bout[it-1],tuneBIO,w[,it-1],
                           t2[it-1],s2[it-1],xiIn,xiOut,nuIn,
                           nuOut,CAUCHY=0,RN,RNBIO, RNB1, tuneB1, B1[,it-1], xiB1, nuB1, C, ncov)
      }
      X[[it]] <- Draws[[1]]
      Bin[it] <- Draws[[2]]
      Bout[it] <- Draws[[3]]
      B1[,it] <- Draws[[4]]
      AccRate <- AccRate + Draws[[5]]
      
      if(it==burnin){
        Xit0 <- t(X[[it]][,1,])
        for(tt in 2:TT) Xit0 <- rbind(Xit0,t(X[[it]][,tt,]))
      }
      if(it>burnin){
        XitCentered <- t(X[[it]][,1,])
        for(tt in 2:TT) XitCentered <- rbind(XitCentered,t(X[[it]][,tt,]))
        procr <- vegan::procrustes(X=Xit0,Y=XitCentered,scale=FALSE)$Yrot
        for(tt in 1:TT){
          X[[it]][,tt,] <- t(procr[((tt-1)*n+1):(tt*n),])
        }
      }
      if(it < N) X[[it+1]] <- X[[it]]
      
      #------------------Step 2--------------------------------
      Draws1 <- c.t2s2Parms(X[[it]],c(n,p,TT,1),shapeT2,
                            shapeS2,scaleT2,scaleS2)
      t2[it] <- rinvgamma(1,shape=Draws1[[1]],scale=Draws1[[2]])
      s2[it] <- rinvgamma(1,shape=Draws1[[3]],scale=Draws1[[4]])
      
      for(q in 1:ncov){
        Drawsphisig <- c.phisigmazParms(X[[it]], c(n, p, TT, 1), thetaZ, phiZ, xiPhi, nuPhi, sigmaZ[q,it-1], C[,q,], phi[q,it-1], xi0, nu0, theta0, phi0, sigma0[q,it-1], mu0[q, it-1])
        sigmaZ[q, it] <- rinvgamma(1, shape=Drawsphisig[[1]], scale=Drawsphisig[[2]])
        phi[q, it] <- rnorm(1, mean=Drawsphisig[[3]], sd=sqrt(Drawsphisig[[4]]))
        if(is.nan(phi[q, it])){
          phi[q, it] <- phi[q, it-1]
          print("phi is nan")
        }
        if(is.nan(sigmaZ[q, it]) | is.infinite(sigmaZ[q,it])){
          sigmaZ[q, it] <- sigmaZ[q, it-1]
          print("sigmaz is nan")
        }
        sigma0[q, it] <- rinvgamma(1, shape = Drawsphisig[[5]], scale = Drawsphisig[[6]])
        mu0[q,it] <- rnorm(1, mean=Drawsphisig[[7]], sd=sqrt(Drawsphisig[[8]]))
        # print(Drawsphisig)
      }
      
      #------------------Step 3--------------------------------
      
      w[,it] <- rdirichlet(1,alpha=Kappa*w[,it-1])
      if(llApprox){
        Draws2 <- c.WAccProb2(X[[it]],c(n,p,TT,dinmax,doutmax),Y,
                              Bin[it],Bout[it],Kappa,w[,it-1],w[,it],
                              ELOUT,ELIN,SUBSEQ,DEGREE)
      }else{
        Draws2 <- c.WAccProb1(X[[it]],c(n,p,TT,1),Y,
                              Bin[it],Bout[it],Kappa,w[,it-1],w[,it], C, B1[,it], ncov)
      }
      w[,it] <- Draws2[[1]]
      AccRate[ncov + 2 + 1] <- AccRate[ncov + 2 + 1] + Draws2[[2]]
      
      #------------------Step 4--------------------------------
      
      if(MissData){
        for(tt in 1:TT){
          Yimp <- c.missing(X[[it]],c(n,p,TT),MMM=Missing[[tt]]-1,Y,Ttt=tt,
                            BETAIN=Bin[it],BETAOUT=Bout[it],WW=w[,it], B1[,it-1], C, ncov)
          Y[mindex] = Yimp[mindex]
        }
        
        
        Cimp <- c.missingcov(c(n,p,TT), C, phi[,it-1], sigmaZ[, it-1], ncov, mu0[,it-1], sigma0[,it-1], rnorm(n*TT*ncov))
        C[cmindex] = Cimp[cmindex]
        C[C>cr[2]] = cr[2]
        C[C<cr[1]] = cr[1]
      }
      
      
      if(it%%100000==0){
        chaindat <- list(Bin = Bin, Bout = Bout,
                         B1 = B1, s2 = s2,
                         t2 = t2, mu0 = mu0,
                         sigma0 = sigma0, phi = phi,  sigmaZ = sigmaZ,
                         AccRate = AccRate, X = X, w = w, C = C, Y = Y)
        save(chaindat, file = paste0("sdres_newpr_centerts.RData"))
      }
      setTxtProgressBar(pb,it)
    }
    close(pb)
  }
})
AccRate[1:10]/(it-1)
summary(AccRate[-c(1:3)]/(it-1))

load("sdres1.RData")
plot(B1[1, ], type = "l")
plot(B1[2,], type = "l")


plot(Bin, type = "l")
plot(Bout, type = "l")

plot(t2, type = "l")
plot(s2, type = "l")
plot(mu0[1,], type = 'l')
plot(mu0[2,], type = 'l')
plot(sigma0[1,], type = 'l')
plot(sigma0[2,], type = 'l')
plot(phi[1,], type = "l")
plot(phi[2,], type = "l")
plot(sigmaZ[1,], type = "l")
plot(sigmaZ[2,], type = "l")

rowMeans(phi)
rowMeans(sigmaZ)
rowMeans(mu0[,5000:10000])
rowMeans(sigma0[,5000:1000])
apply(C[,,1], 2, mean)
apply(C[,,1], 2, var)
rowMeans(B1)
mean(Bin)
mean(Bout)

nst = 50000
ntemp = 250000
geweke.diag(Bin[nst:ntemp][seq(1, ntemp-nst, 10)], 0.1, 0.5)
geweke.diag(Bout[nst:ntemp][seq(1, ntemp-nst, 10)], 0.1, 0.5)
geweke.diag(B1[1,nst:ntemp][seq(1, ntemp-nst, 10)], 0.1, 0.5)
geweke.diag(B1[2,nst:ntemp][seq(1, ntemp-nst, 10)], 0.1, 0.5)
geweke.diag(s2[nst:ntemp][seq(1, ntemp-nst, 10)], 0.1, 0.5)
geweke.diag(t2[nst:ntemp][seq(1, ntemp-nst, 10)], 0.1, 0.5)


