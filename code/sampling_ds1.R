# 
# Extension of Sewell & Chen (2015) to include dynamic non-network covariate evolution. Based on the dynamic latent space model code from https://publish.illinois.edu/yuguo/software/ and https://sewell.lab.uiowa.edu/software
# 
# Reference: Sewell, D. K., & Chen, Y. (2015). Latent space models for dynamic networks. Journal of the american statistical association, 110(512), 1646-1657.

### data loading

load("guangzhou.RData")
Y.pre = guangzhou$net[2:31]
n = 47
TT = 30
Yold = array(as.numeric(unlist(Y.pre)), dim=c(n, n, TT))
Yold[Yold<3] <-0
Yold[Yold>=3] <-1

for (t in 1:TT){
  print(mean(is.na(Yold[,,t])))
}

ncov = 3
Cold = array(NA, dim = c(n, ncov, TT))

for (t in 1:TT) {
  c1 <- scale(guangzhou$covar[[t+1]]$social, center = F, scale = F)
  c2 <- scale(guangzhou$covar[[t+1]]$happy, center = F, scale = F)
  c3 <- scale(guangzhou$covar[[t+1]]$pressure, center = F, scale = F)
  for (i in 1:n){
    Cold[i, 1, t] = c1[i]
    Cold[i, 2, t] = c2[i]
    Cold[i, 3, t] = c3[i]
  }
}

# Cold[,1,] <- scale(Cold[,1,], center = T, scale = F)
# Cold[,2,] <- scale(Cold[,2,], center = T, scale = F)
# Cold[,3,] <- scale(Cold[,3,], center = T, scale = F)

indrm <- c(1, 11, 21, 31, 41, 14, 19, 23, 25, 26, 42)
n = 47 - length(indrm)
Y = array(0, dim=c(n, n, TT))
C = array(NA, dim = c(n, ncov, TT))

for (t in 1:TT){
  Y[,,t] = Yold[-indrm, -indrm, t]
  C[,,t] = Cold[-indrm,,t]
}

C[,1,] = C[,1,] - mean(C[,1,], na.rm = T)
C[,2,] = C[,2,] - mean(C[,2,], na.rm = T)
C[,3,] = C[,3,] - mean(C[,3,], na.rm = T)
# Y[is.na(Y)] <-0



### setup

#Number of MCMC iterations
N = 100000
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
tuneBIO <- 0.1 # 0.5 -> 0.1
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
thin = 10




###Run MCMC
# 
# nt <- 20000
# Bin[1:nt] <- chaindat$Bin[(1:nt)*thin]
# Bout[1:nt] <- chaindat$Bout[(1:nt)*thin]
# B1[, 1:nt] <- chaindat$B1[,(1:nt)*thin]
# s2[1:nt] <- chaindat$s2[(1:nt)*thin]
# t2[1:nt] <- chaindat$t2[(1:nt)*thin]
# mu0[,1:nt] <- chaindat$mu0[,(1:nt)*thin]
# sigma0[,1:nt] <- chaindat$sigma0[,(1:nt)*thin]
# phi[,1:nt] <- chaindat$phi[,(1:nt)*thin]
# sigmaZ[,1:nt] <- chaindat$sigmaZ[,(1:nt)*thin]
# AccRate <- chaindat$AccRate
# X <- chaindat$X[(1:nt)*thin]
# w[,1:nt] <- chaindat$w[,(1:nt)*thin]
# C = chaindat$C
# Y = chaindat$Y


# pb <- txtProgressBar(min=2,max=500,style=3)
system.time({
  for (chain in 1:nchain){
    if (chain == 2){
      source("initialize2.R")
    }else{
      source("initialize.R")
    }
    
    # load("/Users/ziqianxu/Desktop/gzresalt_newpr1.RData")
    # nt <- 50000
    # Bin[1:nt] <- chaindat$Bin
    # Bout[1:nt] <- chaindat$Bout
    # B1[, 1:nt] <- chaindat$B1
    # s2[1:nt] <- chaindat$s2
    # t2[1:nt] <- chaindat$t2
    # mu0[,1:nt] <- chaindat$mu0
    # sigma0[,1:nt] <- chaindat$sigma0
    # phi[,1:nt] <- chaindat$phi
    # sigmaZ[,1:nt] <- chaindat$sigmaZ
    # AccRate <- chaindat$AccRate
    # X <- chaindat$X
    # w[,1:nt] <- chaindat$w
    # C = chaindat$C
    # Y = chaindat$Y
    
    
    # set.seed(1*chain)
    pb <- txtProgressBar(min=2,max=N,style=3)
    for(its in 106150:(N*thin)){ #(2*thin):(N*thin)
      it <- floor(its/thin)
      RN <- rnorm(n*TT*p)
      RNBIO <- rnorm(2)
      RNB1 <- rnorm(ncov)
      RNA <- 0.02*rnorm(1)
      
      
      if(its %% thin == 0){
        Draws <- c.update1(X[[it-1]],c(n,p,TT,1),tuneX,Y,
                           Bin[it-1],Bout[it-1],tuneBIO,w[,it-1],
                           t2[it-1],s2[it-1],xiIn,xiOut,nuIn,
                           nuOut,CAUCHY=0,RN,RNBIO, RNB1, tuneB1, B1[,it-1], xiB1, nuB1, C, ncov)
      }else{
        Draws <- c.update1(X[[it]],c(n,p,TT,1),tuneX,Y,
                           Bin[it],Bout[it],tuneBIO,w[,it],
                           t2[it],s2[it],xiIn,xiOut,nuIn,
                           nuOut,CAUCHY=0,RN,RNBIO, RNB1, tuneB1, B1[,it], xiB1, nuB1, C, ncov)
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
        if (its %% thin == 0){
          Drawsphisig <- c.phisigmazParms(X[[it]], c(n, p, TT, 1), thetaZ, phiZ, xiPhi, nuPhi, sigmaZ[q,it-1], C[,q,], phi[q,it-1], xi0, nu0, theta0, phi0, sigma0[q,it-1], mu0[q, it-1])
        }else{
          Drawsphisig <- c.phisigmazParms(X[[it]], c(n, p, TT, 1), thetaZ, phiZ, xiPhi, nuPhi, sigmaZ[q,it], C[,q,], phi[q,it], xi0, nu0, theta0, phi0, sigma0[q,it], mu0[q, it])
        }
        
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
        
        
        
        
      }
      
      #------------------Step 3--------------------------------
      
      
      if (its %% thin == 0){
        w[,it] <- rdirichlet(1,alpha=Kappa*w[,it-1])
        Draws2 <- c.WAccProb1(X[[it]],c(n,p,TT,1),Y,
                              Bin[it],Bout[it],Kappa,w[,it-1],w[,it], C, B1[,it], ncov)
      }else{
        Draws2 <- c.WAccProb1(X[[it]],c(n,p,TT,1),Y,
                              Bin[it],Bout[it],Kappa,w[,it],rdirichlet(1,alpha=Kappa*w[,it]), C, B1[,it], ncov)
      }
      
      
      w[,it] <- Draws2[[1]]
      AccRate[ncov + 2 + 1] <- AccRate[ncov + 2 + 1] + Draws2[[2]]
      
      
      
      #------------------Step 4--------------------------------
      
      if(MissData){
        for(tt in 1:TT){
          Yimp <- c.missing(X[[it]],c(n,p,TT),MMM=Missing[[tt]]-1,Y,Ttt=tt,
                            BETAIN=Bin[it],BETAOUT=Bout[it],WW=w[,it], B1[,it], C, ncov)
          Y[mindex] = Yimp[mindex]
        }
        
        
        Cimp <- c.missingcov(c(n,p,TT), C, phi[,it], sigmaZ[, it], ncov, mu0[,it], sigma0[,it], rnorm(n*TT*ncov))
        C[cmindex] = Cimp[cmindex]
        C[C>cr[2]] = cr[2]
        C[C<cr[1]] = cr[1]
      }
      
      chain = 1
      if(its%%100000==0){
        chaindat <- list(Bin = Bin, Bout = Bout,
                         B1 = B1, s2 = s2,
                         t2 = t2, mu0 = mu0,
                         sigma0 = sigma0, phi = phi,  sigmaZ = sigmaZ,
                         AccRate = AccRate, X = X, w = w, C = C, Y = Y)
        save(chaindat, file = paste0("/Users/ziqianxu/Desktop/gzresalt_center", chain, ".RData"))
      }
      setTxtProgressBar(pb,it)
    }
    close(pb)
  }
})
AccRate[1:10]/(it-1)
summary(AccRate[-c(1:3)]/(it-1))

# s2 initial 0.001
plot(B1[1,], type = "l")
plot(B1[2,], type = "l")
plot(B1[3,], type = "l")

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

ntemp = 55500
nst = 15000
thinst = 5
geweke.diag(Bin[nst:ntemp][seq(1, ntemp-nst, thinst)])
geweke.diag(Bout[nst:ntemp][seq(1, ntemp-nst, thinst)])
geweke.diag(t2[nst:ntemp][seq(1, ntemp-nst, thinst)])
geweke.diag(s2[nst:ntemp][seq(1, ntemp-nst, thinst)])
geweke.diag(B1[1,nst:ntemp][seq(1, ntemp-nst, thinst)])
geweke.diag(B1[2,nst:ntemp][seq(1, ntemp-nst, thinst)])
geweke.diag(B1[3,nst:ntemp][seq(1, ntemp-nst, thinst)])

