

###### visualization of raw data
load("guangzhou.RData")
Y.pre = guangzhou$net[2:31]
n = 47
TT = 30
Yold = array(as.numeric(unlist(Y.pre)), dim=c(n, n, TT))
Yold[Yold<3] <-0
Yold[Yold>=3] <-1

ncov = 3
Cold = array(NA, dim = c(n, ncov, TT))
for (t in 1:TT) {
  c1 <- scale(guangzhou$covar[[t+1]]$social, center = T, scale = F)
  c2 <- scale(guangzhou$covar[[t+1]]$happy, center = T, scale = F)
  c3 <- scale(guangzhou$covar[[t+1]]$pressure, center = T, scale = F)
  for (i in 1:n){
    Cold[i, 1, t] = c1[i]
    Cold[i, 2, t] = c2[i]
    Cold[i, 3, t] = c3[i]
  }
}

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


for (t in 1:TT) {
  cat("\nCorrelation matrix at time", t, ":\n")
  corr_mat <- cor(C[,,t], use = "pairwise.complete.obs")
  print(corr_mat)
  
  # Exclude diagonal by setting it to NA
  diag(corr_mat) <- NA
  
  # Find min and max
  max_corr <- max(corr_mat, na.rm = TRUE)
  min_corr <- min(corr_mat, na.rm = TRUE)
  
  cat("Max correlation (off-diagonal):", round(max_corr, 3), "\n")
  cat("Min correlation (off-diagonal):", round(min_corr, 3), "\n")
}

library(reshape2)

library(ggplot2)
library(reshape2)

names = c('social', 'happy', 'pressure')
for (cov_idx in 1:ncov) {
  cov_mat <- C[, cov_idx, ]
  colnames(cov_mat) <- paste0("T", 1:TT)
  
  # Compute temporal correlation
  corr_mat <- cor(cov_mat, use = "pairwise.complete.obs")
  cat("\nTemporal correlation for covariate", cov_idx, ":\n")
  print(corr_mat)
  
  # Melt for ggplot
  corr_melt <- melt(corr_mat)
  colnames(corr_melt) <- c("Time1", "Time2", "Correlation")
  
  # Plot
  ggplot(corr_melt, aes(x = Time1, y = Time2, fill = Correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                         limits = c(-1,1), name = "Correlation") +
    theme_minimal() +
    labs(title = paste("Temporal correlation heatmap for covariate", names[cov_idx]),
         x = "Time", y = "Time") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

library(car)

for (t in 1:TT) {
  df <- as.data.frame(C[,,t])
  colnames(df) <- c("lone", "extroversion")
  cat("\nVIF at time", t, ":\n")
  vif_model <- lm(rnorm(n) ~ ., data = df)  # dummy dependent variable
  print(vif(vif_model))
}




# cov visualization

pdf(file="figs/covar.pdf", width = 8, height = 4)
par(mfrow = c(1, 1))
c1mean <- apply(C[,1,], 2, mean, na.rm = T)
plot(c1mean, type = 'l' , col = 'black', lwd = 2,
     xlab = 'week', ylab = 'covariate value', ylim = c(-1, 1))
c2mean <- apply(C[,2,], 2, mean, na.rm = T)
lines(c2mean, type = 'l', col = '#5a95a6', lwd = 2)
c3mean <- apply(C[,3,], 2, mean, na.rm = T)
lines(c3mean, type = 'l', col = '#913c8d', lwd = 2)
legend("bottomright", legend = c('social', "happiness", "pressure"), 
       col = c('black', '#5a95a6', '#913c8d'), lwd = 2, lty = 1)
dev.off()






######################### data 
load("gzresalt_centerts.RData")  # tune BIO = 0.5 tune b1 = 0.5
Bin = chaindat$Bin; Bout = chaindat$Bout; w = chaindat$w; X = chaindat$X; phi = chaindat$phi; s2 = chaindat$s2; t2 = chaindat$t2; sigmaZ = chaindat$sigmaZ; B1 = chaindat$B1; alpha = chaindat$alpha; sigma0 = chaindat$sigma0; mu0 = chaindat$mu0; C = chaindat$C; Y = chaindat$Y
N = 70000
burnin = 20000
thin = 5
select = seq(burnin, N, by = thin)
n =  36
TT = 30
p = 2



effectiveSize(Bin[select])
effectiveSize(as.mcmc(Bout[select]))
effectiveSize(as.mcmc(t2[select]))
effectiveSize(as.mcmc(s2[select]))
effectiveSize(as.mcmc(B1[1, select]))
effectiveSize(as.mcmc(B1[2, select]))
effectiveSize(as.mcmc(B1[3, select]))
effectiveSize(as.mcmc(phi[1, select]))
########### geweke
library(coda)
f1 = 0.1
f2 = 0.5

geweke_results <- list(
  Bin = geweke.diag(Bin[select], f1, f2),
  Bout = geweke.diag(Bout[select], f1, f2),
  s2 = geweke.diag(s2[select], f1, f2),
  t2 = geweke.diag(t2[select], f1, f2),
  B1_1 = geweke.diag(B1[1, select], f1, f2),
  B1_2 = geweke.diag(B1[2, select], f1, f2),
  B1_3 = geweke.diag(B1[3, select], f1, f2),
  phi_1 = geweke.diag(phi[1, select], f1, f2),
  phi_2 = geweke.diag(phi[2, select], f1, f2),
  phi_3 = geweke.diag(phi[3, select], f1, f2),
  sigmaZ_1 = geweke.diag(sigmaZ[1, select], f1, f2),
  sigmaZ_2 = geweke.diag(sigmaZ[2, select], f1, f2),
  sigmaZ_3 = geweke.diag(sigmaZ[3, select], f1, f2)
)

geweke_df <- data.frame(
  Parameter = names(geweke_results),
  Z_score = sapply(geweke_results, function(x) x$z) # Assuming geweke.diag returns a list with an element 'z'
)
geweke_df
library(knitr)
kable(geweke_df, booktabs = T, caption = "means", format = 'latex', digits = 3)

### autocor
autocorr(as.mcmc(Bin[select]))
autocorr(as.mcmc(s2[select]))

###### welch
heidel.diag(Bin[select], 0.3)
heidel.diag(Bout[select],0.3)
heidel.diag(s2[select])
heidel.diag(t2[select])
heidel.diag(B1[1,select])
heidel.diag(B1[2,select])
heidel.diag(B1[3,select])
heidel.diag(phi[1, select])
heidel.diag(phi[2, select])
heidel.diag(phi[3, select])
heidel.diag(sigmaZ[1, select])
heidel.diag(sigmaZ[2, select])
heidel.diag(sigmaZ[3, select])





#################### MCMC  plots 
xlabel = seq(0, 70000, 10000)
### auto
for (var in c("Bin", "Bout","s2", "t2")){
  pdf(paste0("figs/", var, ".pdf"), width = 8, height = 4)
  par(mfrow = c(1, 1))
  plot(chaindat[[var]][1:N], type = 'l', col = 'black', lwd = 1, xaxt = "n", xlab = 'iteration', ylab = var)
  axis(1, at = xlabel, labels = xlabel * 10)
  # plot(density(chaindat[[var]][select]), col = 'black', lwd = 2, xlab = var, main = '')
  dev.off()
}

ncov = 3
for (var in c("phi", "sigmaZ", "B1")){
  for (i in 1:ncov){
    pdf(paste0("figs/", var, i, ".pdf"), width = 8, height = 4)
    par(mfrow = c(1, 1))
    varn = var
    if (var == 'B1') {varn = 'B'}
    if (var == 'sigmaZ') {varn = 'sigmac'}
    plot(chaindat[[var]][i,1:N], type = 'l', col = 'black', lwd = 1,  xaxt = "n",xlab = 'iteration', ylab = paste0(varn, i))
    axis(1, at = xlabel, labels = xlabel * 10)
    # plot(density(chaindat[[var]][i,select]), col = 'black', lwd = 2, xlab = var, main = '')
    dev.off()
  }
}
### handmade
pdf("figs/Bin.pdf", width = 8, height = 4)
par(mfrow = c(1, 1))
plot(Bin, type = 'l', col = 'black', lwd = 0.5, xlab = 'iteration', ylab = var, ylim = c(-1.5,0))
# plot(density(Bin), col = 'black', lwd = 2, xlab = var, main = '')
dev.off()


#################### means & credible intervals
EBin <- mean(Bin[select])
EBout <- mean(Bout[select])
Ew <- apply(w[,][,select],1,mean)
EX <- array(0,dim=c(p,TT,n))
for(it in select){
  EX <- EX + X[[it]]
}
EX <- EX/length(select)
Ephi <- rowMeans(phi[,select])
Es2 <- mean(s2[select])
Et2 <- mean(t2[select])
EsigmaZ <- rowMeans(sigmaZ[,select])
EB1 <- rowMeans(B1[,select])
EA <- mean(alpha[select])
Es0 <- rowMeans(sigma0[,select])
Em0 <- rowMeans(mu0[,select])

ql = 0.025
qu = 0.975
qBin <- quantile(Bin[select],c(ql, qu))
qBout <- quantile(Bout[select],c(ql, qu))
qB11 <- quantile(B1[1,][select],c(ql, qu))
qB12 <- quantile(B1[2,][select],c(ql, qu))
qB13 <- quantile(B1[3,][select],c(ql, qu))
qphi1 <- quantile(phi[1,][select],c(ql, qu))
qphi2 <- quantile(phi[2,][select],c(ql, qu))
qphi3 <- quantile(phi[3,][select],c(ql, qu))
qs2 <- quantile(s2[select],c(ql, qu))
qt2 <- quantile(t2[select],c(ql, qu))
qsigmaZ1 <- quantile(sigmaZ[1,][select],c(ql, qu))
qsigmaZ2 <- quantile(sigmaZ[2,][select],c(ql, qu))
qsigmaZ3 <- quantile(sigmaZ[3,][select],c(ql, qu))
qA <- quantile(alpha[select],c(ql, qu))
qs01 <- quantile(sigma0[1,][select],c(ql, qu))
qs02 <- quantile(sigma0[2,][select],c(ql, qu))
qs03 <- quantile(sigma0[3,][select],c(ql, qu))
qm01 <- quantile(mu0[1,][select],c(ql, qu))
qm02 <- quantile(mu0[2,][select],c(ql, qu))
qm03 <- quantile(mu0[3,][select],c(ql, qu))
Es0;Em0



es <- data.frame(mean = c( EB1[1], EB1[2], EB1[3], EBin, EBout, Es2, Et2,
                           Ephi[1], Ephi[2], Ephi[3], 
                           EsigmaZ[1], EsigmaZ[2], EsigmaZ[3]))
es$lower <- c( qB11[1], qB12[1], qB13[1], 
               qBin[1], qBout[1], qs2[1], qt2[1], 
               qphi1[1], qphi2[1], qphi3[1], qsigmaZ1[1], qsigmaZ2[1], qsigmaZ3[1])
es$upper <- c( qB11[2], qB12[2], qB13[2], 
               qBin[2], qBout[2], qs2[2], qt2[2], 
               qphi1[2], qphi2[2], qphi3[2], qsigmaZ1[2], qsigmaZ2[2], qsigmaZ3[2])
es$var <- c('B1_1', 'B1_2', 'B1_3', 'Bin', 'Bout', 's2', 't2', 'phi_1', 'phi_2', 'phi_3', 'sigmaZ_1', 'sigmaZ_2', 'sigmaZ_3')
es

kable(es[c(4, 1:3)], booktabs = T, caption = "means", format = 'latex', digits = 3)
library(kableExtra)
kable(es[,c(4, 1:3)], format = "html", caption = "Formatted Table for PowerPoint") %>%
  kable_styling(full_width = FALSE, position = "center") %>%
  row_spec(0, bold = TRUE) # Make the header bold



##################### ar1 residual autocor


eps <- array(NA, dim = c(n, ncov, TT))

for (k in 1:ncov) {
  for (i in 1:n) {
    eps[i, k, 1] <- C[i, k, 1] - Em0[k]
    for (t in 2:TT) {
      eps[i, k, t] <- C[i, k, t] - Ephi[k] * C[i, k, t-1]
    }
  }
}

par(mfrow = c(ncov, 1))
names = c('social', 'happy', 'pressure')
for (k in 1:ncov) {
  acf(as.vector(eps[, k, ]), main = paste("ACF residuals for ", names[k]))
}


################## plot of predictors against Y

library(ggplot2)
dists <- array(NA, dim = c(n, n, TT))
d1 <- array(NA, dim = c(n, n, TT))
d2 <- array(NA, dim = c(n, n, TT))
for(tt in 1:TT){
  for(i in 1:n){
    for(j in 1:n){
      dists[i,j,tt] <- sqrt(t(EX[,tt,i]-EX[,tt,j])%*%(EX[,tt,i]-EX[,tt,j]))
      d1[i,j,tt] = 1-dists[i,j,tt]/Ew[i]
      d2[i,j,tt] = 1-dists[i,j,tt]/Ew[j]
    }
  }
}



diag(dists[,,1]) <- NA
dist.flat <- array(dists[,,1], dim = c(n*n))
Y.flat <- array(Y[,,1], dim = c(n*n))
C.flat1 <- array(C[,1,1], dim = c(n*n))
C.flat2 <- array(C[,2,1], dim = c(n*n))
C.flat3 <- array(C[,3,1], dim = c(n*n))
d1.flat <- array(d1[,,1], dim = c(n*n))
d2.flat <- array(d2[,,1], dim = c(n*n))
data <- data.frame(dist = dist.flat, 
                   Y = Y.flat, 
                   C1 = C.flat1, 
                   C2 = C.flat2, 
                   C3 = C.flat3, 
                   d1= d1.flat, 
                   d2 = d2.flat)
pdf('figsnew/eda1.pdf', width = 8, height = 6)
library(RColorBrewer)
cbp1 <- c("black", "#5a95a6")
data$friend = ifelse(data$Y == 1, "yes", "no")

ggplot(data, aes(x=dist))+
  geom_density(aes(color=friend)) + 
  theme_minimal(base_size = 18) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
xlab("distance between latent positions") +
  ylab("density")
dev.off()


pdf('figsnew/eda2.pdf', width = 8, height = 6)
ggplot(data, aes(x=C2))+
  geom_density(aes(color=friend), bw = 0.8) + 
  theme_minimal(base_size = 18) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_color_manual(values = cbp1) + 
  xlab("sender's happiness score - receiver's score") +
  ylab("density ")
dev.off()

pdf('figsnew/eda3.pdf', width = 8, height = 6)
ggplot(data, aes(x=C3))+
  geom_density(aes(color=friend), bw = 1) + 
  theme_minimal(base_size = 18) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_color_manual(values = cbp1) + 
  xlab("sender's pressure score - receiver's score") +
  ylab("density ")
dev.off()

pdf('figsnew/eda4.pdf', width = 8, height = 6)
cbp1 <- c("black", "#5a95a6")
data$friend = ifelse(data$Y == 1, "yes", "no")
ggplot(data, aes(x=d1))+
  geom_density(aes(color=friend), binwith = 1) + 
  theme_minimal(base_size = 18) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_color_manual(values = cbp1) + 
  xlab(expression(1-d[ijt]/r[i])) +
  ylab("density") + xlim (-2,2)
dev.off()

pdf('figsnew/eda5.pdf', width = 8, height = 6)
cbp1 <- c("black", "#5a95a6")
data$friend = ifelse(data$Y == 1, "yes", "no")
ggplot(data, aes(x=d2))+
  geom_density(aes(color=friend), binwith = 1) + 
  theme_minimal(base_size = 18) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  scale_color_manual(values = cbp1) + 
  xlab(expression(1-d[ijt]/r[j])) +
  ylab("density") + xlim (-2,2)
dev.off()



xl <- c(1.15*min(EX[1,,]),1.15*max(EX[1,,]))
yl <- c(1.15*min(EX[2,,]),1.15*max(EX[2,,]))
lims <- range(c(xl,yl))

pdf(file="figs/lps_wp.pdf",width=8,height=8)
par(mar=c(2,2,4,2))
color_scale <- colorRampPalette(c("#addfed", "#0f414f"))
C_min <- min(C[,1,], na.rm=TRUE)
C_max <- max(C[,1,], na.rm=TRUE)
C_scaled <- (C[,1,] - C_min) / (C_max - C_min + 1e-6) # Avoid division by zero
colors <- color_scale(100) # Generate 100 colors

for(tt in 1:TT){
  # Determine colors based on C values
  col_vec <- colors[round(C_scaled[ , tt] * 99) + 1]  # Map to color index
  
  plot(t(EX[,tt,]), xlim=lims, ylim=lims, xlab="", ylab="", 
       pch=16, cex=1, col=col_vec, xaxt="n", yaxt="n", main="")
  
  # if(tt > 1) arrows(EX[1,tt-1,], EX[2,tt-1,], EX[1,tt,], EX[2,tt,], 
  #                   length=0, col='blue')
  
  if(tt < TT) par(new=TRUE)
}
for(tt in 1:TT){
  if(tt > 1){
    arrows(EX[1,tt-1,], EX[2,tt-1,], EX[1,tt,], EX[2,tt,], 
           length=0, col='blue')
  } 
  if (tt == TT){
    arrows(EX[1,tt-1,], EX[2,tt-1,], EX[1,tt,], EX[2,tt,], 
           length=0.1, col='blue')
  }
  if(tt < TT) par(new=TRUE)
}

dev.off()





pdf(file="figs/lpsff_wp.pdf",width=8,height=8)
par(mar=c(2,2,4,2))
col_vec <- colors[round(C_scaled[ , 1] * 99) + 1] 
plot(t(EX[,1,]),xlim=lims,ylim=lims,xlab="",ylab="",
     pch=16,cex=1,xaxt="n",yaxt="n",col = col_vec,
     main="")
par(new=TRUE)
col_vec <- colors[round(C_scaled[, TT] * 99) + 1] 
plot(t(EX[,TT,]),xlim=lims,ylim=lims,xlab="",ylab="",
     pch=16,cex=1,xaxt="n",yaxt="n", col = col_vec,
     main="")
arrows(EX[1,1,],EX[2,1,],EX[1,TT,],EX[2,TT,],length=0.05, col = 'blue', code = 2)
dev.off()



################## AUC calculation #################
# Y1 = Y
# C1 = C
# load("guangzhou.RData")
# Y.pre = guangzhou$net[2:31]
# n = 47
# TT = 30
# Y = array(as.numeric(unlist(Y.pre)), dim=c(n, n, TT))
# Y[Y<3] <-0
# Y[Y>=3] <-1
# # Y[is.na(Y)] <-0
# ncov = 3
# C = array(NA, dim = c(n, ncov, TT))
# for (t in 1:TT) {
#   for (i in 1:n){
#     # C[i, 1, t] = scale(guangzhou$covar[[t+1]]$social[i], scale = T)
#     # C[i, 2, t] = scale(guangzhou$covar[[t+1]]$happy[i], scale = T)
#     # C[i, 3, t] = scale(guangzhou$covar[[t+1]]$friends_hours[i], scale = T)
#     C[i, 1, t] = guangzhou$covar[[t+1]]$social[i]
#     C[i, 2, t] = guangzhou$covar[[t+1]]$happy[i]
#     C[i, 3, t] = guangzhou$covar[[t+1]]$pressure[i]
#   }
# }
# Missing <- list() #Enter as list manually, or if NAs are in data run the following:
# mindex = which(is.na(Y),arr.ind=TRUE)
# for(tt in 1:dim(Y)[3]){
#   Missing[[tt]] = unique(mindex[which(mindex[,3]==tt),1])
# }


require("ROCR")
predY <- array(0,dim=dim(Y))
for(tt in 1:TT){
  for(i in 1:n){
    for(j in c(1:n)[-i]){
      dijt <- sqrt(t(EX[,tt,i]-EX[,tt,j])%*%(EX[,tt,i]-EX[,tt,j]))
      expon <-  EB1[1]*(C[i, 1, tt]-C[j,1,tt]) + EB1[2]*(C[i, 2, tt]-C[j,2,tt]) + EB1[3]*(C[i, 3, tt]-C[j,3,tt]) + EBin*(1-dijt/Ew[j])+EBout*(1-dijt/Ew[i])
      # expon <-  EBin*(1-dijt/Ew[j])+EBout*(1-dijt/Ew[i])
      predY[i,j,tt] <- 1/(1+exp(-expon))
    }
  }
}

# future timepoint
colMeans(Cfuture[,,1])
apply(C[,,3],2,mean)
apply(C[,,2],2,mean)
apply(C[,,1],2,mean)
apply(Y, 3, mean)
mean(Yfuture)

MissData = F
if(MissData){
  for(tt in 1:TT){
    ROCR_pred <- ROCR_Y <- numeric(0)
    ROCR_pred <- c(ROCR_pred,upper.tri(predY[-Missing[[tt]],-Missing[[tt]],tt]))
    ROCR_pred <- c(ROCR_pred,lower.tri(predY[-Missing[[tt]],-Missing[[tt]],tt]))
    ROCR_Y <- c(ROCR_Y,upper.tri(Y[-Missing[[tt]],-Missing[[tt]],tt]))
    ROCR_Y <- c(ROCR_Y,lower.tri(Y[-Missing[[tt]],-Missing[[tt]],tt]))
  }
  pred <- prediction( predictions=ROC_pred,labels=ROC_Y )
  performance( pred, "auc")
  AUCPMean <- slot(performance( pred, "auc"),"y.values")[[1]]
}else{
  ROC_pred <- ROC_Y <- numeric(0)
  for(tt in 1:TT){
    ROC_pred <- c(ROC_pred,c(predY[,,tt][upper.tri(predY[,,tt])]))
    ROC_pred <- c(ROC_pred,c(predY[,,tt][lower.tri(predY[,,tt])]))
    ROC_Y <- c(ROC_Y,c(Y[,,tt][upper.tri(Y[,,tt])]))
    ROC_Y <- c(ROC_Y,c(Y[,,tt][lower.tri(Y[,,tt])]))
  }
  pred <- prediction( predictions=ROC_pred,labels=ROC_Y )
  performance( pred, "auc")
  AUCPMean <- slot(performance( pred, "auc"),"y.values")[[1]]  
}
print(AUCPMean) #0.9777968 0.9777008
pr_curve <- performance(pred, "prec", "rec")
plot(pr_curve, colorize = TRUE)
precision <- pr_curve@y.values[[1]]
recall    <- pr_curve@x.values[[1]]
threshold <- pr_curve@alpha.values[[1]]
f1 <- 2 * precision * recall / (precision + recall)

# remove NA values
valid <- which(!is.na(f1))
best_idx <- valid[which.max(f1[valid])]
best_threshold <- threshold[best_idx]
best_precision <- precision[best_idx]
best_recall    <- recall[best_idx]

best_threshold
best_precision
best_recall
plot(pr_curve, colorize = TRUE)
abline(v=best_recall, col='#5a95a6')
abline(h=best_precision, col='#5a95a6')
text(0.85, 0.91, paste0('cutoff = ', round(best_threshold,3)), cex = 1)

set.seed(100)
predYbin <- ifelse(predY>runif(Y), 1,0) 
mean(predYbin == Y) #0.8907047  0.8904029

perf <- performance(pred, "tpr", "fpr")
cutoffs <- data.frame(cut=perf@alpha.values[[1]], fpr=perf@x.values[[1]], 
                      tpr=perf@y.values[[1]])
cutoffs <- cutoffs[order(cutoffs$tpr+cutoffs$fpr , decreasing=TRUE),]
cf = subset(cutoffs, fpr < 0.1)[1,]


pdf("figs/ROC.pdf", width = 6, height = 6)
par(mar=c(5,5,2,2))
plot(perf, lwd = 2, cex.lab = 1.5, cex.axis = 1,
     xlab = 'false positive rate', ylab = 'true positive rate')
text(0.32, 0.95, paste0('cutoff = ', round(cf[1],3)), cex = 1)
abline(v=cf[2], col='#5a95a6')
abline(h=cf[3], col='#5a95a6')
dev.off()

#----------------------------------------------------
predYbin <- ifelse(predY>0.3, 1,0)
mean(predYbin[,,3] == Y[,,3])

dim(EX)
euc <- function(x, y){
  sqrt(sum((x-y)^2))
}
moves = array(0, dim = c(47, 29))
for (i in 1:47){
  for (t in 1:29){
    move = euc(EX[,t,i], EX[, t+1, i])
    moves[i, t] = move
  }
}
plot(moves[1,], type = 'l')
for (i in 2:47){
  lines(moves[i,], col = 'black')
}
pdf("figs/moves.pdf", width = 10, height = 6)
par(mar=c(5,5,2,2))
boxplot(moves, outline= F, xlab = 'week', ylab = 'moving distance from week to week', col = '#5a95a6', cex.lab = 2, cex.axis = 1.5)
dev.off()

distssame <- array(0, dim = c(47, 47, 30))
for(i in 1:47){
  for (j in 1:47){
    for (t in 1:30){
      distssame[i, j, t] = euc(EX[,t,i], EX[, t, j])
    }
  }
}
distsdist <- apply(distssame, c(1,3), mean)
boxplot(distsdist)

library(RColorBrewer)
cols = brewer.pal(4, "Blues")
# Use the following line with RColorBrewer
pal = colorRampPalette(cols)
# Rank variable for colour assignment
order = findInterval(C[,3, t], sort(C[,3, t]))

for(t in 1:TT){
  pdf(file = paste0("figsnew/changex/pos",t,".pdf"), width = 5, height = 5)
  par(mar = c(5, 5, 2, 2))
  dat <- data.frame(x1 = EX[1,t, ], x2 = EX[2,t,])
  plot(dat$x1, dat$x2, cex = C[,2, t]* 0.7,  col = pal(length(C[,3, t]))[order], pch = 20, xlab = expression(X[1]), ylab = expression(X[2]), main = paste("t = ", t), cex.lab = 1, cex.axis = 1, cex.main = 1.5,xlim = c(-0.03, 0.03), ylim = c(-0.03, 0.03))
  dev.off()
}

# Ew
hist(Ew)
reaches = array(0, dim = c(47, 47, 30))
for(i in 1:n){
  for (j in 1:n){
    for (t in 1:TT){
      if (distssame[i,j,t] < min(Ew[i], Ew[j])){
        reaches[i,j,t] <- "more"
      }
      if (distssame[i,j,t] > max(Ew[i], Ew[j])){
        reaches[i,j,t] <- "less"
      }
    }
  }
}

### c dist
cdiff1 <- array(0, dim = c(47, 47, 30))
cdiff2 <- array(0, dim = c(47, 47, 30))
cdiff3 <- array(0, dim = c(47, 47, 30))
for (i in 1:n){
  for (j in 1:n){
    for (t in 1:TT){
      cdiff1[i,j,t] <- C[i,1,t] - C[j,1,t]
      cdiff2[i,j,t] <- C[i,2,t] - C[j,2,t]
      cdiff3[i,j,t] <- C[i,3,t] - C[j,3,t]
    }
  }
}

t = 1
plot(cdiff1, Y)
plot(cdiff2[,,t], Y[,,t])
plot(cdiff3[,,t], Y[,,t])
mean(cdiff1[Y == 0])
mean(cdiff1[Y == 1])
mean(cdiff2[Y == 0])
mean(cdiff2[Y == 1])
mean(cdiff3[Y == 0])
mean(cdiff3[Y == 1])

