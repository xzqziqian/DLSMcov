# 
# Extension of Sewell & Chen (2015) to include dynamic non-network covariate evolution. Based on the dynamic latent space model code from https://publish.illinois.edu/yuguo/software/ and https://sewell.lab.uiowa.edu/software
# 
# Reference: Sewell, D. K., & Chen, Y. (2015). Latent space models for dynamic networks. Journal of the american statistical association, 110(512), 1646-1657.

c.update1 <-  cxxfunction(
  signature(Xitm1="numeric",DIMS="integer",TUNEX="numeric",Yy="numeric",
            BETAIN="numeric",BETAOUT="numeric",TUNEBIO="numeric",
            WW="numeric",t2X="numeric",s2X="numeric",
            xiBIN="numeric",xiBOUT="numeric",nuBIN="numeric",
            nuBOUT="numeric",CAUCHY="integer",
            RNORMS="numeric",RNORMSBIO="numeric",
            RNORMSB1 = "numeric", TUNEB1 = "numeric", BETA1 = "numeric", XIB1 = "numeric", NUB1 = "numeric",
            Cc = "numeric", NCOV = "numeric"),
  body='
/*Dims is c(n,p,TT,K)
  Z nxT; X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  int Cauchy = Rcpp::as<int>(CAUCHY);
  
  Rcpp::NumericVector XOLD(Xitm1);
  arma::cube Xold(XOLD.begin(),dims[1],dims[2],dims[0]);
  arma::cube Xnew = arma::zeros(dims[1],dims[2],dims[0]);
  for(int i=0;i<dims(0);i++)
{
  Xnew.slice(i)=Xold.slice(i);
}
  arma::mat Xconc = arma::zeros(dims(0)*dims(2),dims(1)); 
  arma::colvec centerings = arma::zeros(dims(1),1);
  
  Rcpp::NumericVector rnormsVec(RNORMS);
  arma::cube rnorms(rnormsVec.begin(),dims(1),dims(2),dims(0));
  arma::colvec rnormsBIO = Rcpp::as<arma::colvec>(RNORMSBIO);
  
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);
  
  Rcpp::NumericVector CC(Cc);
  double ncov = as<double>(NCOV);
  arma::cube C(CC.begin(),dims[0],ncov,dims[2]); //using 3 covariates as a test n x ncov x TT
  
  
  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);
  double xiBin = Rcpp::as<double>(xiBIN);
  double xiBout = Rcpp::as<double>(xiBOUT);
  double nuBin = Rcpp::as<double>(nuBIN);
  double nuBout = Rcpp::as<double>(nuBOUT);
  double BinNew =0, BoutNew =0;
  double tunex = Rcpp::as<double>(TUNEX);
  double tuneBIO = Rcpp::as<double>(TUNEBIO);
  
    
  //Rcpp::NumericVector B1(Beta1);
  //Rcpp::NumericVector B1New(Beta1);
  arma::colvec B1 = Rcpp::as<arma::colvec>(BETA1);
  arma::colvec B1new = Rcpp::as<arma::colvec>(BETA1);
  
  //Rcpp::NumericVector xiB1(XIB1);
  //Rcpp::NumericVector xiB1(NUB1);
  //Rcpp::NumericVector tuneB1(TUNEB1);
  //Rcpp::NumericVector rnormsB1(RNORMSB1);
  
  arma::colvec xiB1 = Rcpp::as<arma::colvec>(XIB1);
  arma::colvec nuB1 = Rcpp::as<arma::colvec>(NUB1);
  arma::colvec tuneB1 = Rcpp::as<arma::colvec>(TUNEB1);
  arma::colvec rnormsB1 = Rcpp::as<arma::colvec>(RNORMSB1);
  

  double t2 = Rcpp::as<double>(t2X);
  double s2 = Rcpp::as<double>(s2X);
  
  arma::colvec ww = Rcpp::as<arma::colvec>(WW);
  
  double AccProb =0;
  double dz=0, dx=0, uu=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  arma::colvec AccRate = arma::zeros(dims(0)*dims(2)+ ncov + 3 ,1);
  
  
  /*-------------------- Latent Positions-------------------*/
  
  for(int tt=0;tt < dims[2]; tt++)
{
  for(int i=0;i<dims[0];i++)
{
  AccProb=0;
  if(Cauchy<0.5)
{
  /*---Normal Random Walk---*/
  //  Xnew.slice(i).col(tt) = Xold.slice(i).col(tt) + tunex*arma::randn(dims(1),1);
  Xnew.slice(i).col(tt) = Xold.slice(i).col(tt) + tunex*rnorms.slice(i).col(tt);
}else{
  /*---Cauchy Random Walk---*/
  for(int ell=0;ell<dims(1);ell++)
{
  uu = arma::randu();
  double PI = 3.14;
  Xnew.slice(i)(ell,tt) = Xold.slice(i)(ell,tt) + tunex*tan(PI*(uu-0.5));
}
}
  
  
  for(int j=0;j<dims[0];j++)
{
  if(j != i){
  dz = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(j).col(tt),2);
  dx = arma::norm(Xold.slice(i).col(tt)-Xnew.slice(j).col(tt),2);
  AccProb += (dx-dz)*(Y.slice(tt)(j,i)*(BIN/ww(i)+BOUT/ww(j))+
  Y.slice(tt)(i,j)*(BIN/ww(j)+BOUT/ww(i)));
  double b1tmp = 0;
  for(int k=0;k<ncov;k++){
     b1tmp += B1(k)*(C.slice(tt)(i,k)-C.slice(tt)(j,k));
  }
  AccProb += log(1+exp(b1tmp + BIN*(1-dx/ww(i))+BOUT*(1-dx/ww(j))));
  AccProb += log(1+exp(b1tmp + BIN*(1-dx/ww(j))+BOUT*(1-dx/ww(i))));
  AccProb -= log(1+exp(b1tmp + BIN*(1-dz/ww(i))+BOUT*(1-dz/ww(j))));
  AccProb -= log(1+exp(b1tmp + BIN*(1-dz/ww(j))+BOUT*(1-dz/ww(i))));
  }
}
  if(tt==0)
{
  insides = trans(Xnew.slice(i).col(tt))*(Xnew.slice(i).col(tt))/t2;
  AccProb += -0.5*insides(0,0);
  insides = trans(Xold.slice(i).col(tt))*(Xold.slice(i).col(tt))/t2;
  AccProb -= -0.5*insides(0,0);
}
  if(tt>0)
{
  muit = Xnew.slice(i).col(tt-1);
  insides = trans(Xnew.slice(i).col(tt)-muit)*(Xnew.slice(i).col(tt)-muit)/s2;
  AccProb += -0.5*insides(0,0);
  insides = trans(Xold.slice(i).col(tt)-muit)*(Xold.slice(i).col(tt)-muit)/s2;
  AccProb -= -0.5*insides(0,0);  
}
  if(tt <dims[2]-1)
{
  muit = Xnew.slice(i).col(tt);
  insides = trans(Xnew.slice(i).col(tt+1)-muit)*(Xnew.slice(i).col(tt+1)-muit)/s2;
  AccProb += -0.5*insides(0,0);
  muit = Xold.slice(i).col(tt);
  insides = trans(Xnew.slice(i).col(tt+1)-muit)*(Xnew.slice(i).col(tt+1)-muit)/s2;
  AccProb -= -0.5*insides(0,0);  
}
  
  uu= arma::randu();
  if(uu<exp(AccProb))
{
  AccRate(ncov + 3 + tt*dims(0)+i) = 1;
}else
{
  Xnew.slice(i).col(tt) = Xold.slice(i).col(tt);
}
  
}
}
  
  /*---Centering---*/
  for(int i=0;i<dims(0);i++)
{
  for(int tt=0;tt<dims(2);tt++)
{
  Xconc.row(i*dims(2)+tt) = trans(Xnew.slice(i).col(tt));
}
}
  for(int ell=0;ell<dims(1);ell++)
{
  centerings(ell) = sum(Xconc.col(ell))/(dims(0)*dims(2));
}
  for(int i=0;i<dims(0);i++)
{
  for(int tt=0;tt<dims(2);tt++)
{
  for(int ell=0;ell<dims(1);ell++)
{
  Xnew.slice(i)(ell,tt) = Xnew.slice(i)(ell,tt) - centerings(ell);
}
}
}

 
  


  /*-----Beta1----*/
  

  for(int k=0;k<ncov;k++){
   B1new(k) = B1(k) + tuneB1(k)*rnormsB1(k);
  }
  
  for(int k=0;k<ncov;k++){
    AccProb=0;
    for(int tt = 0; tt < dims(2); tt++){
     for(int i = 0; i < dims(0); i++){
       for (int j = 0; j < dims(0); j++){
         if(i != j){
            dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(j).col(tt),2); //dijt
            
              double b1tmp = 0;
              for(int q=0;q<ncov;q++){
                 b1tmp += B1(q)*(C.slice(tt)(i,q)-C.slice(tt)(j,q));
              }
              
              double b1tmpnew = 0;
              for(int q=0;q<ncov;q++){
                if (q == k){
                  b1tmpnew += B1new(q)*(C.slice(tt)(i,q)-C.slice(tt)(j,q));
                }else{
                  b1tmpnew += B1(q)*(C.slice(tt)(i,q)-C.slice(tt)(j,q));
                }
              }
              
              
            
              
            AccProb += Y.slice(tt)(i,j)*(B1new(k)-B1(k))*(C.slice(tt)(i,k)-C.slice(tt)(j,k)) + 
                       log(1 + exp(b1tmp + BIN*(1-dx/ww(j)) + BOUT*(1-dx/ww(i)))) -
                       log(1 + exp(b1tmpnew + BIN*(1-dx/ww(j))+BOUT*(1-dx/ww(i))));
         }
       }
     }
    }
     AccProb += -0.5*(B1new(k)-nuB1(k))*(B1new(k)-nuB1(k))/xiB1(k);
     AccProb -= -0.5*(B1(k)-nuB1(k))*(B1(k)-nuB1(k))/xiB1(k);
       uu = arma::randu();
      if(uu < exp(AccProb)){
        AccRate(2+k) = 1;
      }else{
        B1new(k) = B1(k);
      }
  }

  
  
  
  
 
  

  
  
  
  /*-------------------- BetaIn and BetaOut-------------------*/
  AccProb=0;
  if(Cauchy<0.5)
  {
  //  BinNew = BIN + tuneBIO*arma::randn();
      BinNew = BIN + tuneBIO*rnormsBIO(0);
  }else{
      uu = arma::randu();
      double PI = 3.14;
      BinNew = BIN + tuneBIO*tan(PI*(uu-0.5));
  }
  
  for(int tt=0;tt<dims(2);tt++){
    for(int i = 0; i < dims(0); i++){
      for(int j = 0; j < dims(0); j++){
        if(i != j){
          dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(j).col(tt),2);
             
              
              double b1tmpnew = 0;
              for(int k=0;k<ncov;k++){
                 b1tmpnew += B1new(k)*(C.slice(tt)(i,k)-C.slice(tt)(j,k));
              }
              
              
              
          AccProb += Y.slice(tt)(i,j)*(BinNew-BIN)*(1-dx/ww(j)) + 
          log(1+exp(b1tmpnew + BIN*(1-dx/ww(j))+BOUT*(1-dx/ww(i)))) -
          log(1+exp(b1tmpnew + BinNew*(1-dx/ww(j))+BOUT*(1-dx/ww(i))));
        }
      }
    }
  }
  
  AccProb += -0.5*(BinNew-nuBin)*(BinNew-nuBin)/xiBin;
  AccProb -= -0.5*(BIN-nuBin)*(BIN-nuBin)/xiBin;
  
  uu= arma::randu();
  if(uu<exp(AccProb)){
      AccRate(0) = 1;
  }else{
    BinNew = BIN;
  }
  
  AccProb=0;
  if(Cauchy<0.5){
  //  BoutNew = BOUT + tuneBIO*arma::randn();
      BoutNew = BOUT + tuneBIO*rnormsBIO(1);
  }else{
    uu = arma::randu();
    double PI = 3.14;
    BoutNew = BOUT + tuneBIO*tan(PI*(uu-0.5));
  }
  
  
  for(int tt=0;tt<dims(2);tt++){
    for(int i = 0; i < dims(0); i++){
      for(int j = 0; j < dims(0); j++){
        if(i != j){
          dx = arma::norm(Xnew.slice(i).col(tt)-Xnew.slice(j).col(tt),2);
          
          double b1tmpnew = 0;
              for(int k=0;k<ncov;k++){
                 b1tmpnew += B1new(k)*(C.slice(tt)(i,k)-C.slice(tt)(j,k));
              }
              
              
          
          AccProb += Y.slice(tt)(i,j)*(BoutNew-BOUT)*(1-dx/ww(i)) + 
          log(1+exp(b1tmpnew + BinNew*(1-dx/ww(j))+BOUT*(1-dx/ww(i)))) -
          log(1+exp(b1tmpnew + BinNew*(1-dx/ww(j))+BoutNew*(1-dx/ww(i))));
        }
      }
    }
  }
  
  AccProb += -0.5*(BoutNew-nuBout)*(BoutNew-nuBout)/xiBout;
  AccProb -= -0.5*(BOUT-nuBout)*(BOUT-nuBout)/xiBout;
  
  uu= arma::randu();
  if(uu<exp(AccProb)){
    AccRate(1) = 1;
  }else{
    BoutNew = BOUT;
  }
  
  
  
  
  
  return Rcpp::List::create(Xnew,BinNew,BoutNew, B1new, AccRate);
  
  ', plugin="RcppArmadillo")



c.t2s2Parms <-  cxxfunction(
  signature(DATA="numeric",DIMS="integer",THETAT="numeric",
            THETAS="numeric",PHIT="numeric",PHIS="numeric"),
  body='
// Dims is c(n,p,TT,K)
  Rcpp::IntegerVector dims(DIMS);
  
  Rcpp::NumericVector Xvec(DATA);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  double thetaT=Rcpp::as<double>(THETAT);
  double phiT=Rcpp::as<double>(PHIT);
  double thetaS=Rcpp::as<double>(THETAS);
  double phiS=Rcpp::as<double>(PHIS);
  double shapeT=0, scaleT=0, shapeS=0, scaleS=0;
  
  arma::mat insides = arma::zeros(1,1);
  
  //---------------------------------------
  
  shapeT = thetaT + 0.5*dims(0)*dims(1);
  scaleT = phiT;
  shapeS = thetaS + 0.5*dims(0)*dims(1)*(dims(2)-1);
  scaleS = phiS;
  for(int i =0;i<dims(0);i++){
    insides = 0.5*trans(X.slice(i).col(0))*(X.slice(i).col(0));
    scaleT += insides(0,0);
    for(int tt=1;tt<dims(2);tt++){
      insides = 0.5*trans(X.slice(i).col(tt)-X.slice(i).col(tt-1))*
      (X.slice(i).col(tt)-X.slice(i).col(tt-1));
      scaleS += insides(0,0);
    }
  }
  
  
  return Rcpp::List::create(shapeT,scaleT,shapeS,scaleS);
  
  ', plugin="RcppArmadillo")

c.phisigmazParms <- cxxfunction(
  signature(DATA="numeric",DIMS="integer",THETAZ="numeric",
            PHIZ="numeric", XIPHI ="numeric", NUPHI="numeric", SIGMAZ = "numeric", Cc = "numeric",
            PHI = 'numeric',
            XI0 = "numeric", NU0 = "numeric", 
            THETAZ0 = 'numeric', PHIZ0 = 'numeric',
            SZ0 = 'numeric', MU = 'numeric'),
  body='
  
  double xi0 = Rcpp::as<double>(XI0);
  double nu0 = Rcpp::as<double>(NU0);
  double phiZ0 = Rcpp::as<double>(PHIZ0);
  double sz0 = Rcpp::as<double>(SZ0);
  double mu = Rcpp::as<double>(MU);
  double thetaZ0 = Rcpp::as<double>(THETAZ0);
  
  // Dims is c(n,p,TT,K)
  Rcpp::IntegerVector dims(DIMS);
  Rcpp::NumericVector Xvec(DATA);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  double thetaZ=Rcpp::as<double>(THETAZ);
  double phiZ=Rcpp::as<double>(PHIZ);
  double xiphi=Rcpp::as<double>(XIPHI);
  double nuphi=Rcpp::as<double>(NUPHI);
  double shapeZ = 0, scaleZ = 0, muphi = 0, sphi = 0, smu = 0 , mu0 = 0, shapeZ0 = 0, scaleZ0 = 0;
  double sz = Rcpp::as<double>(SIGMAZ); 
  double phi = Rcpp::as<double>(PHI);
  
  Rcpp::NumericVector CC(Cc);
  arma::mat C(CC.begin(),dims[0],dims[2]);
  
  /*arma::mat insides = arma::zeros(1,1);*/
  double insides = 0;
  
  /*----- sigmaz0 -----*/
  shapeZ0 = thetaZ0 + 0.5*dims(0);
  scaleZ0 = phiZ0;
  for (int i = 0; i < dims(0); i++){
    insides = 0.5*(C(i,0)-mu)*(C(i,0)-mu);
    scaleZ0 += insides;
  }
  
  /* ----- mu0 ----- */
  smu = 1/(1/xi0 + dims(0)/sz0);
  mu0 = nu0/xi0;
  for(int i =0; i<dims(0); i++){
    insides = C(i,0)/sz0;
    mu0 += insides;
  }
  mu0 = mu0*smu;
  
  
  /*-----sigmaz----*/
  
    shapeZ = thetaZ + 0.5*dims(0)*(dims(2)-1);
    scaleZ = phiZ;
    for(int i =0; i<dims(0); i++){
      for(int tt=1; tt<dims(2); tt++){
        insides = 0.5*(C(i,tt)-phi*C(i,tt-1))*(C(i,tt)-phi*C(i,tt-1));
        scaleZ += insides;
     }
   }
  
  
  /*-----phi----*/
  sphi = 0;
  double tmp = 1/xiphi;
  for(int i =0; i<dims(0); i++){
    for(int tt=1; tt<dims(2); tt++){
      insides = C(i,tt-1)*C(i,tt-1);
      tmp += 1/sz*insides;
    }
  }
  sphi = 1/(tmp);
  
  muphi = 0;
  tmp = nuphi/xiphi;
  for(int i =0; i<dims(0); i++){
    for(int tt=1; tt<dims(2); tt++){
      insides = C(i, tt)*C(i, tt-1);
      tmp += 1/sz*insides;
    }
  }
  muphi = sphi*tmp;
  
 
  
  return Rcpp::List::create(shapeZ,scaleZ,muphi,sphi, shapeZ0, scaleZ0, mu0, smu);
  
    ', plugin="RcppArmadillo")

c.t2s2Parms <-  cxxfunction(
  signature(DATA="numeric",DIMS="integer",THETAT="numeric",
            THETAS="numeric",PHIT="numeric",PHIS="numeric"),
  body='
// Dims is c(n,p,TT,K)
  Rcpp::IntegerVector dims(DIMS);
  
  Rcpp::NumericVector Xvec(DATA);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  double thetaT=Rcpp::as<double>(THETAT);
  double phiT=Rcpp::as<double>(PHIT);
  double thetaS=Rcpp::as<double>(THETAS);
  double phiS=Rcpp::as<double>(PHIS);
  double shapeT=0, scaleT=0, shapeS=0, scaleS=0;
  
  arma::mat insides = arma::zeros(1,1);
  
  //---------------------------------------
  
  shapeT = thetaT + 0.5*dims(0)*dims(1);
  scaleT = phiT;
  shapeS = thetaS + 0.5*dims(0)*dims(1)*(dims(2)-1);
  scaleS = phiS;
  for(int i =0;i<dims(0);i++)
{
  insides = 0.5*trans(X.slice(i).col(0))*(X.slice(i).col(0));
  scaleT += insides(0,0);
  for(int tt=1;tt<dims(2);tt++)
{
  insides = 0.5*trans(X.slice(i).col(tt)-X.slice(i).col(tt-1))*
  (X.slice(i).col(tt)-X.slice(i).col(tt-1));
  scaleS += insides(0,0);
}
}
  
  
  return Rcpp::List::create(shapeT,scaleT,shapeS,scaleS);
  
  ', plugin="RcppArmadillo")








c.initialize1 <-  cxxfunction(
  signature(Data="numeric",DIMS="integer",Yy="numeric",
            XSCALE="numeric",BETAIN="numeric",
            BETAOUT="numeric",WW="numeric"),
  body='
/*Dims is c(n,p,TT)
  Z nxT; X pxTxn; Y pxpxT;
  */
  Rcpp::IntegerVector dims(DIMS);

  Rcpp::NumericVector Xvec(Data);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);

  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);

  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);

  arma::colvec ww = Rcpp::as<arma::colvec>(WW);

  double Xscale= Rcpp::as<double>(XSCALE);

  double ret =0,dx=0, eta=0;

  /*---------------------------------------*/

  for(int tt=0;tt<dims(2);tt++)
{
  for(int i = 0; i < dims(0); i++)
{
  for(int j = 0; j < dims(0); j++)
{
  if(i != j)
{
  dx = Xscale*arma::norm(X.slice(i).col(tt)-X.slice(j).col(tt),2);
  eta = (BIN*(1-dx/ww(j))+BOUT*(1-dx/ww(i)));
  ret += Y.slice(tt)(i,j)*eta-
              log(1+exp(eta));
}
}
}
}

  return wrap(ret);

  ', plugin="RcppArmadillo")



c.initialize1.grad <-  cxxfunction(
  signature(Data="numeric",DIMS="integer",Yy="numeric",
            XSCALE="numeric",BETAIN="numeric",
            BETAOUT="numeric",WW="numeric"),
  body='
/*Dims is c(n,p,TT)
  Z nxT; X pxTxn; Y pxpxT;
  */
  Rcpp::IntegerVector dims(DIMS);

  Rcpp::NumericVector Xvec(Data);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);

  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);

  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);

  arma::colvec ww = Rcpp::as<arma::colvec>(WW);

  double Xscale= Rcpp::as<double>(XSCALE);

  double dx=0, eta=0;
  Rcpp::NumericVector ret(3);

  /*---------------------------------------*/

  for(int tt=0;tt<dims(2);tt++)
{
  for(int i = 0; i < dims(0); i++)
{
  for(int j = 0; j < dims(0); j++)
{
  if(i != j)
{
  dx = arma::norm(X.slice(i).col(tt)-X.slice(j).col(tt),2);
  eta = BIN*(1-Xscale*dx/ww(j))+BOUT*(1-Xscale*dx/ww(i));
  ret(0) = ret(0)+ dx*(BIN/ww(j)+BOUT/ww(i))*(1/(1+exp(-eta))-Y(i,j,tt));
  ret(1) = ret(1)+ (1-Xscale*dx/ww(j))*(Y(i,j,tt)-1/(1+exp(-eta)));
  ret(2) = ret(2)+ (1-Xscale*dx/ww(i))*(Y(i,j,tt)-1/(1+exp(-eta)));
}
}
}
}

  return wrap(ret);

  ', plugin="RcppArmadillo")


c.WAccProb1 <-  cxxfunction(
  signature(Data="numeric",DIMS="integer",Yy="numeric",
            BETAIN="numeric",BETAOUT="numeric",TUNEW="numeric",
            WWOld="numeric",WWNew="numeric",  Cc = "numeric", BETA1 = "numeric", NCOV = "integer"),
  body='
/*Dims is c(n,p,TT,K)
  Z nxT; X pxTxn; Y pxpxT; 
  */
  Rcpp::IntegerVector dims(DIMS);
  
  Rcpp::NumericVector Xvec(Data);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);
  
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);
  
  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);
  arma::colvec B1 = Rcpp::as<arma::colvec>(BETA1);
  
  arma::colvec wwOld = Rcpp::as<arma::colvec>(WWOld);
  arma::colvec wwNew = Rcpp::as<arma::colvec>(WWNew);
  double tuneW = Rcpp::as<double>(TUNEW);
  
  double AccProb =0,dx=0, uu=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);
  int AccRate = 0;
  
  Rcpp::NumericVector CC(Cc);
  double ncov = as<double>(NCOV);
  arma::cube C(CC.begin(),dims[0],ncov,dims[2]); //using 3 covariates as a test n x ncov x TT
  
  
  /*---------------------------------------*/
  
  for(int tt=0;tt<dims(2);tt++)
{
  for(int i = 0; i < dims(0); i++)
{
  for(int j = 0; j < dims(0); j++)
{
  if(i != j)
{
  dx = arma::norm(X.slice(i).col(tt)-X.slice(j).col(tt),2);
  
   double b1tmp = 0;
  for(int k=0;k<ncov;k++){
     b1tmp += B1(k)*(C.slice(tt)(i,k)-C.slice(tt)(j,k));
  }
  
  
  AccProb += Y.slice(tt)(i,j)*dx*(BIN*(1/wwOld(j)-1/wwNew(j)) + 
  BOUT*(1/wwOld(i)-1/wwNew(i))) +
  log(1+exp(b1tmp + BIN*(1-dx/wwOld(j))+BOUT*(1-dx/wwOld(i)))) -
  log(1+exp(b1tmp + BIN*(1-dx/wwNew(j))+BOUT*(1-dx/wwNew(i))));
}
}
}
}
  
  for(int i =0;i<dims(0);i++)
{
  AccProb += (tuneW*wwNew(i)-1)*log(wwOld(i))- (tuneW*wwOld(i)-1)*log(wwNew(i)) -
  (tuneW*wwNew(i)-0.5)*log(tuneW*wwNew(i))-tuneW*wwNew(i) +
  (tuneW*wwOld(i)-0.5)*log(tuneW*wwOld(i))+tuneW*wwOld(i);
}
  
  uu= arma::randu();
  if(uu<exp(AccProb))
{
  AccRate = 1;
}else
{
  wwNew = wwOld;
}
   
  return Rcpp::List::create(wwNew,AccRate);
  
  ', plugin="RcppArmadillo")

c.missingcov <- cxxfunction(
  signature(DIMS="integer", Cc = "numeric", PHI = 'numeric', SIGMAZ = 'numeric', NCOV = 'integer', MU0 = 'numeric', S0 = 'numeric', RN = 'numeric'),
  body='
  
  Rcpp::IntegerVector dims(DIMS);
  
  /* cov */
   Rcpp::NumericVector CC(Cc);
  double ncov = as<double>(NCOV);
  arma::cube C(CC.begin(),dims[0],ncov,dims[2]); //using 3 covariates as a test n x ncov x TT
  
  Rcpp::NumericVector rn1(RN);
  arma::cube rn(rn1.begin(),dims[0],ncov,dims[2]);
  

  arma::colvec phi = Rcpp::as<arma::colvec>(PHI);
  arma::colvec sigmaZ = Rcpp::as<arma::colvec>(SIGMAZ);
  arma::colvec mu0 = Rcpp::as<arma::colvec>(MU0);
  arma::colvec s0 = Rcpp::as<arma::colvec>(S0);
  
  
 for(int i = 0; i < dims(0); i++){
  for (int k = 0; k< ncov; k++){
    for (int ttt = 0; ttt < dims(2); ttt++){
      if (ttt > 0){
        C.slice(ttt)(i,k) = phi(k)*C.slice(ttt-1)(i,k) + sigmaZ(k)*rn.slice(ttt)(i,k);
      }else{
        C.slice(ttt)(i,k) = mu0(k) + s0(k)*rn.slice(ttt)(i,k);
      }
    }
  }
 }

  
    return Rcpp::wrap(C);

  ', plugin="RcppArmadillo")



c.missing <-  cxxfunction(
  signature(Data="numeric",DIMS="integer",MMM="integer",Yy="numeric",Ttt="integer",
            BETAIN="numeric",BETAOUT="numeric",WW="numeric", BETA1 = "numeric", 
            Cc = "numeric", NCOV = "integer"),
  body='
/*Dims is c(n,p,TT,K);
  Z nxT; X pxTxn; Y pxpxT;
  */
  Rcpp::IntegerVector dims(DIMS);
  int ttt = Rcpp::as<int>(Ttt)-1;

  Rcpp::NumericVector Xvec(Data);
  arma::cube X(Xvec.begin(),dims[1],dims[2],dims[0]);

  Rcpp::IntegerVector MM(MMM);
  Rcpp::NumericVector YY(Yy);
  arma::cube Y(YY.begin(),dims[0],dims[0],dims[2]);

  double BIN = Rcpp::as<double>(BETAIN);
  double BOUT = Rcpp::as<double>(BETAOUT);

  arma::colvec ww = Rcpp::as<arma::colvec>(WW);

  double dx=0, uu=0, Prob=0;
  arma::mat insides = arma::zeros(1,1);
  arma::colvec muit = arma::zeros(dims[1],1);

   Rcpp::NumericVector CC(Cc);
  double ncov = as<double>(NCOV);
  arma::cube C(CC.begin(),dims[0],ncov,dims[2]); //using 3 covariates as a test n x ncov x TT
  
  arma::colvec B1 = Rcpp::as<arma::colvec>(BETA1);
  

  
  /*---------------------------------------*/

  for(int i = 0; i < dims(0); i++)
{
  if(std::find(MM.begin(),MM.end(),i) !=MM.end())
{
  for(int j = 0; j < dims(0); j++)
{
  if(i != j)
{
  dx = arma::norm(X.slice(i).col(ttt)-X.slice(j).col(ttt),2);
  Prob = BIN*(1-dx/ww(j))+BOUT*(1-dx/ww(i));
  Prob = 1/(1+exp(-Prob));
  uu= arma::randu();
  if(uu<Prob)
{
  Y(i,j,ttt) = 1;
}else{
  Y(i,j,ttt) = 0;
}
  /*  if(std::find(MM.begin(),MM.end(),j) ==MM.end())
{
  double b1tmp = 0;
  for(int k=0;k<ncov;k++){
     b1tmp += B1(k)*(C.slice(tt)(i,k)-C.slice(tt)(j,k));
  }
  Prob = b1tmp + BIN*(1-dx/ww(i))+BOUT*(1-dx/ww(j));
  Prob = 1/(1+exp(-Prob));
  uu= arma::randu();
  if(uu<Prob)
{
  Y(j,i,ttt) = 1;
}else{
  Y(j,i,ttt) = 0;
}
}*/

}
}
}
}

  return Rcpp::wrap(Y);

  ', plugin="RcppArmadillo")








