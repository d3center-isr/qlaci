

########################################################################
qlaci <- function( H10, H11, A1, Y1, H20, H21, A2, Y2,S,
                   c1=NULL, c2=NULL, alpha=0.05, gridscale=2.58, ngrid=10, lambda=log(log(NROW(S))), nb=1000){
  # qlaci performs two-stage q learning as well as computing adative confidence interval for stage 1/2 coefficients.
  # S: c(n,), 
  # H10: c(n, p10), H11: c(n, p11), A1: c(n,), Y1: c(n,), 
  # H20: c(n, p20), H21: c(n, p21), A2: c(n,), Y2:(n,), 
  # A1 and A2 are 1/-1 coded. 
  # c1: c((p10+p11),nc1),c2: c((p20+p21),nc2) specifiy contrasts for stage 1 and 2 coefficients
  # each contrast is a column. 
  # if c1 and/or c2 is NULL, no confidence interval is computed. 
  
  #output: point estimate beta1, beta2 for \beta_1 and \beta_2
  #        confidence interval for c*\beta_1, c * \beta_2
  
  #todo: input check and conversion; 
  
  if (any(c(NROW(S), NROW(H10), NROW(H11), NROW(H20), NROW(H21), NROW(A1),
           NROW(A2), NROW(Y1), NROW(Y2))!=NROW(S))) {
    stop("Input matrices and vectors should have the same number of rows");
  }
  
  if (mode(S)!="logical" && !all(unique(S) %in% c(0,1))) {
    stop("S should be either logical or 0-1 valued");
  }
  
  
  ok<- complete.cases(H10, H11, A1, Y1,Y2,  S);
  if (!all(ok)) {
    stop("Cases should be complete, i.e. have no missing values"); 
  }
  
  ok<- complete.cases(  H20[S==1,], H21[S==1,], A2[S==1], Y2[S==1]);
  if (!all(ok)) {
    stop("Cases should be complete, i.e. have no missing values"); 
  }

  
  if (!all(unique(A1) %in% c(-1,1)) || !all(unique(A2[S==1]) %in% c(-1,1))) {
    stop("A1 and A2 must be -1/+1 coded");  
  }
  
  
  

  H10<-as.matrix(H10);
  H11<-as.matrix(H11);
  H20<-as.matrix(H20);
  H21<-as.matrix(H21);
  
  
  if (!is.null(c1)){
    c1<-as.matrix(c1);
  }
  
  if (!is.null(c2)){
    c2<-as.matrix(c2);
  }
  
  
  A1<-as.vector(A1); A2<-as.vector(A2);
  Y1<-as.vector(Y1); Y2<-as.vector(Y2);
  
  S<- as.logical(S); 

  # a few checks that try to catch some common input errors:

  n<- NROW(S);
  
  p10<- NCOL(H10);
  p11<- NCOL(H11);
  p20<- NCOL(H20);
  p21<- NCOL(H21);
  n2 <- sum(S);
  
  
  #call QL function to get point estimate for \beta_1 and \beta_2
  #and also the estimate for value of action Y
  
  #browser();
  ret<-QL(H10, H11, A1, Y1, H20, H21, A2, Y2, S); 
  

  
  #todo: check return code, continue only if it is an OK flag
  
  ci1=NULL;
  ci2=NULL;
  
  if (!is.null(c2)) {
    ci2=compute_ci2(c2, ret$beta2h, H10, H11, A1, Y1, H20, H21, A2, Y2, S, nb, alpha);  
  }

  if (is.null(c1)){
    return( list(stg1coeff=ret$beta1h, stg2coeff=ret$beta2h, ci1=ci1, ci2=ci2) );
  }
    
  
  
  beta1h<- ret$beta1h; 
  beta2h<- ret$beta2h; 
  Yh <- ret$Yh;
  beta10h <- beta1h[1:p10];
  beta11h <- beta1h[(p10+1):(p10+p11)];
  beta20h <- beta2h[1:p20];
  beta21h <- beta2h[(p20+1):(p20+p21)];
  
  #compute \Sigma_{2,n}^{(2,2)}
  Sigma2n <- cov_beta_2(H10, H11, A1, Y1, H20, H21, A2, Y2, S, nb)
  Sigma2n22 <- Sigma2n[(p20+1):(p20+p21), (p20+1):(p20+p21), drop=F];

   
  #hack:
  #precompute H21*Sigma2n22* H21';
  hsh= vector('numeric', length=n);
  for (i in 1:n) {
	  hsh[i]= H21[i,] %*% Sigma2n22 %*% H21[i,] ;  
  }
  
  
  #generate grid \Tao= 
  center<- sqrt(n)* beta21h;  
  r<- gridscale * sqrt(max(n2* diag(as.matrix(Sigma2n22)))); 
  Gamma<- compute_Gamma(center, r,  ngrid);
  
  #boot package is a little too much, diy
  
  BS <- generate_bootstrap(S, nb);
  
  #z<- matrix(0, nrow=NCOL(BS),  ncol=NROW(Gamma));
  
  z<- array(dim=c(dim(c1)[2], NCOL(BS), NROW(Gamma) )); 
  
  for (i in 1:NCOL(BS)) {
#     datab <- get_bootstrap(BS[,i], H10, H11, A1, Y1, 
#               H20, H21, A2, Y2, S);
#     H10b<-datab$H10b;
#     H11b<-datab$H11b;
#     A1b<- datab$A1b;
#     Y1b<- datab$Y1b;
#     H20b<- datab$H20b;
#     H21b<- datab$H21b;
#     A2b <- datab$A2b;
#     Y2b <- datab$Y2b;
#     Sb<-datab$Sb;
        
    bs<- BS[,i];
    H10b<-H10[bs,,drop=F];
    H11b<-H11[bs,,drop=F];
    A1b<- A1[bs];
    Y1b<- Y1[bs];
    H20b<- H20[bs,,drop=F];
    H21b<- H21[bs,,drop=F];
    A2b <- A2[bs];
    Y2b <- Y2[bs];
    Sb<-S[bs];
    
    n2b<- sum(Sb);
    
    B1b<- cbind(H10b, A1b* H11b);
    B2b<- cbind(H20b, A2b* H21b);
    
    Sigma1nb <- compute_Sigma1n(B1b);
	
    retb<- QL(H10b, H11b, A1b, Y1b, H20b, H21b, A2b, Y2b, Sb);
    beta1b<- retb$beta1h;
    beta2b<- retb$beta2h;
    Yb <- ret$Yh;
    beta10b <- beta1b[1:p10];
    beta11b <- beta1b[(p10+1):(p10+p11)];
    beta20b <- beta2b[1:p20];
    beta21b <- beta2b[(p20+1):(p20+p21)];
    
    #function(beta21, H21,S, lambda, T)
    hshb<- hsh[BS[,i]];
    I1b <- compute_I1(beta21b, H21b,Sb, lambda,hshb);
    I2b <- !I1b;
    
    Yh<- Sb* ( H20b %*% beta20h) + Sb * abs(H21b %*% beta21h) + (1-Sb) * Y2b + Y1b;
    
    Sigma1nIb = solve(Sigma1nb);
    
    
    
    #for g
    n <- length(Sb);
    H21bS <- H21b[Sb, ,drop=F];      #as.matrix(H21b)[Sb, ,drop=F];
    B1bS<- B1b[Sb, , drop=F];
    T0<-t(c1) %*% Sigma1nIb %*% t(B1bS);
    
    #for z
    n<- length(Sb);
    T00<- sqrt(n) * t(c1)  %*% Sigma1nIb;
    beta20h<- beta2h[1:p20];
    beta21h<- beta2h[(p20+1): (p20+p21)];
    
    beta20b <- beta2b[1:p20];
    beta21b <- beta2b[(p20+1): (p20+p21)];
    B1bS <- B1b[Sb, , drop=F];
    H20bS <- H20b[Sb, , drop=F];
    H21bS <- H21b[Sb, , drop=F]; #as.matrix(H21b)[Sb, , drop=F];
    T000=T00 %*% (t(B1b) %*% (Yh-B1b %*% beta1h) + t(B1bS) %*% H20bS %*% (beta20b-beta20h) +
      t(B1bS) %*% ( (abs(H21bS %*%  beta21b) -abs(H21bS %*%  beta21h))*I1b))/n;
    
    
    for (j in 1:NROW(Gamma)){
      gamma=Gamma[j,];
      
      #g <- compute_g(gamma, Sigma1nIb, B1b, beta21b, beta21h, 
      #               I2b, Sb, H21b, c1); 
      
      
      
      g<- T0 %*% 
        ((abs(H21bS %*% (sqrt(n)* (beta21b-beta21h) + gamma))-abs(H21bS %*% gamma))*I2b)  /n; 
      
      
      
      #z[i,j]<- compute_z(Sigma1nIb, B1b, beta1b, beta2b, Yh, I1b, 
      #                   g, beta1h, beta2h, Sb, H10b, H11b, H20b, H21b, c1);
      
      z[,i,j] = T000 +g;
      
      
      
    } #inner for j
    #cat(i);
  } #outer for i
  
  
  zmax<- apply(z,c(1,2),max);
  zmin<- apply(z,c(1,2), min);
  
  u<- apply(zmax, 1, function(row){quantile(row, 1-alpha/2)});
  l<- apply(zmin, 1, function(row){quantile(row, alpha/2)});
  
  Estimates=t(c1)%*%beta1h
  ci1 <- cbind(Estimates,t(c1)%*%beta1h-u/sqrt(n), t(c1)%*%beta1h -l/sqrt(n) ); 

 colnames(ci1)<-c("est","low","upp")
 


  
  return( list(stg1coeff=beta1h, stg2coeff=beta2h, ci1=ci1, ci2=ci2) );
  
}
########################################################################

### alg 2
########################################################################
QL<- function(H10, H11, A1, Y1, H20, H21, A2, Y2, S){

  B2<- cbind(H20, A2*H21); 
  
  Y2S<- Y2[S]; B2S<- B2[S, , drop=F];
  
  ls2 <- lsfit(B2S, Y2S, intercept=FALSE); # intercept should already be included in B2S; 
  
  beta2h<- ls2$coef; 
  
  p20<- NCOL(H20); p21<-NCOL(H21);
  
  beta20h<- beta2h[1:p20];
  beta21h<- beta2h[(p20+1):length(beta2h)];
  
  B1<- cbind(H10, A1*H11); 

#   tryCatch( { Yh<- S* ( as.matrix(H20) %*% beta20h) + S * abs(as.matrix(H21) %*% beta21h) + (1-S) * Y2 + Y1
# }, error=function(e){print(str(H20)); print(str(beta20h))})
  
  Yh<- S* (H20 %*% beta20h) + S * abs(H21 %*% beta21h) + (1-S) * Y2 + Y1;
  
  ls1<- lsfit(B1, Yh,intercept=FALSE); 
  
  
  beta1h <- ls1$coef; 
  
  #check conditional number of ls1 and ls2
  if (kappa(ls2$qr)>70){
    warning("The second stage regression may have multicollinearity/singularity");
  } 
  
  if (kappa(ls1$qr)>70){
    warning("The first stage regression may have multicollinearity/singularity");
  }
  
  return(list(beta1h=beta1h, beta2h=beta2h, Yh=Yh, ls1=ls1, ls2=ls2))
  
}
########################################################################

# alg 3:
########################################################################
generate_bootstrap<-function (S, nb=1000, minS=20) {
  n<-NROW(S);  
  R<-matrix(0, nrow=n, ncol=nb);
  for (i in 1:nb) {
    r <- sample.int(n,n, replace=TRUE);
    if (sum(S[r])<minS) {
      stop('stage 2 sample is too small for valid bootstrap');
    }
    R[,i]<-r
  }
  
  return(R)
}
########################################################################

#alg 4:
########################################################################
# hsbn is for speed-up the calculation of H21S %*% Sigma2n22 % H21S
compute_I1<- function(beta21, H21,S, lambda, hsbh){
  #n2<- sum(S);
  H21S <- H21[S, ,drop=F];
  
  T1 <- (H21S %*% beta21)^2;
  
  T2 <- hsbh[S];
  
  I1<- T1 > (lambda * T2); 
  
}
########################################################################


#alg 7:
########################################################################
compute_Sigma1n<- function(B1){
  n<- NROW(B1);
  Sigma1n <- t(B1) %*% B1 /n; 
  
  return(Sigma1n)
}
########################################################################

#alg 8:
########################################################################
cov_beta_2<- function(H10, H11, A1, Y1, H20, H21, A2,Y2, S, nb=1000){
  BS<- generate_bootstrap(S,nb);
  Beta2b <- matrix(NA, nb, NCOL(H20)+ NCOL(H21));
  
  for (i in 1:nb) {
    bs<- BS[,i];
    H10b<-H10[bs, , drop=F];
    H11b<-H11[bs, , drop=F];
    A1b<- A1[bs];
    Y1b<- Y1[bs];
    H20b<- H20[bs, , drop=F];
    H21b<- H21[bs, , drop=F];
    A2b <- A2[bs];
    Y2b <- Y2[bs];
    Sb<-S[bs];
    
    ret<-QL(H10b, H11b, A1b, Y1b, H20b, H21b, A2b, Y2b, Sb); 
    beta2b<-ret$beta2h; 
    
    Beta2b[i,]<- beta2b; 
    
  }

  Sigma2n <- cov(Beta2b);
  
  return(Sigma2n)
}
########################################################################

#alg 9:
########################################################################
compute_Gamma<- function (c, r, ngrid=10){
  dimc <- length(c);
  clow <- c-r; 
  cup <- c+r; 
  
  clowup = list(); 
  
  for (i in 1:dimc){
    clowup[[i]]<- seq(clow[i], cup[i], length.out=ngrid);
  }
  
  Gamma<- expand.grid(clowup); 
  
  Gamma<- as.matrix(Gamma); 
  
  return(Gamma);
}
########################################################################

#alg 10:
compute_ci2<- function(c2, beta2, H10, H11, A1, Y1, H20, H21, A2, Y2, S, nb=1000, alpha=0.05){
  BS <- generate_bootstrap(S, nb);
  
  #z<- matrix(0, nrow=NCOL(BS),  ncol=NROW(Gamma));
  
  z<- array(dim=c(dim(c2)[2], NCOL(BS))); 
  
  for (i in 1:NCOL(BS)) {
    
    bs<- BS[,i];
    H10b<-H10[bs,,drop=F];
    H11b<-H11[bs,,drop=F];
    A1b<- A1[bs];
    Y1b<- Y1[bs];
    H20b<- H20[bs,,drop=F];
    H21b<- H21[bs,,drop=F];
    A2b <- A2[bs];
    Y2b <- Y2[bs];
    Sb<-S[bs];
    n <- length(Sb);
    
    retb<- QL(H10b, H11b, A1b, Y1b, H20b, H21b, A2b, Y2b, Sb);  
    
    beta2b= retb$beta2;
    z[,i]= t(c2) %*% (beta2b-beta2) *sqrt(n);
  }
  
  u<- apply(z , 1, function(row){quantile(row, 1-alpha/2)});
  l<- apply(z, 1, function(row){quantile(row, alpha/2)});
  
  c2u = t(c2) %*% beta2 -u/sqrt(n); 
  c2l = t(c2) %*% beta2 - l/sqrt(n);
 Estimates2=t(c2)%*%beta2 
 outc2<-cbind(est=Estimates2,low=c2u, upp=c2l)
 colnames(outc2)<-c("est","low","upp")
 
 
#  return(cbind(est=Estimates2,low=c2u, upp=c2l));
  return(outc2);
}
