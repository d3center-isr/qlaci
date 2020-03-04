 % File src/library/base/man/load.Rd
\name{qlaci}
\alias{qlaci} 
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Q-Learning: A Data Analysis Method for Constructing Adaptive Interventions 
}
\description{
The \code{"qlaci"} package  implements a generalization of Q-learning, a method developed in computer science, to inform the development of adaptive interventions. This package uses the Adaptive Confidence Interval (ACI) method developed by Laber et al. [2010] to construct the confidence intervals.


}
\usage{
qlaci(H10, H11, A1, Y1, H20, H21, A2, Y2, S, c1, c2, alpha, gridscale,
 ngrid, nb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{H10}{
A matrix where each column represents a main effect term involved in the first stage regression (usually includes a column of ones as the intercept).
}
  \item{H11}{
A matrix where each column represents a term interacting with the first stage intervention A1 (usually includes a column of ones as the main effect of A1).
}
  \item{A1}{
A vector representing the first stage intervention (binary, coded as -1 and +1).}
  \item{Y1}{
Intermediate Outcome (continuous).
}
  \item{H20}{
A matrix where each column represents a main effect term involved in the second stage regression (usually includes a column of ones as the intercept).
}
  \item{H21}{
A matrix where each column represents a term interacting with the second stage intervention A2 (usually includes a column of ones as the main effect of A2).
}
  \item{A2}{
A vector representing the second stage intervention (binary, coded as -1 and +1).
}
  \item{Y2}{
Final Outcome (continuous).
}
  \item{S}{
Rerandomization status at stage 2 (binary, coded 0 and 1).
}
  \item{c1}{
The contrast matrix that is used for the stage 1 regression (default = NULL).
}
  \item{c2}{
The contrast matrix that is used for the stage 2 regression (default = NULL).
}
  \item{alpha}{
Significance value of the confidence intervals (default = 0.05).
}
  \item{gridscale}{
Scale for grid range (default = 2.58).

}
  \item{ngrid}{
Number of grid points (default = 10).
}
  \item{nb}{
Number of bootstraps to be used in stage 1 & 2 (default = 1000).

}
}
\details{
Q-learning fits a linear regression model at each decision stage. The \code{"qlaci"} package accomodates designs with two stages.  The regressions are fit sequentially starting from the stage 2, working backward to the stage 1. At each stage, the purpose of the regression is to evaluate usefulness of candidate tailoring variables and to find the treatment option at that stage that maximizes the outcome. 
To fit the stage 2 model, only data from participants who are re-randomized is used. The primary outcome is regressed on the intermediate outcomes, assigned treatments at stages 1 & 2 and baseline covariates to fit the stage 2 model. Using the estimated parameters in the stage 2 model, the stage 2 treatment option that optimizes the expected value of the primary outcome can be derived. Next, data from all participants is used to fit the first-stage model; this model is fit by regressing the appropriate outcome on the baseline covariates and the treatment assigned at the stage 1. 

Missing covariates and/or missing primary outcomes are not permitted. See the user's guide for more detail.
}
\value{
%%  ~Describe the value returned
An object 'qlaci' is a list of containing the following components:
  \item{stg1coeff }{The stage 1 regression coefficients.}
  \item{stg2coeff }{The stage 2 regression coefficients.}
  \item{ci1 }{Estimated  \code{"contrast"} of the stage 1 regression parameters.}
  \item{ci2 }{Estimated  \code{"contrast"} of the stage 2 regression parameters.}
%% ...
}
\references{

A. Ertefaie, D. Almirall, L. Huang,  J. J. Dziak, A.T. Wagner, & S. A. Murphy (2012). SAS PROC QLEARN users' guide (Version 1.0.0). University Park: The Methodology Center, Penn State. Available from http://methodology.psu.edu.

E. Laber, M. Qian, D.J. Lizotte, and S.A. Murphy. Statistical inference in dynamic treatment regimes. Arxiv preprint arXiv:1006.5831, 2010.

I. Nahum-Shani, M. Qian, D. Almiral, W.. Pelham, B. Gnagy, G. Fabiano, J. Waxmonsky, J. Yu and S.A. Murphy (2010a). Experimental Design and Primary Data Analysis Methods for Comparing Adaptive Interventions. Technical Report 10-108, The Methodology Center, The Pennsylvania State University, University Park, PA.

I. Nahum-Shani, M. Qian, D. Almiral, W.. Pelham, B. Gnagy, G. Fabiano, J. Waxmonsky, J. Yu and S.A. Murphy (2010b). Q-Learning: a data analysis method for constructing adaptive interventions. Technical Report 10-107, The Methodology Center, The Pennsylvania State University, University Park, PA.
}
\author{
Kun Deng, Ashkan Ertefaie and Susan Murphy.

Maintainer: Kun Deng <kundeng@umich.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## load qlearning package
   library(qlaci)

## check the usage of the qlearning package
## this will help you understand the arguments to the function
   ?qlaci
##  (An html help page about the qlearning package will pop up.)

## Example 1: Autism simulated data:
## load the example data set provided in the qlearning package
   data(ASDdat)

## take a look at the first few observations in the data
   head(ASDdat)
 

 ## We are going to use the qlearning() function to fit the following two models:
    ## stage 1 model: E(Y | A1, O11) = b1*H10 + b2*A1*H11
    ## stage 2 model: E(Y | A1, O11, A2, O21) = b3*H20 + b4*A2*H21
    ## note:
    ## b1, b2, b3 and b4 are vectors. The stage 2 model is fitted among non-responders to     
    ## A1=-1.
    ## construct covariates used in the first-stage and the second-stage regression
  attach(ASDdat); 
## Since the confidence intervals are obtained using a bootstrap method, we use 
## the "set.seed" comment to reproduce the same intervals.
  set.seed(300); 
  H10<- cbind(1,O11);
  colnames(H10)<-c("int","O11");

  H11<- cbind(1, O11);
  colnames(H11)<-c("A1","A1O11");
  H20<- cbind(1,O11,O21);
  colnames(H20)<-c("int","O11","O21");
  H21<- cbind(1,O21); 
  colnames(H21)<-c("A2","A2O21");
  S<- as.logical(S); 


## Construct contrast matrices
row1<-c(1,2,1,2) # mean  under JAE+EMT for children w/O11 = 2 = low; 
row2<-c(1,2,-1,-2) #mean  under JAE+AAC for children w/O11 = 2 = low; 
row3<-c(0,0,2,4) #mean difference in  for JAE+EMT - JAE+AAC for children w/O11 = 2 = low; 
row4<-c(0,0,2,8) #mean difference in  for JAE+EMT - JAE+AAC for children w/O11 = 4 = hi; 
c1<-rbind(row1, row2, row3, row4) 
c2<-NULL       
    ## Run qlearning function to get estimates and confidence intervals for the contrasts
    
result1<-qlaci(H10, H11, A1, Y1, H20, H21, A2, Y2, S,c1=t(c1) ,nb=1000);
print(result1) 
    detach(ASDdat)     

## Example 2:
## load the example data set provided in the qlearning package
   data(dat2)

## take a look at the first few observations in the data
 head(dat2)
 

 ## We are going to use the qlearning() function to fit the following two models:
    ## stage 1 model: E(Y | A1, O11) = b1*H10 + b2*A1*H11
    ## stage 2 model: E(Y | A1, O11, A2, O21) = b3*H20 + R*A2*H21+(1-R)*A2*H21
    ## note:
    ## b1, b2, b3 and b4 are vectors and R is the responder indicator which is 1 if responder and     
    ##zero otherwise. Since in this example all the individuals are re-
    ##randomized at the stage 2, the stage 2 model is fitted among all the               
    ##indivduals.     
    ## construct covariates used in the first-stage and the second-stage regression
  attach(dat2);
  set.seed(300);
  H10<- cbind(1,O11);
  colnames(H10)<-c("int","O11");

  H11<- cbind(1, O11);
  colnames(H11)<-c("A1","A1O11");
  Y1<-rep(0,200);
  H20<- cbind(1,O11,O21); 
  colnames(H20)<-c("int","O11","O21");
  H21<- cbind(R,R*O21,1-R,(1-R)*O21);  # R=1=responders
  colnames(H21)<-c("A2R","A2RO21","A2(1-R)","A2(1-R)O21"); 
  Y2<- Y; 
  S<- rep(1,200); # everyone is randomized at stage 2
  
    ## Construct contrast matrices
    c1<-diag(4); # number of rows must be equal to the number of parameters in the stage 1 model.
    c2<-diag(7); # number of rows must be equal to the number of parameters in the stage 2 model.
  
    ## Run qlearning function to get estimates and confidence intervals for the contrasts
    
result2<-qlaci(H10, H11, A1, Y1, H20, H21, A2, Y2, S,c1=t(c1),c2=t(c2) ,nb=1000);
print(result2);
  detach(dat2);

## Example 3: A simulated data with no embedded tailoring variables:
## load the example data set provided in the qlearning package
   data(dat3)

## take a look at the first few observations in the data
 head(dat3)
 

 ## We are going to use the qlearning() function to fit the following two models:
    ## stage 1 model: E(Y | A1, O11) = b1*H10 + b2*A1*H11
    ## stage 2 model: E(Y | A1, O11, A2, O21, O22) = b3*H20 + b4*A2*H21
    ## note:
    ## b1, b2, b3 and b4 are vectors. Since in this example all the individuals are re-
    ##randomized at the stage 2, the stage 2 model is fitted among all the               
    ##individuals. 
    ## construct covariates used in the first-stage and the second-stage regression
  attach(dat3);
  set.seed(300); 
  H10<- cbind(1, O11);
  colnames(H10)<-c("int","O11");

  H11<- cbind(1, O11);
  colnames(H11)<-c("A1","A1O11");
Y1<- rep(0,200); # there is no Y1 in this simulated data
H20<- cbind(1,O11,O21,O22);
colnames(H20)<-c("int","O11","O21","O22");
H21<- cbind(1,O21);  #O21 is a candidate tailoring variable for stage 2
colnames(H21)<-c("A2","A2O21"); 
Y2<- Y; 
S<- rep(1,200); # everyone is randomized at stage 2
  
    ## Construct contrast matrices
    c1<-diag(4); 
    c2<-diag(6); 
  
    ## Run qlearning function to get estimates and confidence intervals for the contrasts

result3<-qlaci(H10, H11, A1, Y1, H20, H21, A2, Y2, S,c1=t(c1),c2=t(c2) ,nb=1000);
print(result3);
  detach(dat3);



    
   } 

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{  }
%\keyword{  }% __ONLY ONE__ keyword per line