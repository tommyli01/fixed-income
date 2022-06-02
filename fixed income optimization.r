
###Q2
###assume spot rate curve:
ttm.yrs<-seq(0.5,5,0.5)
realspot<-c(0.03,0.035,0.04,0.0425,0.045,0.0475,0.05,0.05125,0.0525,0.0525)
curve.table<-data.frame(ttm.yrs,realspot)
###instruction: In the Rendleman-Bartter model (Bloomberg's "fair value model"), 
###the one-period spot rate follows the same process 
###as the Black-Scholes model for stock prices.
###These spot rates are based on semi-annual compounding. 
###Use optimx to find the two parameters, (lambda,sigma), 
###to calibrate the Rendleman-Bartter model
###packages:
library("optimx")
###functions
prt.tree <- function(tree,digit=2) {
  nt <- nrow(tree)
  # transpose tree
  trantree <- t(tree)
  nt1 <- 2*nt-1
  bintree <- matrix(rep("",nt1*nt),nrow=nt1,ncol=nt)
  # convert to bin tree
  i1 <- nt
  for (j in 1:nt) {
    i1 <- nt-j+1
    for (i in 1:j) {
      bintree[i1,j] <- as.character(round(trantree[i,j],digit))
      i1 <- i1 + 2
    }
  }
  rownames(bintree) <- rep("",nt1)
  colnames(bintree) <- rep("",nt)
  return(noquote(bintree))
}

# Calculate price of a zero coupon bond
#           paying $1 at time T (in years)
#           using up probability tree qtree
#           and   instantaneous interest rate tree ztree
#
BM_ZBP <- function(T,qtree,ztree,delt) {
  # Create price tree for zero-coupon bond paying $1 in year T
  # using a given instantaneous spot rate tree and up probability tree
  #
  # Args: 
  #   T:     The maturity of the zero coupon bond (in years) 
  #   qtree: The up probability tree
  #   ztree: The instantaneous spot rate tree
  #   delt:  length of a time step (in years)
  #
  # Returns:
  #   Price tree of the zero coupon bond
  
  N <- floor(T/delt)+1
  
  # make sure N < nrow(ztree)
  if (N>nrow(ztree)) return(matrix(NA,nrow=N,ncol=N))
  if (N>nrow(qtree)) return(matrix(NA,nrow=N,ncol=N))
  
  ptree <- matrix(0,nrow=N,ncol=N)
  ptree[N,1:N] <- 1    
  for (i in (N-1):1) {
    i1 <- i+1
    ptree[i,1:i] <- exp(-delt*ztree[i,1:i])*
      ( qtree[i,1:i]*ptree[i+1,1:i]+ (1-qtree[i,1:i])*ptree[i+1,2:i1])
  }
  return(ptree)
}

###setup
tmat <- seq(0.5,5,0.5)
r0 <- 0.03
sigma <- 0.01
lambda<- 0
delt <- 0.5
m <-2
parm <- c(lambda,sigma)
u <- exp(sigma*sqrt(delt))
d <- exp(-sigma*sqrt(delt))
q <- (exp(lambda*delt)-d)/(u-d)
deltm <- delt*m
### create one-period spot tree
N    <- 10   # number of steps
###ztree function
BM_ztree <- function(r0,u,d,N) {
  ztree  <- matrix(0, nrow=N+1, ncol=N+1)
  ztree[1,1] <- r0
  for (i in 2:N) {
    i1 <- i-1
    ztree[i,1:i1] <- ztree[i-1,1:i1]*u
    ztree[i,i] <- ztree[i-1,i-1]*d
  }
  return(ztree)
}
prt.tree(BM_ztree(r0,u,d,N),4)
ztree<-BM_ztree(r0,u,d,N)
###qtree function
BM_qtree <- function(r0,q,N) {
  qtree  <- matrix(0, nrow=N+1, ncol=N+1)
  qtree[1,1]<- q
  for (i in 2:N) {
    i1 <- i-1
    qtree[i,1:i1] <- qtree[i-1,1:i1]*q
    qtree[i,i] <- qtree[i-1,i-1]*(1-q)
  }
  return(qtree)
}
prt.tree(BM_qtree(r0,q,N),5)
qtree<-BM_qtree(r0,q,N)
### one period price
P1 = matrix(0, nrow=N, ncol=N)
for (i in 1:(N)) {
  P1[i,1:i] = 100/(1+ztree[i,1:i]/m)^deltm
}
prt.tree(P1)
### two period price
P2 = matrix(0, nrow=N-1, ncol=N-1)
for (i in 1:(N-1)) {
  i1=i+1
  P2[i,1:i] = (qtree[i,1:i]*P1[i+1,1:i] + 
                 (1-qtree[i,1:i])*P1[i+1,2:i1])/(1+ztree[i,1:i]/m)^deltm
}
prt.tree(P2)
### three period price
P3 = matrix(0, nrow=N-2, ncol=N-2)
for (i in 1:(N-2)) {
  i1=i+1
  P3[i,1:i] = (qtree[i,1:i]*P2[i+1,1:i]+
                 (1-qtree[i,1:i])*P2[i+1,2:i1])/(1+ztree[i,1:i]/m)^deltm
}
prt.tree(P3)
### four period price
P4 = matrix(0, nrow=N-3, ncol=N-3)
for (i in 1:(N-3)) {
  i1=i+1
  P4[i,1:i] = (qtree[i,1:i]*P3[i+1,1:i]+
                 (1-qtree[i,1:i])*P3[i+1,2:i1])/(1+ztree[i,1:i]/m)^deltm
}
prt.tree(P4)
### five period price
P5 = matrix(0, nrow=N-4, ncol=N-4)
for (i in 1:(N-4)) {
  i1=i+1
  P5[i,1:i] = (qtree[i,1:i]*P4[i+1,1:i]+
                 (1-qtree[i,1:i])*P4[i+1,2:i1])/(1+ztree[i,1:i]/m)^deltm
}
prt.tree(P5)
### six period price
P6 = matrix(0, nrow=N-5, ncol=N-5)
for (i in 1:(N-5)) {
  i1=i+1
  P6[i,1:i] = (qtree[i,1:i]*P5[i+1,1:i]+
                 (1-qtree[i,1:i])*P5[i+1,2:i1])/(1+ztree[i,1:i]/m)^deltm
}
prt.tree(P6)
### seven period price
P7 = matrix(0, nrow=N-6, ncol=N-6)
for (i in 1:(N-6)) {
  i1=i+1
  P7[i,1:i] = (qtree[i,1:i]*P6[i+1,1:i]+
                 (1-qtree[i,1:i])*P6[i+1,2:i1])/(1+ztree[i,1:i]/m)^deltm
}
prt.tree(P7)
### eight period price
P8 = matrix(0, nrow=N-7, ncol=N-7)
for (i in 1:(N-7)) {
  i1=i+1
  P8[i,1:i] = (qtree[i,1:i]*P7[i+1,1:i]+
                 (1-qtree[i,1:i])*P7[i+1,2:i1])/(1+ztree[i,1:i]/m)^deltm
}
prt.tree(P8)
### nine period price
P9 = matrix(0, nrow=N-8, ncol=N-8)
for (i in 1:(N-8)) {
  i1=i+1
  P9[i,1:i] = (qtree[i,1:i]*P8[i+1,1:i]+
                 (1-qtree[i,1:i])*P8[i+1,2:i1])/(1+ztree[i,1:i]/m)^deltm
}
prt.tree(P9)
### ten period price
P10 = matrix(0, nrow=N-9, ncol=N-9)
for (i in 1:(N-9)) {
  i1=i+1
  P10[i,1:i] = (qtree[i,1:i]*P9[i+1,1:i]+
                  (1-qtree[i,1:i])*P9[i+1,2:i1])/(1+ztree[i,1:i]/m)^deltm
}
prt.tree(P10)
###spot rate
price <- c(P1[1,1],P2[1,1],P3[1,1],P4[1,1],P5[1,1],
           P6[1,1],P7[1,1],P8[1,1],P9[1,1],P10[1,1])
price
spot <- c()
for (i in 1:length(price)){
  spot<- c(spot,2*((100/price[i])^(1/i)-1))
}
###calculating sse
spot
realspot<-c(0.03,0.035,0.04,0.0425,0.045,0.0475,0.05,0.05125,0.0525,0.0525)
data <- data.frame(cbind(ttm = tmat,realspot,modelspot = spot))
data
data$err <- round(data$realspot - data$modelspot,8)
data
round(sum(data$err^2),8)
###minimizing
parm
sse <- function(parm) {
  z = data[,2]
  t = data[,1]
  N = length(z)
  delt = 0.5
  r0=r0
  sigma = parm[2]
  lambda = parm[1]
  u <- exp(sigma*sqrt(delt))
  d <- exp(-sigma*sqrt(delt))
  q <- (exp(lambda*delt)-d)/(u-d)
  zmod = rep(r0,N)
  qtree <- BM_qtree(r0,q,N)
  ztree <- BM_ztree(r0,u,d,N)
  for (imat in 1:N) {
    yrs <- t[imat]
    ZBP   <- BM_ZBP(yrs,qtree,ztree, delt)
    zmod[imat] <- 2*((1/ZBP[1,1])^(1/(t[imat]*2))-1)
  }
  
  sse <- sum( (z-zmod)^2 )*1000000 
  # scale sse up by 1000000 for use in optimization
  
  return(sse)
  
}
ZBP<-c()
prt.tree(ztree,5)
for (imat in 1:N) {
  yrs <- tmat[imat]
  ZBP   <- c(ZBP, BM_ZBP(yrs,qtree,ztree, delt)[1,1])
}
ZBP
price
round(sse(parm)/1000000,8)
lb <- c(-Inf,0)
opt2 <- optimx(parm,sse,method=c("Nelder-Mead"),control=list(maxit=5000))
opt2
parm2 <- c(opt2$p1,opt2$p2)
round(parm2,4)
round(sse(parm2)/1000000,8)

table<-cbind(tmat,z = data[,2],zmod = spot,err = round(data[,2]-spot,4), price)
table

###end
