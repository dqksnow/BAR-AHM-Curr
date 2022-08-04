library(parallel)
options(scipen=200)
N=100 # replication times 
p=8 # number of covariates
bar  <- function(q){
  library(MASS)
  library(survival)
  library(mvtnorm)
  n <- 1500 #number of sample
  p <- 8 #number of covariaties
  lambda_0 <- 2
  lambda_c0 <- 0.1
  #COX =matrix(0,2,2)                ## gerenate the multinormal distribution variables matrices##
  #diag(COX) = c(rep(12,dim(COX)[1]))
  #beta_0 <- c(rep(0.7,2),rep(0,p-4),rep(0.7,2)) 
  beta_0 <- c(1, 0, 1, 0, 1, 0, 1, 0)
  #beta_0 <- c(rep(beta,3),rep(0,12))
  #gamma_0<- rep(0,p)
  MSE.ALASSO <- TP.ALASSO <- FP.ALASSO <- c() 
  #Z <- mvrnorm(2*n,rep(0,p),VarCovS)
  #Z=matrix(runif(n*p,0,sqrt(2)),n,p)
  Z = matrix(NA,n,p)
  
  repeat { 
    Z[,1] <- rnorm(n,0,1)
    Z[,2] <- rnorm(n,0,1)
    
    Z[,3] <- rpois(n,1)
    Z[,4] <- rexp(n,1)
    
    Z[,5] <- rpois(n,1)
    Z[,6] <- rexp(n,1)
    
    Z[,7] <- rpois(n,1)
    Z[,8] <- rexp(n,1)
    if( floor(min(Z[,1:8]%*%beta_0) )>= -lambda_0) {
      break
    }
  }
  pred0   = Z %*% beta_0
  pred    = pred0 + lambda_0
  set.seed(q)
  time <- rexp(n,pred)#参数he越小，time越大
  set.seed(q)
  #pred_c = exp(pred0)*lambda_c0
  #set.seed(q)
  #mc<- rexp(n,pred_c)
  mc<- runif(n,0,0.33)
  #mc<- runif(n,0,1)
  #delta_c <- rep(1,n) #所有C都是已知的、观测到的
  delta <- rep(0,n)
  delta[which(time>=mc)] = 1
  
  
  ####################################
  #partial likelihood function
  logl  <- function(beta_e){
    A=exp(-Z%*%beta_e*mc)
    B=rep(NA,n)
    for (j in 1:n){
      B[j]= sum( exp(-(Z[mc>=mc[j], ]%*%beta_e*mc[j])) )
    }
    return(sum(delta*log(A/B))) #log(prod( (A/B)^delta) )
  }
  ###############################
  logl.BAR  <- function(beta_e){
    A=exp(-Z%*%beta_e*mc)
    B=rep(NA,n)
    for (j in 1:n){
      B[j]= sum( exp(-(Z[mc>=mc[j], ]%*%beta_e*mc[j])) )
    }
    #pen <- c((beta_e^2)/(bet_iteraive^2))
    #print(pen)
    return(sum(delta*log(A/B))-n*sum(lamb*(beta_e^2)/(bet_iteraive^2)))
  }
  n.lamb <- 5   # Total number of lambda's to be considered
  lamb.seq <- seq(0.00040,0.00055,,n.lamb)
  BIC  <- rep(NA,n.lamb)
  bet_BIC <- matrix(0, n.lamb, p)
  for (i in 1:n.lamb){
    lamb <- lamb.seq[i]
    set.seed(1)
    bet <- rep(0,p)
    bet_iteraive <- optim(bet, logl, method = 'BFGS',  control=c(fnscale=-1))$par #initial value
    bet.all <- bet_iteraive
    for (m in 1:10) {
      bet_iteraive <- optim(bet_iteraive, logl.BAR,method = 'BFGS' , control=c(fnscale=-1))$par
      bet.all <- rbind(bet.all,bet_iteraive) #bet.all involves bet_iteraive in every iteration
      if (abs(mean(bet.all[m,]-bet_iteraive))<1e-4) { break }
    }
    
    bet = bet_iteraive
    bet[abs(bet) < 0.01] <- 0
    BIC[i]=-2*logl(bet)+length(which(bet!=0))*log(n)
    bet_BIC[i,] <- bet
  }
  lamb <- lamb.seq[which(BIC==min(BIC))[1]] #in order to avoid the same BIC in different lamb, choose BIC==min(BIC))[1]
  bet_optim <- bet_BIC[which(BIC==min(BIC))[1], ]
  
  x<-rep(1,p)
  VarCovS<-diag(x)
  MSE.LASSO <- matrix((bet_optim-beta_0),1,p)%*%VarCovS%*%t(matrix((bet_optim-beta_0),1,p)) #MSE
  TP.LASSO <- sum(as.numeric(bet_optim[c(1,3,5,7)]!=0)) 
  FP.LASSO <- sum(as.numeric(bet_optim[c(2,4,6,8)]!=0)) 
  select_time=rep(0,p)
  select_time[which(bet_optim!=0)]=1
  return(list(bet_optim, MSE.LASSO, TP.LASSO, FP.LASSO, lamb, select_time))
}
cl <- makeCluster(12)
t1=proc.time()
results <- parLapply(cl, 1:N , bar)
t2=proc.time()
time_cost=t2-t1
print(paste0('cost time：',time_cost[3][[1]]/60,'mins'))

MSE.LASSO <- TP.LASSO <- FP.LASSO <- lamb_BIC <-  c() #results
select_times=matrix(0,N,p)
bet_optim=matrix(0,N,p)
for(i in 1:N){
  bet_optim[i,]=results[i][[1]][[1]]
  MSE.LASSO[i]=results[i][[1]][[2]]
  TP.LASSO[i]=results[i][[1]][[3]]
  FP.LASSO[i]=results[i][[1]][[4]]
  lamb_BIC[i]=results[i][[1]][[5]] #the lambda in the penality function choosen by BIC
  select_times[i,]=results[i][[1]][[6]]
}
print(paste('median(MSE.LASSO):',median(MSE.LASSO)))
print(paste('sd(MSE.LASSO):',sd(MSE.LASSO)))
print(paste('mean(TP.LASSO):',mean(TP.LASSO)))
print(paste('mean(FP.LASSO):',mean(FP.LASSO)))
#

bet_mean = colMeans(bet_optim)
bet_sd=rep(0,p)
for (i in 1:p){bet_sd[i]=sd(bet_optim[,i])} #
bet_mean
bet_sd
lamb_BIC
apply(select_times,2,sum)/N
boxplot(MSE.LASSO)
#setwd("E:/onedrive/当前状态数据可加风险模型的变量选择问题0311/0620/C_U_MSE")
#write.csv( MSE.LASSO, file = "MSE_bar.csv",row.names = FALSE)
