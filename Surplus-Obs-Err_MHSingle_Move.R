#*************************************************************************************
# Hybrid Gibbs and Metropis Sampling method for fitting 
#  Surplus production model with observation error
# 
# Author: Dr. Zhenming Su
#
# suz@michigan.gov
# IFR, DNR, State of Michigan
# Sept. 2007
#
#*************************************************************************************

"hakecpue.catch" <-
structure(list(Year = as.integer(c(1965, 1966, 1967, 1968, 1969,
1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980,
1981, 1982, 1983, 1984, 1985, 1986, 1987)), Cpue = c(1.78, 1.31,
0.91, 0.96, 0.88, 0.9, 0.87, 0.72, 0.57, 0.45, 0.42, 0.42, 0.49,
0.43, 0.4, 0.45, 0.55, 0.53, 0.58, 0.64, 0.66, 0.65, 0.63), Catch = as.integer(c(94,
212, 195, 383, 320, 402, 366, 606, 378, 319, 309, 389, 277, 254,
170, 97, 91, 177, 216, 229, 211, 231, 223))), .Names = c("Year",
"Cpue", "Catch"), class = "data.frame", row.names = c("1", "2",
"3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
"15", "16", "17", "18", "19", "20", "21", "22", "23"))

schaefer.model.Ypre <- function(para)
{
# get the catch data
#Catch
  Catch <- hakecpue.catch[,3]
  Y <- length(Catch)
  B <- numeric(Y)
  U <- numeric(Y)

  r <- para[1]
  K <- para[2]
  q <- para[3]
  
  B[1] <-  K
  for (t in 1:(Y-1))
  {
      U[t] <- (q * B[t])
      B[t+1] <- B[t] + r * B[t] * (1 - B[t]/K) - Catch[t]
      if (B[t+1] < 0) B[t+1] <- 1
  }
  U[Y] <- q * B[Y]
  U
}

sumsq <- function(para,logerror = TRUE)
{
   cpue <- hakecpue.catch[,2]
   ypre <- schaefer.model.Ypre(para)
   sumsq <- sum((log(cpue) - log(ypre))^2)
}

sigma2.update <- function(nd,theta)
{
  ssq <- sumsq(theta, TRUE)
  sigma2  <- 1 / rgamma(1, shape=(0.001 + nd / 2), scale = 1/(0.001 + ssq / 2))
}

theta.update <- function(theta, sigma2, p.jump, sd.jump, LB, UB)
{
 np <- length(theta)
 theta.old <- theta
 for (i in 1:np)
 {
  theta.star <- rnorm(1, theta.old[i], sd.jump[i])
  if (theta.star < LB[i] | theta.star > UB[i])
    p.jump[i] <- 0
  else {
    log.post.old <- (-sumsq(theta.old, TRUE)/(2*sigma2)) #-log(theta.old[i]))
    theta.old[i] <- theta.star
    log.post.star <- (-sumsq(theta.old, TRUE)/(2*sigma2)) #-log(theta.old[i]))
    r <- exp(log.post.star - log.post.old)
    theta[i] <- ifelse (runif(1) < r, theta.star, theta[i])
    p.jump[i] <- min(r,1)
  }
 }
 list(theta=theta, p.jump = p.jump)
} 

gibbs <- function(theta.vec, sigma2,niter){

nd <- length(hakecpue.catch[,1])
sd.jump <- 2.4 * c(0.043732, 600, 0.000049)

p.jump <- numeric(3)

theta <- theta.vec

theta.matrix <- matrix(data = NA, nrow = niter, ncol = 3, byrow = TRUE)
p.jump.matrix <- matrix(data = NA, nrow = niter, ncol = 3, byrow = TRUE)
sig2 <- numeric(niter)
for (i in 1: niter)
{
  #update siga
  sigma2 <- sigma2.update(nd,theta)
  res <- theta.update(theta, sigma2, p.jump, sd.jump, c(0,0,0), c(1,10000,1))
  theta <- res$theta
  p.jump <- res$p.jump
  theta.matrix[i,] <- theta
  p.jump.matrix[i,] <- p.jump
  sig2[i] <- sigma2
} #Next iter

list(theta = theta.matrix, sig2 = sig2, mean.sig2 = mean(sig2), mean = apply(theta.matrix,2, mean),
    sd = apply(theta.matrix,2, sd), p.jump.mean = apply(p.jump.matrix,2, mean))
}

res <- gibbs(c(0.4, 2700, 0.0004), 1, 5000)
res$p.jump.mean
res$mean
res$sd
pairs((res$theta))
par(mfrow =c(2,2))
plot(density(res$theta[,1]))
plot(density(res$theta[,2]))
plot(density(res$theta[,3]))



