#*************************************************************************************
#
# R functions and script for performing fish stock-recruitment analysis
# by Dr. Zhenming Su
#
# suz@michigan.gov
# IFR, DNR, State of Michigan
# 2005-2007(R)
#
#*************************************************************************************

# Data
lnr.s <- structure(list(lnRS = c(-0.490177076, 0.469642607, 0.490973641,
0.077094548, 0.086549933, -0.912923449, 1.292167433, 0.97449055,
-1.944904111, 0.932226166, 0.283556949, -1.461697353, 0.959053681,
0.360778969, 0.3990024, 0.545237575, 0.527538564, -0.653310053,
0.92012115, -2.379617605, 2.601512275, -0.862384841, 1.106699081,
-0.456328691, 1.279392573, -2.654092655, 2.239375003, -0.849244046,
1.112050985, -1.538486548, 0.504073192, 0.823506581, 0.742526149,
-1.055453923, -0.28554639, 1.142742903, -0.096385854, 1.415007779,
-2.14918348, 1.934049889), spawner = c(80.47, 49.29, 63.07, 82.44,
71.24, 62.14, 19.95, 72.64, 192.48, 27.53, 69.92, 92.84, 21.52,
44.93, 51.56, 61.47, 84.83, 115.02, 59.85, 150.19, 13.91, 150,
50.66, 122.57, 62.13, 178.65, 12.57, 118.01, 50.48, 122.78, 21.09,
27.93, 50.91, 85.58, 29.78, 22.39, 70.19, 50.99, 167.93, 15.66
)), .Names = c("lnRS", "spawner"), class = "data.frame", row.names = c("1",
"2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
"14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
"25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35",
"36", "37", "38", "39", "40"))


gibbs<-function(start.values, lnr.s, MaxIterations){
  #MaxIterations	<- 7000	#number of MCMC iterations
  #BurnIn<-500	#no. of iterations to allow mixing of sequences
  #MinBurn<-200	minimum #iterations before sampling	
  #TestIterations<-100	#a check on convergence before sampling	
  #SampleSize<-5000	#posterior sample size you want	
  #Rcrit	1.2	#convergence criterion	
  #nchains	3	#number of MCMC sequences (fix at 3)	
  #npars	3	#number of model parameters (a,b,sigma^2)	
  #mobs	40	#number of observations	

  # Std Ricker model
  # lnRS = alpha + beta * spawner
  spawner <- lnr.s$spawner
  lnRS <- lnr.s$lnRS

  alpha <- start.values[1]
  beta <- start.values[2]
  s2 <- start.values[3]
  
  mobs <- length(lnRS)
  
  # prior of alpha: a ~ N(0, V.A)
  # prior of beta: b ~ N(0, V.B)
  V.A <- 10^6
  V.B <- 10^6
  # prior of s^2: s2 ~ IG(0.001, 0.001)
  theta.matrix <- matrix(NA, nr= MaxIterations, nc= length(start.values), byrow =TRUE)
  for (iter in 1: MaxIterations)
  {
    #Update stock-spec. parameters in-turn
    #  by their full conditional distributions.

    #update alpha
    mu.alpha <- mean(lnRS - beta * spawner) 
    VARA <- s2 / mobs
    mean.alpha <- mu.alpha * V.A / (V.A + VARA)
    v.alpha <- V.A * VARA / (V.A + VARA)
  
    alpha <- rnorm(1, mean.alpha, sqrt(v.alpha))
    
    #update beta
    sum.s2 <-sum(spawner * spawner)  
    mu.beta <- sum((lnRS - alpha) * spawner) / sum.s2
    VARB <- s2 / sum.s2

    mean.beta <- mu.beta * V.B / (V.B + VARB)
    v.beta <- V.B * VARB / (V.B + VARB)
  
    beta <- rnorm(1, mean.beta, sqrt(v.beta))

    #update siga
    ssq <- sum((lnRS -(alpha + beta * spawner))^2)
    s2  <- 1 / rgamma(1, shape=(0.001 + mobs / 2), scale = 1/(0.001 + ssq / 2))
    theta.matrix[iter,] <- c(alpha, beta, s2)
  } #Next iter
  theta.matrix
}

res <- gibbs(c(1, 0.1, 1),lnr.s,5000)
# trace plots
par(mfrow=c(3,1))
plot(res[,1], type='l', col = "red")
plot(res[,2], type='l', col = "green")
plot(res[,3], type='l', col = "blue")
#discard first 100 values
res <- res[100:length(res[,1]),]
#summaries
apply(res, 2, mean)
apply(res, 2, sd)
#correlations
cor(res)
#scatter plot
pairs(res)


