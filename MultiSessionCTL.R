#----------------------------------------------------------------------------------#
#Multi-session Control file as per Cahill et al. Lower Bow River Rainbow Trout CJFAS#
#Script simulates multi-session data, and then calls stan to analyze trend#
#See also Kery and Royle 2016 (Applied Hierarchical Modeling), or
#Royle et al. 2014 (Spatial Capture-Recapture) for BUGS Examples

#Run using R version 3.5.0 (2018-04-23) -- "Joy in Playing" &
#rstanarm (Version 2.17.4, packaged: 2018-04-13 01:51:52 UTC)
#Coded by Chris Cahill 16 May 2018#
#Contact: christopher.cahill@ucalgary.ca

#PLEASE NOTE!!!
#If stan has never been installed on your machine or you've had to update the 
#version of stan visit https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Windows
#and follow the instructions to configure stan on your computer before running script
#including the "optional configuration" re: the Makevars file

#Similarly, I had to update Makevars file as per last URL, restart R, and run script
#to get this running with the rstan 2.17.4
#----------------------------------------------------------------------------------#

#Set your local directory--'trend.stan' must be saved to this directory
setwd("C:/Users/Chris Cahill/Documents/GitHub/Bayesian-Mark-Recapture-Models/SubmissionRuns/Revision Code")

#Step 1: Simulate Multi-Session Lambdas & N's with a trend:
set.seed(7335)
G         <-  7 #Number of groups or sessions
Years     <-  1:G 
beta0     <-  5.5 #Intercept
beta1     <- -0.5 #Trend term
J         <-  4 #Number of ocassions for cap-recap each year

#Generate lambdas and simulate abundances:
lambdas   <- exp( beta0 + beta1*(Years - Years[4]) )
N         <- rpois(G, lambdas)
N

#Step 2: Generate capture histories based on N's:

#Use the N's to simulate observed data from model Mth
#This function is adapted from code in Kery and Schaub 2011 p.154
GetCapHist <- function(N=2000, J=4, mean.p=0.3, time.effects=rep(0, J), sd=1.25){
	yfull <- yobs <- p <- array(NA, dim=c(N, J))
	mean.lp <- log(mean.p/(1-mean.p)) 
	eps <- rnorm(N, 0, sd)
	
	for(j in 1:J){
		pp <- p[,j] <- plogis(mean.lp + time.effects[j] + eps)
		yfull[,j] <- rbinom(n=N, size=1, prob=pp)
	}
	
	ever.detected <- apply(yfull, 1, max)
	C <- sum(ever.detected)
	yobs <- yfull[ever.detected==1,]
	
  return(list(N=N, mean.lp=mean.lp, time.effects=time.effects, sd=sd, 
  	     eps=eps, C=C, J=J, yfull=yfull, yobs=yobs))
}

M            <- 2000 #Data augmentation variable
CapHistArray <- array(NA, dim=c(M, J, G))

#Simulate and augment capture history data:
for(i in 1:length(N)){
	SimData <- GetCapHist(N=N[i])$yobs
	CapHistArray[1:M ,1:J , i] <- 
	rbind(SimData,array(0, dim=c((M-nrow(SimData)[[1]]), dim(SimData)[[2]]) ) )
}

#------------------------------------------------------------------------------------------#
#Step 3: Load packages, set up data, and parameters to monitor, and then call Stan:
#------------------------------------------------------------------------------------------#

list.of.packages <- c("shinystan", "arm", "rstan", "rstanarm")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
require(shinystan)
require(arm)
require(rstan)
require(rstanarm)

#Option to detect and use multiple cores:
options(mc.cores = parallel::detectCores())

#Stan run parameters:
km = c( rep(4.32, 4), 4.22, 4.47, 4.51 ) #distances in km of river e-fished
data <- list(M = M, J = J, nyear = G, y_array = CapHistArray, X = 1:7, km = km) 
trend_pars=c("beta0", "beta1", "G_N", "mean_p") #can also monitor "psi" "G_N_km" "sigma"

#Stan run parameters:
iter   <- 1000 
thin   <- 1 
chains <- 4 #Should be fine for most laptops

if(FALSE){
     #Takes 45 minutes to run 1,000 iterations on CLC's desktop  --> Okay results
     #Takes 3.2 hours to run 10,000 iterations on CLC's desktop  --> Best results
     fit = stan(file = "trend.stan", model_name = "Trend", 
                data = data, iter = iter, thin = thin, chains = chains,
	        warmup = floor(iter/2), seed = 1,
                verbose = T, pars=trend_pars)
}

#Step 4: Ensure model is behaving and inspect results (no divergent transitions & Rhat < 1.1):

sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
divergent <- sapply(sampler_params, function(x) max(x[, "divergent__"]))

if(any(divergent == 1)){message("Model has divergent transitions and results should not be trusted!")}
if(any(summary(fit)$summary[,"Rhat"] > 1.1)){message("Chains may not have converged to a common distribution--Gelman-Rubin Stat > 1.1!")}

#Inspect results and compare to values used in simulation:
fit
beta0
beta1
N
#capture probabilities were all 0.3

#Launch shinystan in web browser to explore mixing of chains, etc.:
launch_shinystan(fit)
#------------------------------------------------------------------------------------------#
                                        #End End End#
#------------------------------------------------------------------------------------------#
