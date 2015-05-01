#'Rolling example

#'The beta-binomial distribution can also be motivated via an urn model for positive integer values of $\alpha$ and $\beta$ , known as the Polya urn model. Specifically, imagine an urn containing $\alpha$ red balls and $\beta$ black balls, where random draws are made. If a red ball is observed, then two red balls are returned to the urn. Likewise, if a black ball is drawn, then two black balls are returned to the urn. If this is repeated n times, then the probability of observing k red balls follows a beta-binomial distribution with parameters n, $\alpha$ and $\beta$.

#'Note that if the random draws are with simple replacement (no balls over and above the observed ball are added to the urn), then the distribution follows a binomial distribution and if the random draws are made without replacement, the distribution follows a hypergeometric distribution.

#+
    ## Not run:
        # Data
library(DPpackage)
data(rolling)
#first column include the number of success and the second the number of trials
y <- cbind(rolling$y1,rolling$y2)
head(y)
hist(y)
                                        # Prior information
prior<-list(alpha=1,
            a1=1,
            b1=1)
                                        # Initial state
state <- NULL
                                        # MCMC parameters
mcmc <- list(nburn=5000,
             nsave=10000,
             nskip=3,
             ndisplay=20000)
#                                       Fitting the model
fit <- DPbetabinom(y=y,ngrid=100,
                   prior=prior,
                   mcmc=mcmc,
                   state=state,
                   status=TRUE)


#' This function fits a semiparametric Beta-Binomial model

#'$y_i|n_i,p_i \sim Binom(n_i, p_i), i =1, \ldots n$
#' 
#' $p_i|G \sim G$
#' 
#' $G|\alpha, G_0 \sim DP(\alpha G_0)$
#' 
#' where the baseline distribution is
#' 
#' $G_0 = Beta (a_1,b_1)$
#' 
#' To complete the model specification, $G_0$, is a conjugate prior in this model specification.

#+
fit
summary(fit)
                                        # density estimate
plot(fit,output="density")
                                        # parameters
plot(fit,output="param")
## End(Not run)


##############################################################################
# Density estimation examples
# DPdensity fits univariate and multivariate versions of 
# the conjugate model of Escobar and West (1995).
# PTdensity fits marginalized (Hanson and Johnson, 2002) and finite Polya
# tree MPTs (Hanson, 2006) to univariate and multivariate data.  For 
# multivariate, only the marginalized version is fit.
##############################################################################

#+
library(DPpackage)
help(DPdensity) # look at documentation
help(PTdensity)
data(galaxy)
help(galaxy)
attach(galaxy)
speeds=galaxy/1000

# MCMC specifications and prior for DP mixture
mcmc=list(nburn=1000,nsave=2000,nskip=19,ndisplay=10000)
m=mean(speeds[,1])
s=var(speeds[,1])
prior=list(a0=1,b0=1,nu1=3,nu2=1,s2=s,m2=m,psiinv2=solve(0.1*s),tau1=0.2,tau2=2)

# Fit DP mixture, plot, and obtain LPML
fit.dpm=DPdensity(y=speeds,prior=prior,mcmc=mcmc,state=NULL,status=TRUE)
summary(fit.dpm)
plot(fit.dpm)
plot(fit.dpm,output="param")
sum(log(fit.dpm$cpo)) # LPML statistic

# MCMC specifications and prior for MPT
mcmc=list(nburn=1000,nsave=2000,nskip=19,ndisplay=100,tune1=0.15,tune2=1.1,tune3=1.1)
prior=list(a0=5,b0=1,M=5)

# Fit MPT, plot, and obtain LPML
fit.mpt=PTdensity(y=speeds,prior=prior,mcmc=mcmc,state=state,status=TRUE)
summary(fit.mpt)
plot(fit.mpt)
plot(fit.mpt,output="param")
sum(log(fit.mpt$cpo)) # LPML statistic

##############################################################################
# Generalized additive model example: enthanol data.
# pord=order of penalty prior, pord=1 gives random walk prior in Chapter 15.
# degree=degree of B-spline basis functions, d=3 gives cubic.
# beta0 and Sbeta0 are mean and covariance of fixed effects: intercept and
# any other linear terms in the model.
# Penalty *prior* is the same across additive transformations lambda_i ~ gamma(taub1,taub2).
# For Gaussian responses, precision tau ~ gamma(tau1,tau2).
##############################################################################

#+
library(lattice)
data(ethanol)
attach(ethanol)
help(ethanol)
plot(ethanol)
help(PSgam)

# Model with only transformed E ( included
prior=list(taub1=0.01,taub2=0.01,beta0=rep(0,1),Sbeta0=diag(100,1),tau1=0.001,tau2=0.001)
mcmc=list(nburn=2000,nsave=2000,nskip=9,ndisplay=10000)
fit1=PSgam(formula=ethanol$NOx~ps(ethanol$E,k=20,degree=3,pord=1),
          family=gaussian(identity),prior=prior,mcmc=mcmc,ngrid=50,state=NULL,status=TRUE)
plot(fit1)
sum(log(fit1$cpo))

# Additive model with additive E and C functions
fit2=PSgam(formula=ethanol$NOx~ps(ethanol$E,ethanol$C,k=20,degree=3,pord=1),
          family=gaussian(identity),prior=prior,mcmc=mcmc,ngrid=50,state=NULL,status=TRUE)
plot(fit2)
sum(log(fit2$cpo))

# Model with transformed E but linear C
prior=list(taub1=0.01,taub2=0.01,beta0=rep(0,2),Sbeta0=diag(100,2),tau1=0.001,tau2=0.001)
fit3=PSgam(formula=ethanol$NOx~ps(ethanol$E,k=20,degree=3,pord=1)+ethanol$C,
          family=gaussian(identity),prior=prior,mcmc=mcmc,ngrid=50,state=NULL,status=TRUE)
plot(fit3)
summary(fit3)
sum(log(fit3$cpo))

# Model with transformed E but categorical C
prior=list(taub1=0.01,taub2=0.01,beta0=rep(0,5),Sbeta0=diag(100,5),tau1=0.001,tau2=0.001)
fit4=PSgam(formula=ethanol$NOx~ps(ethanol$E,k=20,degree=3,pord=1)+factor(ethanol$C),
          family=gaussian(identity),prior=prior,mcmc=mcmc,ngrid=50,state=NULL,status=TRUE)
plot(fit4)
summary(fit4)
sum(log(fit4$cpo))

##############################################################################
# Generalized additive model example: O-Ring data of Table 8.1
# B-spline seems to overfit unless one is careful.  
# MCMC mixing can also be problematic.
##############################################################################

#+
failure=c(1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0)
temp=c(53,57,58,63,66,67,67,67,68,69,70,70,70,70,72,73,75,75,76,76,78,79,81)

mcmc=list(nburn=2000,nsave=2000,nskip=9,ndisplay=10000,tune=1.1)
prior=list(beta0=rep(0,2),Sbeta0=diag(1000,2))
fit.logit=  Pbinary(failure~temp,link="logit",  prior=prior,mcmc=mcmc,state=NULL,status=TRUE)
fit.cloglog=Pbinary(failure~temp,link="cloglog",prior=prior,mcmc=mcmc,state=NULL,status=TRUE)
fit.probit= Pbinary(failure~temp,link="probit", prior=prior,mcmc=mcmc,state=NULL,status=TRUE)
sum(log(fit.logit$cpo))
sum(log(fit.cloglog$cpo))
sum(log(fit.probit$cpo))

mcmc=list(nburn=5000,nsave=5000,nskip=19,ndisplay=10000)
prior=list(taub1=0.01,taub2=0.01,beta0=rep(0,1),Sbeta0=diag(1000,1))
fit1.additive1=PSgam(formula=failure~ps(temp,k=20,degree=3,pord=1),
          family=binomial(logit),prior=prior,mcmc=mcmc,ngrid=50,state=NULL,status=TRUE)
sum(log(fit1.additive$cpo))


prior=list(taub1=1,taub2=0.01,beta0=rep(0,1),Sbeta0=diag(1000,1))
fit2.additive=PSgam(formula=failure~ps(temp,k=20,degree=3,pord=1),
          family=binomial(logit),prior=prior,mcmc=mcmc,ngrid=50,state=NULL,status=TRUE)
sum(log(fit2.additive$cpo))

plot(fit1.additive)
plot(fit2.additive)





