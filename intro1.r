#'Rolling example

#+
    ## Not run:
        # Data
library(DPpackage)
data(rolling)
y <- cbind(rolling$y1,rolling$y2)
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
                                        # Fitting the model
fit <- DPbetabinom(y=y,ngrid=100,
                   prior=prior,
                   mcmc=mcmc,
                   state=state,
                   status=TRUE)
fit
summary(fit)
                                        # density estimate
plot(fit,output="density")
                                        # parameters
plot(fit,output="param")
## End(Not run)
