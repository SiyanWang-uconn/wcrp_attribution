# GEV helper checks -- interactive RStudio version
source(here::here("Code", "config", "paths.R"))
source(here::here("Code", "helpers", "gev.R"))

# 1. Quantile and survival-probability identity
for(shape in c(-.2,0,.2)) for(p in c(.5,.96,.99,.999)) {
  q <- qgev_local(p,30,2,shape)
  stopifnot(abs(exp(log_survival(q,30,2,shape))-(1-p))<1e-10)
}
stopifnot(log_survival(41,30,2,-.2)==-Inf)

# 2. No GMST effect implies RR = 1
f <- list(mu=30,alpha=0,sigma=2,xi=-.1,g0=0)
stopifnot(max(abs(attribution(f,1,25)$RR-1))<1e-10)
stopifnot(max(abs(attribution(f,1,25)$FAR))<1e-10)

# 3. Positive GMST effect implies RR > 1
f$alpha <- 3
stopifnot(all(attribution(f,1,25)$RR>1))
stopifnot(all(attribution(f,1,25)$FAR>0))

# 4. Recover known parameters from simulated data
set.seed(55)
g <- seq(-.5,1,length.out=2000)
x <- qgev_local(runif(length(g)),30+3*g,1.5,-.1)
fit <- fit_gev(x,g)
stopifnot(abs(fit$alpha-3)<.3,abs(fit$sigma-1.5)<.2,abs(fit$xi+.1)<.08)
cat('PASS: quantile/survival identity, support boundary, null RR, positive warming, simulated parameter recovery\n')
