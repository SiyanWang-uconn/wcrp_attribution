# GEV convention: F(x)=exp(-(1+xi*(x-mu)/sigma)^(-1/xi)).
qgev_local <- function(p,mu,sigma,xi) {
  stopifnot(all(p>0 & p<1),sigma>0)
  z <- -log(p)
  if(abs(xi)<1e-7) mu-sigma*log(z) else mu+sigma*expm1(-xi*log(z))/xi
}
log_survival <- function(x,mu,sigma,xi) {
  z <- (x-mu)/sigma
  if(abs(xi)<1e-7) log_a <- -z else {
    h <- 1+xi*z
    if(h<=0) return(if(xi<0) -Inf else 0)
    log_a <- -log(h)/xi
  }
  if(log_a < -35) return(log_a)
  if(log_a > 35) return(0)
  log(-expm1(-exp(log_a)))
}
fit_gev <- function(x,g) {
  stopifnot(length(x)==length(g),all(is.finite(x)),all(is.finite(g)))
  g0 <- mean(g); gc <- g-g0
  reg <- lm(x~gc); s <- max(sd(residuals(reg))*sqrt(6)/pi,0.05)
  nll <- function(p) {
    sigma <- exp(p[3]); xi <- p[4]; z <- (x-p[1]-p[2]*gc)/sigma
    if(!is.finite(sigma) || sigma<1e-6 || any(!is.finite(z))) return(1e100)
    if(abs(xi)<1e-7) val <- sum(log(sigma)+z+exp(-z)) else {
      h <- 1+xi*z
      if(any(h<=0)) return(1e100)
      val <- sum(log(sigma)+(1+1/xi)*log(h)+exp(-log(h)/xi))
    }
    if(is.finite(val)) val else 1e100
  }
  candidates <- lapply(c(-.25,-.1,0,.1),function(shape) {
    optim(c(unname(coef(reg)[1])-0.5772157*s,unname(coef(reg)[2]),log(s),shape),
      nll,method='Nelder-Mead',control=list(maxit=5000,reltol=1e-10))
  })
  good <- vapply(candidates,function(f) f$convergence==0 && f$value<1e90,TRUE)
  if(!any(good)) stop('GEV optimizer failed')
  f <- candidates[good][[which.min(vapply(candidates[good],function(z) z$value,0.))]]
  if(f$par[4]<=-.5) stop('GEV shape <= -0.5: nonregular fit; inspect data/model')
  list(mu=f$par[1],alpha=f$par[2],sigma=exp(f$par[3]),xi=f$par[4],g0=g0,nll=f$value)
}

attribution <- function(f,gf,gc,dat) {
  muf <- f$mu+f$alpha*(gf-f$g0)
  threshold<- dat$tmax[dat$year==2021]
  #threshold <- qgev_local(1-1/T,muf,f$sigma,f$xi)
  delta <- gf-gc
  log_pc <- vapply(delta,function(d) log_survival(threshold,muf-f$alpha*d,f$sigma,f$xi),0.)
  log_pf <- log_survival(threshold,muf,f$sigma,f$xi)
  log_rr <- log_pf-log_pc
  data.frame(reference=c('1951-like'),delta_G=delta,return_period=exp(-log_pf),
    gmst_factual=gf,threshold_C=threshold,p_factual=exp(log_pf),p_counterfactual=exp(log_pc),
    log_RR=log_rr,RR=exp(log_rr),FAR=1-exp(-log_rr),
    intensity_change=f$alpha*delta)
}
