# Figure 7-style plot -- interactive RStudio version.

# 1. Load paths and helpers
project_root <- "C:/Users/siyan/WCRP/wcrp_common_attribution_assessment"
path_path <- file.path(project_root,"Code", "config", "paths.R")

source(file.path(project_root,"Code", "helpers", "gev.R"))
source(file.path(project_root,"Code", "helpers", "workflow.R"))
source(path_path)


output_dir <- file.path(project_dir,"Results", "figure")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

datasets <- list(
  ERA5=list(dir='Results/ERA5', title='a) ERA5'),
  CPC=list(dir='Results/CPC', title='b) CPC')
  #Berkeley=list(dir='results_berkeley_tx', title='c) Berkeley'),
 # `E-Obs v33e`=list(dir='results_eobs_tx', title='d) E-Obs v33e')
)

return_period <- exp(seq(log(1.01), log(10000), length.out=250))
p <- 1 - 1/return_period
B <- figure7_bootstrap_B
set.seed(2026)

make_panel_data <- function(spec) {
  path <- file.path(project_dir, spec$dir)
  dat <- read.csv(file.path(path, 'annual_extremes_and_gmst.csv'))
  dat <- dat[dat$used_in_fit, ]
  fit <-fit_gev(dat$tmax, dat$gmst)
  rr <- read.csv(list.files(path, pattern='risk_ratios_PROVISIONAL.csv$', full.names=TRUE)[1])
  gf <- unique(rr$gmst_factual)
  gp <- -0.08645833 ### 1951 like
  mu_now <- fit$mu + fit$alpha * (gf - fit$g0)
  mu_past <- fit$mu + fit$alpha * (gp - fit$g0)
  now <- qgev_local(p, mu_now, fit$sigma, fit$xi)
  past <- qgev_local(p, mu_past, fit$sigma, fit$xi)

  # Shift each observed annual maximum to a common GMST before ranking.
  obs_now <- dat$tmax + fit$alpha * (gf - dat$gmst)
  obs_past <- dat$tmax + fit$alpha * (gp - dat$gmst)
  empirical_rp <- function(x) {
    ord <- order(x)
    pp <- ppoints(length(x))
    data.frame(rp=1/(1-pp), value=x[ord])
  }

  bootstrap_curves <- bootstrap_return_level_curves(
    fit, dat, gf, gp, p, B, seed = 20260922
  )
  threshold=dat$tmax[dat$year==2021]
  
  list(title=spec$title, now=now, past=past,
       now_ci=bootstrap_curves$current_ci,
       past_ci=bootstrap_curves$past_ci, obs_now=empirical_rp(obs_now),
       obs_past=empirical_rp(obs_past), threshold=dat$tmax[dat$year==2021],
       past_rp=exp(log_survival(threshold,
                    mu_past,fit$sigma,fit$xi)),
       success=sum(bootstrap_curves$success), n=nrow(dat))
}

panels <- lapply(datasets, make_panel_data)

draw_panel <- function(z) {
  # Set the visible range from fitted curves and observations. Very unstable
  # bootstrap tail fits are clipped at the panel edge rather than flattening
  # the scientifically relevant part of the plot.
  vals <- c(z$now, z$past, z$obs_now$value, z$obs_past$value)
  ylim <- range(vals[is.finite(vals)], na.rm=TRUE)
  pad <- diff(ylim)*.06
  plot(NA, xlim=c(1,15000), ylim=ylim+c(-pad,pad), log='x',
       xlab='Return period (years)', ylab='3-day maximum temperature (°C)',
       main=z$title, xaxt='n')
  axis(1, at=c(1,10,100,1000,10000), labels=c('1','10','100','1000','10000'))
  polygon(c(return_period,rev(return_period)),
          c(z$past_ci[1,],rev(z$past_ci[2,])), border=NA,
          col=adjustcolor('#4866ff',alpha.f=.14))
  polygon(c(return_period,rev(return_period)),
          c(z$now_ci[1,],rev(z$now_ci[2,])), border=NA,
          col=adjustcolor('#d73027',alpha.f=.12))
  lines(return_period,z$past,col='#001eff',lwd=2)
  lines(return_period,z$now,col='#bd241e',lwd=2)
  points(z$obs_past$rp,z$obs_past$value,pch=16,cex=.48,col='#001eff')
  points(z$obs_now$rp,z$obs_now$value,pch=16,cex=.48,col='#bd241e')
  abline(h=z$threshold,col='#ff33ff',lty=2,lwd=1.2)
  segments(25,par('usr')[3],25,par('usr')[3]+.025*diff(par('usr')[3:4]),
           col='#bd241e',lwd=3,xpd=NA)
  if (is.finite(z$past_rp) && z$past_rp <= 15000)
    segments(z$past_rp,par('usr')[3],z$past_rp,par('usr')[3]+.025*diff(par('usr')[3:4]),
             col='#001eff',lwd=3,xpd=NA)
  box()
}

draw_all <- function() {
  par(mfrow=c(2,1), mar=c(4.2,4.5,2.5,1), oma=c(1,1,2,1), las=1)
  draw_panel(panels$ERA5)
  draw_panel(panels$CPC)
  #draw_panel(panels$Berkeley)
  #draw_panel(panels[['E-Obs v33e']])
  mtext('Annual TXx, present vs 1951-like climate',
        outer=TRUE,side=3,font=2,cex=1.05)
}

# 2. Draw interactively first
draw_all()

# 3. Export PNG, PDF and summary
png(file.path(output_dir,'figure7_tx3x_four_datasets.png'),width=1900,height=1500,res=200)
draw_all(); dev.off()
pdf(file.path(output_dir,'figure7_tx3x_four_datasets.pdf'),width=9.5,height=7.5)
draw_all(); dev.off()

summary <- do.call(rbind,lapply(names(panels),function(nm)
  data.frame(dataset=nm,fit_years=panels[[nm]]$n,
             bootstrap_success=panels[[nm]]$success,
             bootstrap_requested=B,threshold_25yr=panels[[nm]]$threshold,
             past_return_period=panels[[nm]]$past_rp)))
write.csv(summary,file.path(output_dir,'plot_summary.csv'),row.names=FALSE)
print(summary)
