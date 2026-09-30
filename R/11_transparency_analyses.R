suppressPackageStartupMessages({library(data.table); library(survival); library(splines)})
setDTthreads(1)
z <- readRDS('work/revision_data.rds'); d <- copy(z$data)
out <- 'outputs/final/supplementary'
d[, `:=`(t60=pmin(months,60), e_os=os_event*(months<=60), e_css=css_event*(months<=60))]

# Reverse Kaplan-Meier follow-up, reported separately by histology.
fu <- survfit(Surv(months, 1-os_event) ~ histology, data=d, conf.type='log-log')
fut <- summary(fu)$table
fu_tab <- data.table(histology=sub('histology=','',rownames(fut)), n=as.integer(fut[,'records']),
                     median_followup_months=unname(fut[,'median']),
                     lower95=unname(fut[,'0.95LCL']), upper95=unname(fut[,'0.95UCL']))
fwrite(fu_tab, file.path(out,'TableS14_followup_by_histology.csv'))

# Bootstrap execution audit.
b <- readRDS('work/revision_bootstrap.rds'); bm <- b$replicates
boot_audit <- data.table(planned_replicates=2000L, completed_replicates=nrow(bm),
  complete_replicates=sum(complete.cases(bm)), failed_replicates=sum(!complete.cases(bm)),
  replicates_with_warnings=sum(lengths(b$warnings)>0), total_warnings=sum(lengths(b$warnings)),
  nonfinite_estimates=sum(!is.finite(bm)), seed=b$seed)
fwrite(boot_audit,file.path(out,'TableS15_bootstrap_execution_audit.csv'))

# Primary and landmark conditional-model PH diagnostics with term chi-square, df and p.
term_rows <- list(); summary_rows <- list()
record_zph <- function(fit,outcome,landmark,specification) {
  zz <- cox.zph(fit, transform='km', terms=TRUE, singledf=FALSE)
  tt <- as.data.table(zz$table, keep.rownames='term')
  setnames(tt,names(tt)[2:4],c('chisq','df','p'))
  tt[, `:=`(outcome=outcome, landmark=landmark, specification=specification)]
  term_rows[[length(term_rows)+1L]] <<- tt[,.(outcome,landmark,specification,term,chisq,df,p)]
  non_global <- tt[term!='GLOBAL']
  summary_rows[[length(summary_rows)+1L]] <<- data.table(outcome=outcome,landmark=landmark,
    specification=specification, global_chisq=tt[term=='GLOBAL',chisq],
    global_df=tt[term=='GLOBAL',df], global_p=tt[term=='GLOBAL',p],
    terms_p_lt_0_05=sum(non_global$p<.05,na.rm=TRUE), terms_tested=nrow(non_global),
    violating_terms=paste(non_global[p<.05,term],collapse='; '))
}
for (oo in c('OS','CSS')) {
  ev <- if(oo=='OS') 'e_os' else 'e_css'
  fp <- coxph(as.formula(paste0('Surv(t60,',ev,') ~ ',z$covars,' + strata(histology)')),
              data=d, ties='breslow', x=TRUE)
  record_zph(fp,oo,0,'Primary 0-60-month histology-stratified model')
  raw_ev <- if(oo=='OS') 'os_event' else 'css_event'
  for(L in c(0,12,24)) {
    dl <- d[months>L]
    dl[, `:=`(ltime=pmin(months-L,36), levent=get(raw_ev)*(months<=L+36))]
    f1 <- coxph(as.formula(paste0('Surv(ltime,levent) ~ ',z$covars,' + strata(histology)')),
                data=dl,ties='efron',x=TRUE)
    f2 <- coxph(as.formula(paste0('Surv(ltime,levent) ~ ',z$covars,' + strata(histology,stage)')),
                data=dl,ties='efron',x=TRUE)
    record_zph(f1,oo,L,'Landmark histology-stratified model')
    record_zph(f2,oo,L,'Landmark histology-and-stage-stratified model')
  }
}
fwrite(rbindlist(term_rows),file.path(out,'TableS9_conditional_model_diagnostics.csv'))
fwrite(rbindlist(summary_rows),file.path(out,'TableS9a_conditional_model_diagnostic_summary.csv'))

# Combine primary and alternative conditional survival estimates in one transparent table.
pr <- fread('outputs/final/standardized_conditional_survival_ASC_target.csv')
pr[, outcome := fifelse(outcome=='os_event','OS','CSS')]
pr[, specification := 'Primary 0-60-month histology-stratified model']
pr2 <- pr[,.(outcome,landmark=landmark_month,specification,ASC=standardized_ASC,
             AC_NOS=standardized_AC_NOS,difference,lower95,upper95,interval='Bootstrap 95% CI')]
alt <- fread(file.path(out,'TableS10_conditional_model_sensitivity.csv'))
alt[, outcome := fifelse(outcome=='os_event','OS','CSS')]
alt2 <- alt[,.(outcome,landmark,specification,ASC,AC_NOS,difference,
              lower95=NA_real_,upper95=NA_real_,interval='Point estimate only')]
allm <- rbindlist(list(pr2,alt2),fill=TRUE); setorder(allm,outcome,landmark,specification)
fwrite(allm,file.path(out,'TableS10_all_conditional_model_estimates.csv'))

# Smooth time-varying histology coefficient from scaled Schoenfeld residuals.
smooth_rows <- list()
for (oo in c('OS','CSS')) {
  ev <- if(oo=='OS') 'e_os' else 'e_css'
  ff <- coxph(as.formula(paste0('Surv(t60,',ev,') ~ histology + ',z$covars)),
              data=d,ties='efron',x=TRUE)
  zz <- cox.zph(ff,transform='identity',terms=TRUE,singledf=TRUE)
  nm <- grep('histology',colnames(zz$y),value=TRUE)[1]
  keep <- is.finite(zz$x) & is.finite(zz$y[,nm]) & zz$time>0 & zz$time<=60
  x <- zz$time[keep]; y <- zz$y[keep,nm]
  # Equal-count bins stabilise the display; a weighted natural-spline regression smooths the coefficient.
  br <- unique(quantile(x,probs=seq(0,1,length.out=41),na.rm=TRUE,type=8))
  bin <- cut(x,breaks=br,include.lowest=TRUE,labels=FALSE)
  ag <- data.table(x=x,y=y,bin=bin)[,.(time=median(x),beta=mean(y),n=.N),by=bin]
  sm <- lm(beta ~ ns(time,df=4),data=ag,weights=n)
  grid <- data.table(time=seq(0.5,60,by=0.5))
  pp <- predict(sm,newdata=grid,se.fit=TRUE)
  grid[,`:=`(outcome=oo,log_HR=as.numeric(pp$fit),se=as.numeric(pp$se.fit))]
  grid[,`:=`(HR=exp(log_HR),lower95=exp(log_HR-1.96*se),upper95=exp(log_HR+1.96*se))]
  smooth_rows[[oo]] <- grid
}
smooth <- rbindlist(smooth_rows)
fwrite(smooth,file.path(out,'FigureS3_smooth_time_varying_histology_source.csv'))

# Nature-style R graphics: stage-specific conditional OS and smooth coefficient sensitivity.
cols <- c('Localized'='#0072B2','Regional'='#E69F00','Distant'='#D55E00')
st <- fread(file.path(out,'TableS7_stage_specific_conditional_OS.csv'))
st <- st[histology=='ASC' & stage %chin% names(cols)]
st[,`:=`(landmark=as.numeric(landmark),survival=100*survival,lower95=100*lower95,upper95=100*upper95)]
plot_stage <- function(file,type='png') {
  if(type=='png') png(file,1800,1250,res=240) else if(type=='tiff') tiff(file,1800,1250,res=240,compression='lzw') else svg(file,7.5,5.2)
  par(mar=c(4.4,4.7,1.2,1.2),las=1,mgp=c(2.7,.7,0),tcl=-.25,family='sans')
  plot(NA,xlim=c(-1,25),ylim=c(0,100),xaxt='n',xlab='Time already survived (months)',ylab='Probability of surviving a further 36 months (%)',bty='l')
  axis(1,at=c(0,12,24)); abline(h=seq(0,100,20),col='#E6E6E6',lwd=.8)
  for(s in names(cols)){q<-st[stage==s]; lines(q$landmark,q$survival,col=cols[s],lwd=2); arrows(q$landmark,q$lower95,q$landmark,q$upper95,angle=90,code=3,length=.045,col=cols[s]); points(q$landmark,q$survival,pch=21,bg='white',col=cols[s],lwd=1.5,cex=1.1); text(q$landmark,q$upper95+4,paste0('n=',q$n),cex=.72,col=cols[s])}
  legend('bottomright',names(cols),col=cols,lwd=2,pch=21,pt.bg='white',bty='n',cex=.85)
  mtext('Descriptive, unadjusted Kaplan-Meier estimates',side=3,adj=0,cex=.8,col='#555555')
  dev.off()
}
plot_stage('outputs/final/FigureS2_stage_conditional_OS.png','png'); plot_stage('outputs/final/FigureS2_stage_conditional_OS.tiff','tiff'); plot_stage('outputs/final/FigureS2_stage_conditional_OS.svg','svg')

plot_smooth <- function(file,type='png') {
  if(type=='png') png(file,1800,1250,res=240) else if(type=='tiff') tiff(file,1800,1250,res=240,compression='lzw') else svg(file,7.5,5.2)
  par(mar=c(4.4,4.7,1.2,1.2),las=1,mgp=c(2.7,.7,0),tcl=-.25,family='sans')
  plot(NA,xlim=c(0,60),ylim=c(.35,2.8),log='y',xlab='Months since diagnosis',ylab='Exploratory smoothed coefficient\n(hazard-ratio scale)',yaxt='n',bty='l')
  axis(2,at=c(.5,.75,1,1.5,2,2.5),labels=c('.5','.75','1','1.5','2','2.5')); abline(h=1,lty=2,col='#666666')
  cc<-c(OS='#0072B2',CSS='#D55E00')
  for(o in c('OS','CSS')){q<-smooth[outcome==o]; polygon(c(q$time,rev(q$time)),c(q$lower95,rev(q$upper95)),border=NA,col=adjustcolor(cc[o],alpha.f=.14));lines(q$time,q$HR,col=cc[o],lwd=2)}
  legend('topright',c('All-cause mortality','Cancer-specific mortality'),col=cc,lwd=2,bty='n',cex=.85)
  mtext('Scaled Schoenfeld residual diagnostic; not a directly fitted continuous-time HR',side=3,adj=0,cex=.76,col='#555555')
  dev.off()
}
plot_smooth('outputs/final/FigureS3_smooth_time_varying_histology.png','png'); plot_smooth('outputs/final/FigureS3_smooth_time_varying_histology.tiff','tiff'); plot_smooth('outputs/final/FigureS3_smooth_time_varying_histology.svg','svg')
cat('Transparency analyses completed\n'); print(fu_tab); print(boot_audit); print(rbindlist(summary_rows))
