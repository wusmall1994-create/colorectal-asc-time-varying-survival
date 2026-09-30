suppressPackageStartupMessages({library(data.table);library(survival);library(splines)})
setDTthreads(1)
z<-readRDS('work/revision_data.rds'); d<-z$data
X<-model.matrix(as.formula(paste('~',z$covars)),d)[,-1,drop=FALSE]
H<-as.integer(d$histology); M<-d$months
Y<-lapply(c('os_event','css_event'),function(o)Surv(pmin(M,60),d[[o]]*(M<=60)))
idx<-split(seq_len(nrow(d)),H)
target<-which(H==2)
# Breslow frequency weights equal explicit patient replication; fixed spline knots.
one <- function(w=rep(1,nrow(d))) {
 k<-which(w>0); xx<-X[k,,drop=FALSE]; hh<-H[k]; ww<-w[k]
 ans<-list()
 for(o in 1:2){
  yy<-Y[[o]][k,,drop=FALSE]
  f<-coxph.fit(xx,yy,hh,offset=rep(0,length(k)),init=NULL,control=coxph.control(),weights=ww,method='breslow',rownames=NULL,resid=FALSE)
  if(any(!is.finite(f$coefficients)))stop('Invalid coefficient')
  risk<-exp(drop(xx%*%f$coefficients)); rr<-exp(drop(X[target,,drop=FALSE]%*%f$coefficients))
  g<-data.table(h=hh,t=yy[,1],events=ww*yy[,2],risk=ww*risk)[,.(events=sum(events),risk=sum(risk)),by=.(h,t)]
  setorder(g,h,t);g[,cumhaz:=cumsum(events/rev(cumsum(rev(risk)))),by=h]
  getH<-function(h,t){v<-g[h==..h & t<=..t,cumhaz]; if(length(v))tail(v,1)else 0}
  # Avoid data.table evaluation ambiguity for scalar query names.
  getH<-function(group,time){v<-g[h==group & t<=time,cumhaz];if(length(v))tail(v,1)else 0}
  for(lm in c(0,12,24)){
   sel<-M[target]>lm;tw<-w[target][sel];r<-rr[sel]
   a<-weighted.mean(exp(-(getH(2,lm+36)-getH(2,lm))*r),tw)
   c<-weighted.mean(exp(-(getH(1,lm+36)-getH(1,lm))*r),tw)
   ans[[length(ans)+1]]<-c(a,c,a-c)
  }
 }
 unlist(ans)
}
point<-one(); saveRDS(point,'work/revision_point.rds')
# Verify low-level implementation against public coxph/basehaz predictions.
d[,t60:=pmin(months,60)];d[,e60:=os_event*(months<=60)]
f<-coxph(as.formula(paste('Surv(t60,e60)~',z$covars,'+strata(histology)')),d,ties='breslow',x=TRUE)
bh<-basehaz(f,centered=FALSE);r<-exp(drop(X[target,,drop=FALSE]%*%coef(f)))
ref<-sapply(levels(d$histology),function(g)mean(exp(-tail(bh$hazard[bh$strata==g & bh$time<=36],1)*r)))
stopifnot(max(abs(point[1:2]-rev(ref)))<1e-8)
print(point);print('Public API cross-check passed')
B<-2000L;workers<-4L
cl<-parallel::makePSOCKcluster(workers)
parallel::clusterEvalQ(cl,{library(data.table);library(survival);setDTthreads(1);NULL})
parallel::clusterExport(cl,c('X','H','M','Y','idx','target','d','one'))
parallel::clusterSetRNGStream(cl,20260929)
bootfun<-function(b){
 w<-integer(nrow(d));for(j in idx)w[j]<-tabulate(sample.int(length(j),length(j),replace=TRUE),nbins=length(j))
 warnings<-character()
 value<-tryCatch(withCallingHandlers(one(w),warning=function(e){warnings<<-c(warnings,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)rep(NA_real_,18))
 list(value=value,warnings=warnings)
}
res<-list()
for(start in seq(1,B,by=40)){
 block<-parallel::parLapply(cl,start:min(start+39,B),bootfun)
 res<-c(res,block);saveRDS(res,'work/revision_bootstrap_checkpoint.rds')
 cat('Completed',length(res),'of',B,'at',format(Sys.time()),'\n');flush.console()
}
parallel::stopCluster(cl)
mat<-do.call(rbind,lapply(res,'[[','value'))
stopifnot(sum(complete.cases(mat))>=.99*B)
rows<-list()
for(o in 1:2)for(j in 1:3){i<-(o-1)*9+(j-1)*3
 rows[[length(rows)+1]]<-data.table(outcome=c('os_event','css_event')[o],landmark_month=c(0,12,24)[j],ASC_target_n=sum(M[target]>c(0,12,24)[j]),standardized_ASC=point[i+1],standardized_AC_NOS=point[i+2],difference=point[i+3],lower95=quantile(mat[,i+3],.025,na.rm=T),upper95=quantile(mat[,i+3],.975,na.rm=T),draws=sum(complete.cases(mat)))
}
fwrite(rbindlist(rows),'outputs/final/standardized_conditional_survival_ASC_target.csv')
saveRDS(list(point=point,replicates=mat,warnings=lapply(res,'[[','warnings'),seed=20260929),'work/revision_bootstrap.rds')
print(rbindlist(rows))
