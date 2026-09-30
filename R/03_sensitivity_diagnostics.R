suppressPackageStartupMessages({library(data.table);library(survival);library(splines)})
setDTthreads(1); z<-readRDS('work/revision_data.rds');d<-z$data;covars<-z$covars
out<-'outputs/final/supplementary'; warnings<-list()
capfit<-function(form,q,label,ties='efron')withCallingHandlers(coxph(form,data=q,ties=ties,x=TRUE),warning=function(w){warnings[[length(warnings)+1]]<<-data.table(model=label,message=conditionMessage(w));invokeRestart('muffleWarning')})
# Risk sets and survivor characteristics.
risk<-rbindlist(lapply(list(c(0,6),c(6,12),c(12,24),c(24,60)),function(ab)d[months>ab[1],.(at_risk=.N,all_deaths=sum(os_event*(months<=ab[2])),cancer_deaths=sum(css_event*(months<=ab[2])),unknown_deaths=sum(cod_unknown*(months<=ab[2])),start=ab[1],end=ab[2]),by=histology]))
fwrite(risk,'outputs/final/Table2_risk_sets_events.csv')
comp<-rbindlist(lapply(c(0,12,24),function(lm){q<-d[months>lm];r<-q[,.(variable='Age',level='Mean (SD)',n=.N,value=sprintf('%.1f (%.1f)',mean(age),sd(age)),landmark=lm),by=histology];for(v in c('sex','race','period','stage','grade','subsite','sequence_group'))for(lev in levels(q[[v]]))r<-rbind(r,q[,.(variable=v,level=lev,n=sum(get(v)==lev),value=sprintf('%d (%.1f%%)',sum(get(v)==lev),100*mean(get(v)==lev)),landmark=lm),by=histology]);r}))
fwrite(comp,file.path(out,'TableS5_landmark_composition.csv'))
# Competing risks and stage-specific descriptive conditional OS.
cif<-list();stages<-list()
for(lm in c(0,12,24))for(g in levels(d$histology)){
 q<-copy(d[months>lm & histology==g]);q[,tt:=months-lm];q[,status:=factor(fifelse(css_event==1,'Cancer',fifelse(other_event==1,'Other',fifelse(cod_unknown==1,'Unknown','Censored'))),levels=c('Censored','Cancer','Other','Unknown'))]
 f<-survfit(Surv(tt,status)~1,q);s<-summary(f,times=36,extend=FALSE)
 for(j in 2:4)cif[[length(cif)+1]]<-data.table(landmark=lm,histology=g,n=nrow(q),cause=f$states[j],probability=s$pstate[1,j],lower95=s$lower[1,j],upper95=s$upper[1,j])
 for(st in c('Localized','Regional','Distant')){
  qq<-q[stage==st];ss<-summary(survfit(Surv(tt,os_event)~1,qq),times=36,extend=FALSE)
  stages[[length(stages)+1]]<-data.table(landmark=lm,histology=g,stage=st,n=nrow(qq),deaths36=sum(qq$os_event*(qq$tt<=36)),survival=ss$surv,lower95=ss$lower,upper95=ss$upper,at_risk36=ss$n.risk)
 }
}
fwrite(rbindlist(cif),file.path(out,'TableS6_competing_risk_landmarks.csv'));fwrite(rbindlist(stages),file.path(out,'TableS7_stage_specific_conditional_OS.csv'))
# Correct selection across histologies, before comparator removal.
s<-readLines('work/colorectal_formal_preanalysis.R');s<-s[seq_len(which(s=='cox_results <- list()')-1)]
i<-grep('all <- all[!(histology',s,fixed=TRUE)
s[i]<-'setorder(all,id,year,sequence_order,histology); all <- all[, .SD[1], by=id]'
e<-new.env();eval(parse(text=s),e);earliest<-e$all
ov<-z$pre_cross[, .N,by=id][N>1,id]
fwrite(z$pre_cross[id%in%ov,.(histology,year,sequence_order,months)][order(year,sequence_order)],'work/overlap_records_without_ids.csv')
flow<-data.table(selection=c('ASC preferred','Earliest eligible tumour'),ASC=c(sum(d$histology=='ASC'),sum(earliest$histology=='ASC')),AC_NOS=c(sum(d$histology=='AC_NOS'),sum(earliest$histology=='AC_NOS')))
fwrite(flow,file.path(out,'TableS8_alternative_cohort_counts.csv'))
variant<-list()
for(v in c('Known stage and grade','Earliest tumour overall'))for(o in c('os_event','css_event'))for(ab in list(c(0,6),c(6,12),c(12,24),c(24,60))){
 q<-copy(if(v=='Known stage and grade')droplevels(d[stage!='Unknown' & grade!='Unknown'])else earliest)
 q<-q[months>ab[1]];q[,tt:=pmin(months-ab[1],ab[2]-ab[1])];q[,ee:=get(o)*(months<=ab[2])]
 f<-capfit(as.formula(paste('Surv(tt,ee)~histology+',covars)),q,paste(v,o,ab[1]))
 b<-coef(f)['histologyASC'];se<-sqrt(vcov(f)['histologyASC','histologyASC'])
 variant[[length(variant)+1]]<-data.table(analysis=v,outcome=o,interval_start=ab[1],interval_end=ab[2],n=nrow(q),events=sum(q$ee),HR=exp(b),lower95=exp(b-1.96*se),upper95=exp(b+1.96*se),p=2*pnorm(-abs(b/se)))
}
fwrite(rbindlist(variant),'outputs/final/additional_sensitivity_HR.csv');fwrite(rbindlist(variant),file.path(out,'TableS2b_additional_sensitivity_analyses.csv'))
# Diagnose original sparse interactions; exclude them from inference in revision.
for(o in c('os_event','css_event'))for(mod in c('stage','subsite','period')){
 q<-copy(d);q[,tt:=pmin(months,24)];q[,ee:=get(o)*(months<=24)]
 f<-capfit(as.formula(paste('Surv(tt,ee)~histology+',covars,'+histology:',mod)),q,paste('original interaction',o,mod))
 if(any(abs(coef(f))>10))warnings[[length(warnings)+1]]<-data.table(model=paste(o,mod),message=paste(names(coef(f))[abs(coef(f))>10],collapse=';'))
}
# Conditional model diagnostics and specification sensitivity.
diag<-list();pred<-list()
for(o in c('os_event','css_event'))for(lm in c(0,12,24))for(spec in c('Landmark histology strata','Landmark histology and stage strata')){
 q<-copy(d[months>lm]);q[,tt:=pmin(months-lm,36)];q[,ee:=get(o)*(months<=lm+36)]
 rhs<-if(spec=='Landmark histology strata')paste(covars,'+strata(histology)')else paste(gsub(' + stage','',covars,fixed=T),'+strata(histology,stage)')
 f<-capfit(as.formula(paste('Surv(tt,ee)~',rhs)),q,paste(spec,o,lm),ties='breslow')
 ph<-cox.zph(f);diag[[length(diag)+1]]<-data.table(outcome=o,landmark=lm,specification=spec,term=rownames(ph$table),p=ph$table[,'p'])
 bh<-basehaz(f,centered=FALSE);ta<-copy(q[histology=='ASC']);mm<-f$x[q$histology=='ASC',names(coef(f)),drop=FALSE];rr<-exp(drop(mm%*%coef(f)))
 ps<-sapply(c('ASC','AC_NOS'),function(g){hg<-if(spec=='Landmark histology strata')rep(g,nrow(ta))else paste(g,ta$stage,sep=', ');hs<-vapply(hg,function(h){a<-bh$hazard[as.character(bh$strata)==h & bh$time<=36];if(length(a))tail(a,1)else NA_real_},numeric(1));mean(exp(-hs*rr))})
 pred[[length(pred)+1]]<-data.table(outcome=o,landmark=lm,specification=spec,ASC=ps[1],AC_NOS=ps[2],difference=ps[1]-ps[2])
}
fwrite(rbindlist(diag),file.path(out,'TableS9_conditional_model_diagnostics.csv'));fwrite(rbindlist(pred),file.path(out,'TableS10_conditional_model_sensitivity.csv'))
# Baseline overlap weights evaluated again in later survivor populations.
ps<-glm(as.formula(paste('I(histology=="ASC")~',covars)),d,family=binomial());d[,ow:=ifelse(histology=='ASC',1-predict(ps,type='response'),predict(ps,type='response'))]
xx<-model.matrix(as.formula(paste('~',covars)),d)[,-1,drop=FALSE]; bal<-list();ess<-list()
for(lm in c(0,6,12,24)){
 keep<-d$months>lm;a<-keep & d$histology=='ASC';c<-keep & d$histology=='AC_NOS'
 for(j in seq_len(ncol(xx))){m1<-weighted.mean(xx[a,j],d$ow[a]);m0<-weighted.mean(xx[c,j],d$ow[c]);v1<-weighted.mean((xx[a,j]-m1)^2,d$ow[a]);v0<-weighted.mean((xx[c,j]-m0)^2,d$ow[c]);bal[[length(bal)+1]]<-data.table(landmark=lm,variable=colnames(xx)[j],SMD=abs(m1-m0)/sqrt((v1+v0)/2))}
 ess[[length(ess)+1]]<-d[keep,.(n=.N,ESS=sum(ow)^2/sum(ow^2),landmark=lm),by=histology]
}
fwrite(rbindlist(bal),file.path(out,'TableS11_weighted_survivor_balance.csv'));fwrite(rbindlist(ess),file.path(out,'TableS11b_weighted_survivor_ESS.csv'))
if(length(warnings))fwrite(rbindlist(warnings),file.path(out,'model_warning_audit.csv'))
cat('Audit completed\n');print(flow)

