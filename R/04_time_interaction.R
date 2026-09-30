suppressPackageStartupMessages({library(data.table);library(survival);library(splines)});setDTthreads(1)
z<-readRDS('work/revision_data.rds');d<-z$data;d[,time:=pmin(months,60)];d[,h:=as.integer(histology=='ASC')]
r<-list()
for(o in c('os_event','css_event')){
 d[,event:=get(o)*(months<=60)]
 q<-as.data.table(survSplit(Surv(time,event)~.,data=as.data.frame(d[,c('id','time','event','h','age','sex','race','period','stage','grade','subsite','sequence_group'),with=F]),cut=c(6,12,24),episode='interval'))
 q[,interval:=factor(interval)]
 f0<-coxph(as.formula(paste('Surv(tstart,time,event)~h+',z$covars,'+strata(interval)')),q,ties='efron')
 f1<-coxph(as.formula(paste('Surv(tstart,time,event)~h:interval+',z$covars,'+strata(interval)')),q,ties='efron')
 statistic<-2*(f1$loglik[2]-f0$loglik[2]);df<-length(coef(f1))-length(coef(f0))
 r[[length(r)+1]]<-data.table(outcome=o,LR=statistic,df=df,p=pchisq(statistic,df,lower.tail=F))
 print(r[[length(r)]])
}
fwrite(rbindlist(r),'outputs/final/supplementary/TableS12_histology_time_interaction.csv')
capture.output(sessionInfo(),file='outputs/final/supplementary/R_session_info.txt')
