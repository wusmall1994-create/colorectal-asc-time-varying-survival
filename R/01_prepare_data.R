suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(splines)
})

dir.create("work", showWarnings=FALSE, recursive=TRUE)
dir.create("outputs/final/supplementary", showWarnings=FALSE, recursive=TRUE)

asc_file <- Sys.getenv("SEER_ASC_FILE")
ac_file <- Sys.getenv("SEER_AC_NOS_FILE")
if (!nzchar(asc_file) || !nzchar(ac_file)) {
  stop("Set SEER_ASC_FILE and SEER_AC_NOS_FILE to the authorised local SEER exports.")
}

cols <- c(
  "Patient ID", "Sequence number", "Year of diagnosis", "Sex",
  "Race recode (W, B, AI, API)", "Age recode with single ages and 85+",
  "Primary Site - labeled", "Diagnostic Confirmation",
  "Combined Summary Stage with Expanded Regional Codes (2004+)",
  "SEER historic stage A (1973-2015)", "Survival months",
  "Vital status recode (study cutoff used)",
  "SEER cause-specific death classification", "Grade Recode (thru 2017)",
  "Grade Clinical (2018+)", "Grade Pathological (2018+)",
  "CS tumor size (2004-2015)", "Tumor Size Summary (2016+)",
  "Regional nodes examined (1988+)", "Regional nodes positive (1988+)",
  "RX Summ--Surg Prim Site (1998-2022)", "Radiation recode",
  "Chemotherapy recode (yes, no/unk)"
)

read_histology <- function(path, histology) {
  x <- fread(path, select=cols, colClasses=list(character="Patient ID"), showProgress=FALSE)
  setnames(x, c("id","sequence","year","sex","race_raw","age_label","site_raw",
    "confirmation","summary_stage","historic_stage","months_raw","vital","css_class",
    "grade_old","grade_clin","grade_path","size_old","size_new","nodes_examined_raw",
    "nodes_positive_raw","surgery_raw","radiation_raw","chemotherapy_raw"))
  x[, histology := histology]
  x
}

all <- rbindlist(list(read_histology(asc_file,"ASC"), read_histology(ac_file,"AC_NOS")), use.names=TRUE)
all[, age := as.numeric(sub("[^0-9].*$","",age_label))]
all <- all[age>=18 & year<=2018 & confirmation=="Positive histology" & grepl("^[0-9]+$",months_raw)]
all[, months := as.numeric(months_raw)][months==0, months:=0.5]
all[, sequence_order := suppressWarnings(as.integer(sub("^(\\d+)(st|nd|rd|th).*$","\\1",sequence)))]
all[sequence=="One primary only",sequence_order:=0L][is.na(sequence_order),sequence_order:=99L]
setorder(all,histology,id,year,sequence_order)
all <- all[, .SD[1], by=.(histology,id)]
pre_cross <- copy(all)
asc_ids <- all[histology=="ASC",unique(id)]
all <- all[!(histology=="AC_NOS" & id %chin% asc_ids)]

all[, `:=`(os_event=as.integer(vital=="Dead"),
  css_event=as.integer(css_class=="Dead (attributable to this cancer dx)"),
  cod_unknown=as.integer(css_class=="Dead (missing/unknown COD)"))]
all[, other_event := as.integer(vital=="Dead" & css_event==0 & cod_unknown==0)]
all[, stage := {z<-fifelse(year<2004,historic_stage,summary_stage); fifelse(grepl("^Localized",z),"Localized",fifelse(grepl("^Regional",z),"Regional",fifelse(grepl("^Distant",z),"Distant","Unknown")))}]
all[, grade := {old<-fifelse(grepl("^Well|^Moderately",grade_old),"I-II",fifelse(grepl("^Poorly|^Undifferentiated",grade_old),"III-IV","Unknown")); rc<-fifelse(grade_path %chin% as.character(1:4),grade_path,fifelse(grade_clin %chin% as.character(1:4),grade_clin,"9")); recent<-fifelse(rc %chin% c("1","2"),"I-II",fifelse(rc %chin% c("3","4"),"III-IV","Unknown")); fifelse(year<=2017,old,recent)}]
all[, subsite := {code<-sub("-.*$","",site_raw); fifelse(code %chin% c("C18.0","C18.2","C18.3","C18.4"),"Right colon",fifelse(code %chin% c("C18.5","C18.6","C18.7"),"Left colon",fifelse(code %chin% c("C19.9","C20.9"),"Rectosigmoid/rectum","Other colon")))}]
all[, race := fifelse(race_raw=="White","White",fifelse(race_raw=="Black","Black","Other"))]
all[, period := cut(year,c(1999,2006,2012,2018),labels=c("2000-2006","2007-2012","2013-2018"))]
all[, sequence_group := fifelse(sequence_order<=1,"Only/first","Later primary")]
all[, surgery := factor(fifelse(surgery_raw %chin% c("00","Blank(s)","99"),"No/unknown","Yes"))]
all[, radiation := factor(fifelse(grepl("None/Unknown|Recommended|Refused|Blank",radiation_raw),"No/unknown","Yes"))]
all[, chemotherapy := factor(fifelse(chemotherapy_raw=="Yes","Yes","No/unknown"))]
all[, histology := relevel(factor(histology),"AC_NOS")]
for(v in c("sex","race","period","stage","grade","subsite","sequence_group")) set(all,j=v,value=factor(all[[v]]))
covars <- "ns(age, 4) + sex + race + period + stage + grade + subsite + sequence_group"
saveRDS(list(data=all,pre_cross=pre_cross,covars=covars),"work/revision_data.rds")
rm(asc_ids); invisible(gc())
