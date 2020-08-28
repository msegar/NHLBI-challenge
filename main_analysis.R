# Import database
combined <- readRDS("datasets/combined.Rds")
combined <- droplevels(combined[combined$treatment != 'UF',]) # removes 94 from UF group (856 -> 762)


# ---------------------------- Additional Outcomes ----------------------------
combined$diureticEfficiency <- combined$output72h / combined$totalLasixIV
combined$GVASperD <- (combined$GLVASD - combined$GLVASB) / combined$LOS
combined$JVP <- as.factor(ifelse(combined$JVP == 3 | combined$JVP == 4, 1, 0)) # JVP >12 cm
combined$RALES <- as.factor(ifelse(combined$RALES == 2 | combined$RALES == 3, 1, 0)) # rales > 1/3
combined$ORTHPNEA <- as.factor(ifelse(combined$ORTHPNEA == 3 | combined$ORTHPNEA == 4, 1, 0)) # 3+ pillow orthopnea

is.na(combined$diureticEfficiency)<-sapply(combined$diureticEfficiency, is.infinite) # remove infinite values
combined[is.na(combined)] <- NA

# Filtering: Remove if missing output (removes 5, 762 -> 757)
combined <- combined[!is.na(combined$diureticEfficiency),]
combined <- combined[,c(1:39,46,53,60,67,74,80:105)]
d <- combined

# Remove colinear
library(corrplot)
x <- dplyr::select_if(d, is.numeric)
cormat <- cor(x, use = "complete.obs")
cormat[cormat > -0.7 & cormat < 0.7] = 0
par(xpd=TRUE)
corrplot(cormat, "square", "lower", order = "AOE", tl.col = "black", tl.cex = 0.7, mar = c(0,0,8,0))
rm(x, cormat)

# Remove cystatin C and GFR
d$CL_CYSTC <- NULL
d$CL_GFR <- NULL

# Impute
#visdat::vis_miss(d, sort_miss = T)
d <- d[, which(colMeans(is.na(d)) < 0.2)]
d <- randomForestSRC::impute(data=d, fast=T)
d.imp <- d

# ---------------------------- Analysis ----------------------------
# scale data
vars <- names(d[,c(2:5,8:45)])
d[,vars] <- scale(d[,vars])

# determine the top-20 most important variables
fmla.f <- as.formula(paste("diureticEfficiency", paste(vars, collapse = "+"), sep = "~"))
rf.mod <- rfsrc(fmla.f, d, importance = T, ntree = 1000, do.trace = T)
as.matrix(names(sort(rf.mod$importance, decreasing = T)))
top_vars <- names(sort(rf.mod$importance, decreasing = T))[1:20]

# unsupervised clustering of just the top-20 variables
temp <- d[,c(top_vars)]
res <- VarSelLCM::VarSelCluster(temp, gvals = 2:7) # evaluate between 2-7 clusters
d$cluster <- as.factor(VarSelLCM::fitted(res))
table(VarSelLCM::fitted(res))

# get most discriminative variabesl from VarSelLCM
names(sort(res@criteria@discrim, decreasing = T))[1:20]

# ---------------------------- Outcomes ----------------------------
derivation <- d[d$database == "rose" | d$database == "carress",] # setup derivation cohort
table(derivation$cluster)

# unadjusted outcomes
aovP(derivation$cluster, derivation$diureticEfficiency)
aovP(derivation$cluster, derivation$LOS)
aovP(derivation$cluster, derivation$GVASperD)
tablePCS(derivation$cluster, derivation$PTDIED)

# adjusted outcomes
lmCIG(lm(diureticEfficiency ~ cluster+AGE+FEMALE+RACE+cr0h+BPSYS+BPDIA+BMI+VLVEF,derivation),2)
lmCIG(lm(LOS ~ cluster+AGE+FEMALE+RACE+cr0h+BPSYS+BPDIA+BMI+VLVEF,derivation),2)
lmCIG(lm(GVASperD ~ cluster+AGE+FEMALE+RACE+cr0h+BPSYS+BPDIA+BMI+VLVEF,derivation),2)

# setup validation cohort
validation <- d[d$database == "dose" | d$database == 'escape',]
dose <- d[d$database == 'dose',]
escape <- d[d$database == 'escape',]
table(validation$cluster)

# unadjusted outcomes
aovP(validation$cluster, validation$diureticEfficiency)
aovP(validation$cluster, validation$LOS)
aovP(validation$cluster, validation$GVASperD)

# adjusted outcomes
lmCIG(lm(diureticEfficiency ~ cluster+AGE+FEMALE+RACE+cr0h+BPSYS+BPDIA+BMI+VLVEF,validation),2)
lmCIG(lm(LOS ~ cluster+AGE+FEMALE+RACE+cr0h+BPSYS+BPDIA+BMI+VLVEF,validation),2)
lmCIG(lm(GVASperD ~ cluster+AGE+FEMALE+RACE+cr0h+BPSYS+BPDIA+BMI+VLVEF,validation),2)

# escape long-term outcomes
escape <- d[d$database == 'escape',]
escape$DEIDNUM <- as.integer(as.character(escape$PATNUMB))
patient <- read.csv("../patient.csv") # read in the patient.csv file
escape <- merge(escape, subset(patient, select = c(DEIDNUM, T2DTHHSP, DTHRHSP)), by='DEIDNUM')
coxCIG(coxph(Surv(T2DTHHSP, DTHRHSP) ~ cluster, escape),2)
coxCIG(coxph(Surv(T2DTHHSP, DTHRHSP) ~ cluster+AGE+FEMALE+RACE+cr0h+BPSYS+BPDIA+BMI+VLVEF, escape),2)

# sensitivity analysis
aovP(derivation[derivation$VLVEF < 45,]$cluster, derivation[derivation$VLVEF < 45,]$diureticEfficiency) # HFrEF
aovP(derivation[derivation$VLVEF >= 45,]$cluster, derivation[derivation$VLVEF >= 45,]$diureticEfficiency) # HFpEF
aovP(derivation[derivation$FEMALE == 1,]$cluster, derivation[derivation$FEMALE == 1,]$diureticEfficiency) # female
aovP(derivation[derivation$FEMALE == 0,]$cluster, derivation[derivation$FEMALE == 0,]$diureticEfficiency) # male
aovP(derivation[derivation$RACE == 'black',]$cluster, derivation[derivation$RACE == 'black',]$diureticEfficiency) # black
aovP(derivation[derivation$RACE == 'white',]$cluster, derivation[derivation$RACE == 'white',]$diureticEfficiency) # black

# ---------------------------- Plots ----------------------------
library(ggpubr)
ggboxplot(derivation, x="cluster", y="diureticEfficiency", ylim = c(0,75))
ggboxplot(derivation, x="cluster", y="LOS", ylim = c(0,15))
ggboxplot(derivation, x="cluster", y="GVASperD", ylim = c(0,15))

library(tableone)
t1 <- CreateTableOne(vars=names(derivation), strata = "cluster", data=derivation[-1])
write.csv(print(t1, contDigits = 1, pDigits = 2, noSpaces = T),"baseline_characteristics.csv")

# ---------------------------- ESCAPE Additional Outcomes ----------------------------
hemo <- read.csv("../data/ESCAPE/data/main/analdata/hemo.csv")
escape <- combined[combined$database == "escape",]
escape$DEIDNUM <- as.integer(as.character(escape$PATNUMB))
escape <- merge(escape, hemo, by="DEIDNUM")

for (i in 1:length(escape)){ # convert to factor variables
  if(length(table(escape[,i])) < 8){
    escape[,i] <- as.factor(escape[,i])
  }
}
escape$PATNUMB <- NULL

library(tableone)
t1 <- CreateTableOne(vars=names(escape), strata = "cluster", data=escape[-1])
write.csv(print(t1, contDigits = 1, pDigits = 2, noSpaces = T),"escape_pa_cath.csv")

# echo data
echo <- read.csv("../data/ESCAPE/data/echo/analdata/echo_one.csv")
escape <- merge(escape, echo, by="DEIDNUM") #348 have echo data

t1 <- CreateTableOne(vars=names(escape), strata = "cluster", data=escape[-1])
write.csv(print(t1, contDigits = 1, pDigits = 2, noSpaces = T),"escape_echo.csv")

# evaluate improvement in C-index
library(rms)
cph1 <- cph(Surv(T2DTHHSP, DTHRHSP) ~ cluster+AGE+FEMALE, escape, x=T, y=T, surv = T, time.inc = 5)
cph2 <- cph(Surv(T2DTHHSP, DTHRHSP) ~ cluster+AGE+FEMALE+RACE+cr0h+BPSYS+BPDIA+BMI+VLVEF, escape, x=T, y=T, surv = T, time.inc = 5)
cph1 <- survest(cph1, escape, times = 5)
cph2 <- survest(cph2, escape, times = 5)
compareC::compareC(escape$T2DTHHSP, escape$DTHRHSP, cph1$surv, cph2$surv)

# print ROC curves
library(survivalROC)
cph1 <- survivalROC(Stime = (escape$T2DTHHSP+2), status = escape$DTHRHSP, marker = (1 - cph1$surv), predict.time = 5, method = 'KM')
cph2 <- survivalROC(Stime = (escape$T2DTHHSP+2), status = escape$DTHRHSP, marker = (1 - cph2$surv), predict.time = 5, method = 'KM')