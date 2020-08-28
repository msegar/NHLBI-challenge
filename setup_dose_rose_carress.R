############################################################################
################################### BASE ###################################
############################################################################
dose_base <- haven::read_sas("../data/DOSE/data/analysis/a_base.sas7bdat")
rose_base <- haven::read_sas("../data/ROSE/Data/a_base.sas7bdat")
car_base <- haven::read_sas("../data/CARRESS/main_study/data/sas/a_base.sas7bdat")
dose_base$treatment <- as.factor(ifelse(dose_base$ROUTE == 0 & dose_base$INTENS == 0, "bolus_low",
                                        ifelse(dose_base$ROUTE == 0 & dose_base$INTENS == 1, "bolus_high",
                                               ifelse(dose_base$ROUTE == 1 & dose_base$INTENS == 1, "cont_high", "cont_log"))))
dose_base$FEMALE <- ifelse(dose_base$SEX == 2, 1, 0)
dose_base <- dose_base[,c(1,5,14,7,9,10,13)]
dose_base$database <- "dose"

rose_base$treatment <- as.factor(ifelse(rose_base$TREATMENT == 1, "placebo", ifelse(rose_base$TREATMENT == 2, "dopamine", "nesiritide")))
rose_base$RACE <- ifelse(rose_base$RACE == 99, NA, rose_base$RACE)
rose_base$FEMALE <- ifelse(rose_base$SEX == 2, 1, 0)
rose_base <- rose_base[,c(1,4,13,6,9,10, 12)]
rose_base$database <- "rose"

car_base$treatment <- as.factor(ifelse(car_base$TREATMENT == 1, "standard", "UF"))
car_base$FEMALE <- ifelse(car_base$SEX == 2, 1, 0)
car_base <- car_base[,c(1,3,16,5,7,8,15)]
car_base$database <- "carress"

base <- rbind(dose_base, rose_base, car_base)
base$RACE <- ifelse(base$RACE == 2, "white", ifelse(base$RACE == 1, "other", "black"))
base[,c(1,3,4,7,8)] <- lapply(base[,c(1,3,4,7,8)], as.factor)

saveRDS(base, "datasets/base.Rds")
rm(list = ls())

############################################################################
################################## MEDHX ###################################
############################################################################
dose_medhx <- haven::read_sas("../data/DOSE/data/sasdata/medhist1.sas7bdat")
rose_medhx <- haven::read_sas("../data/ROSE/Data/medhist1.sas7bdat")
car_medhx <- haven::read_sas("../data/CARRESS/main_study/data/sas/medhist1.sas7bdat")
vars <- names(dose_medhx)[names(dose_medhx) %in% names(rose_medhx)]
vars <- vars[c(2,3,7:9,17,19,21,32)]
dose_medhx <- dose_medhx[,vars]
rose_medhx <- rose_medhx[,vars]
car_medhx <- car_medhx[,vars]

dose_medhx2 <- haven::read_sas("../data/DOSE/data/sasdata/medhist2.sas7bdat")
rose_medhx2 <- haven::read_sas("../data/ROSE/Data/medhist2.sas7bdat")
car_medhx2 <- haven::read_sas("../data/CARRESS/main_study/data/sas/medhist2.sas7bdat")
vars <- names(dose_medhx2)[names(dose_medhx2) %in% names(rose_medhx2)]
vars <- vars[c(14,16,18,24,26:28,33:35,38,39)]
dose_medhx2 <- dose_medhx2[,vars]
rose_medhx2 <- rose_medhx2[,vars]
car_medhx2 <- car_medhx2[,vars]

mhx <- rbind(dose_medhx, rose_medhx, car_medhx)
mhx2 <- rbind(dose_medhx2, rose_medhx2, car_medhx2)
mhx <- merge(mhx, mhx2, by="patnumb")
names(mhx)[1] <- "PATNUMB"
mhx[is.na(mhx)] <- 0
mhx[,c(4:20)] <- lapply(mhx[,c(4:20)], as.factor)

saveRDS(mhx, "datasets/medhx.Rds")
rm(list = ls())

############################################################################
################################## OUTPUT ##################################
############################################################################
dose <- haven::read_sas("../data/DOSE/data/sasdata/fluid.sas7bdat")
rose <- haven::read_sas("../data/ROSE/Data/rsfluid.sas7bdat")
car <- haven::read_sas("../data/CARRESS/main_study/data/sas/crfluid.sas7bdat")

# dose
dose_output <- reshape2::dcast(dose, patnumb~FORM, value.var = "TOTALOUT")
dose_output$output72h <- apply(dose_output[,2:4], 1, sum, na.rm=T)
rose_output <- reshape2::dcast(rose, patnumb~FORM, value.var = "RSUROUT")
rose_output$output72h <- apply(rose_output[,2:4], 1, sum, na.rm=T)
car_output <- reshape2::dcast(car, patnumb~FORM, value.var = "CRURINOT")
car_output$output72h <- apply(car_output[,2:4],1,sum,na.rm=T)

output <- rbind(dose_output[,c(1,6)], rose_output[,c(1,5)], car_output[,c(1,9)])
names(output)[1] <- "PATNUMB"
saveRDS(output, "datasets/output.Rds")
rm(list = ls())

############################################################################
################################### ECG ####################################
############################################################################
dose <- haven::read_sas("../data/DOSE/data/sasdata/ecg.sas7bdat")
rose <- haven::read_sas("../data/ROSE/Data/ecg.sas7bdat")
car <- haven::read_sas("../data/CARRESS/main_study/data/sas/ecg.sas7bdat")

dose <- dose[,c(9,5,7)]
rose <- rose[,c(9,5,7)]
car <- car[,c(9,5,7)]

ecg <- rbind(dose, rose, car)
names(ecg)[1] <- "PATNUMB"
ecg$ECGRHYTH <- as.factor(ecg$ECGRHYTH)
saveRDS(ecg, "datasets/ecg.Rds")
rm(list = ls())

############################################################################
################################# DIURETIC #################################
############################################################################
dose <- haven::read_sas("../data/DOSE/data/analysis/a_endpts.sas7bdat")
rose <- haven::read_sas("../data/ROSE/Data/rsdaydiu.sas7bdat")
car <- haven::read_sas("../data/CARRESS/main_study/data/sas/crdaydiu.sas7bdat")

dose$totalLasixIV <- apply(dose[,c(73,75:77)], 1, sum, na.rm = T)
# rose and carress have same format
#------
rose$RSDAILY <- ifelse(rose$RSDAILY == 1, "lasix", ifelse(rose$RSDAILY == 2, "torsemide", ifelse(rose$RSDAILY == 3, "bumetanide", ifelse(rose$RSDAILY == 4, "metolazone", ifelse(rose$RSDAILY == 5, "hctz", "chlorothiazide")))))
rose[is.na(rose)] <- 0
names(rose)[4] <- "med"
rose$totalIV <- rose$RSBOLVAL + rose$RSIVVAL
rose$totalIV <- ifelse(rose$med == "lasix", rose$totalIV, 
                       ifelse(rose$med == "torsemide", (rose$totalIV * 2),
                              ifelse(rose$med == "bumetanide", (rose$totalIV * 40), 0)))
rose$totalPO <- ifelse(rose$med == "lasix", (rose$RSPOVAL / 2), 
                       ifelse(rose$med == "torsemide", (rose$RSPOVAL * 2),
                              ifelse(rose$med == "bumetanide", (rose$RSPOVAL * 40), 0)))
rose$totalLasix <- rose$totalIV + rose$totalPO
time2 <- rose[rose$RSTIMEPT == 2,]
time2 <- aggregate(time2$totalLasix, list(time2$patnumb), sum)
time3 <- rose[rose$RSTIMEPT == 3,]
time3 <- aggregate(time3$totalLasix, list(time3$patnumb), sum)
time4 <- rose[rose$RSTIMEPT == 4,]
time4 <- aggregate(time4$totalLasix, list(time4$patnumb), sum)
time5 <- rose[rose$RSTIMEPT == 5,]
time5 <- aggregate(time5$totalLasix, list(time5$patnumb), sum)
rose_sums <- data.frame(PATNUMB = time2$Group.1, totalLasixIV = (time2$x + time3$x + time4$x + time5$x))
#------
car$CRDAILY <- ifelse(car$CRDAILY == 1, "lasix", ifelse(car$CRDAILY == 2, "toCRemide", ifelse(car$CRDAILY == 3, "bumetanide", ifelse(car$CRDAILY == 4, "metolazone", ifelse(car$CRDAILY == 5, "hctz", "chlorothiazide")))))
car[is.na(car)] <- 0
names(car)[4] <- "med"
car$totalIV <- car$CRBOLVAL + car$CRIVVAL
car$totalIV <- ifelse(car$med == "lasix", car$totalIV, 
                       ifelse(car$med == "toCRemide", (car$totalIV * 2),
                              ifelse(car$med == "bumetanide", (car$totalIV * 40), 0)))
car$totalPO <- ifelse(car$med == "lasix", (car$CRPOVAL / 2), 
                       ifelse(car$med == "toCRemide", (car$CRPOVAL * 2),
                              ifelse(car$med == "bumetanide", (car$CRPOVAL * 40), 0)))
car$totalLasix <- car$totalIV + car$totalPO
time2 <- car[car$CRTIMEPT == 2,]
time2 <- aggregate(time2$totalLasix, list(time2$patnumb), sum)
time3 <- car[car$CRTIMEPT == 3,]
time3 <- aggregate(time3$totalLasix, list(time3$patnumb), sum)
time4 <- car[car$CRTIMEPT == 4,]
time4 <- aggregate(time4$totalLasix, list(time4$patnumb), sum)
time5 <- car[car$CRTIMEPT == 5,]
time5 <- aggregate(time5$totalLasix, list(time5$patnumb), sum)
car_sums <- data.frame(PATNUMB = time2$Group.1, totalLasixIV = (time2$x + time3$x + time4$x + time5$x))

diuretics <- rbind(dose[,c(1,78)], rose_sums, car_sums)
saveRDS(diuretics, "datasets/diuretics.Rds")
rm(list = ls())

############################################################################
################################ ASSESSMENT ################################
############################################################################
dose <- haven::read_sas("../data/DOSE/data/sasdata/assessmt.sas7bdat")
rose <- haven::read_sas("../data/ROSE/Data/assessmt.sas7bdat")
car <- haven::read_sas("../data/CARRESS/main_study/data/sas/assessmt.sas7bdat")
vars <- names(dose)[names(dose) %in% names(rose)]
vars <- vars[c(3,5,6,8,12,14,16,18,20,23,25,26)]

dose <- dose[dose$FORM == "BASELINE",]
rose <- rose[rose$FORM == "BASELINE",]
car <- car[car$FORM == "BASELINE",]
dose <- dose[,vars]
rose <- rose[,vars]
car <- car[,vars]

assessment <- rbind(dose, rose, car)
assessment <- assessment[,c(12,1:11)]
assessment[,c(6:12)] <- lapply(assessment[,c(6:12)], as.factor)

names(assessment)[1] <- "PATNUMB"
saveRDS(assessment, "datasets/assessment.Rds")
rm(list = ls())

############################################################################
################################### LABS ###################################
############################################################################
dose_labs <- haven::read_sas("../data/DOSE/data/analysis/a_visitsumm.sas7bdat")
rose_labs <- haven::read_sas("../data/ROSE/Data/a_visitsumm.sas7bdat")
car_labs <- haven::read_sas("../data/CARRESS/main_study/data/sas/a_visitsumm.sas7bdat")

# creatinine ------------------------------------
dose_labs$screat <- ifelse(is.na(dose_labs$CL_CREAT), dose_labs$LL_CREAT, dose_labs$CL_CREAT)
dose_cr <- reshape2::dcast(dose_labs, PATNUMB~FORM, value.var = "screat")
names(dose_cr) <- c("PATNUMB", "cr24h", "cr48h", "cr72h", "cr96h", "cr0h", "cr60d", "cr7d")
dose_cr <- dose_cr[,c(1,6,2:5,8,7)]
rose_labs$screat <- ifelse(is.na(rose_labs$CL_CREAT), rose_labs$LL_CREAT, rose_labs$CL_CREAT)
rose_cr <- reshape2::dcast(rose_labs, PATNUMB~FORM, value.var = "screat")
names(rose_cr) <- c("PATNUMB", "cr24h", "cr48h", "cr72h", "cr0h","cr7d")
rose_cr$cr96h <- NA
rose_cr$cr60d <- NA
rose_cr <- rose_cr[,c(1,5,2:4,7,6,8)]
car_labs$screat <- ifelse(is.na(car_labs$CL_CREAT), car_labs$LL_CREAT, car_labs$CL_CREAT)
car_cr <- reshape2::dcast(car_labs, PATNUMB~FORM, value.var = "screat")
car_cr <- car_cr[,c(1:5,7,11,10)]
names(car_cr) <- c("PATNUMB","cr0h", "cr24h", "cr48h", "cr72h", "cr96h","cr7d","cr60d")

# bun ------------------------------------
dose_bun <- reshape2::dcast(dose_labs, PATNUMB~FORM, value.var = "LL_BUN")
names(dose_bun) <- c("PATNUMB", "bun24h", "bun48h", "bun72h", "bun96h", "bun0h", "bun60d", "bun7d")
dose_bun <- dose_bun[,c(1,6,2:5,8,7)]
rose_bun <- reshape2::dcast(rose_labs, PATNUMB~FORM, value.var = "LL_BUN")
names(rose_bun) <- c("PATNUMB", "bun24h", "bun48h", "bun72h", "bun0h","bun7d")
rose_bun$bun96h <- NA
rose_bun$bun60d <- NA
rose_bun <- rose_bun[,c(1,5,2:4,7,6,8)]
car_bun <- reshape2::dcast(car_labs, PATNUMB~FORM, value.var = "LL_BUN")
car_bun <- car_bun[,c(1:5,7,11,10)]
names(car_bun) <- names(rose_bun)

# sodium ------------------------------------
dose_na <- reshape2::dcast(dose_labs, PATNUMB~FORM, value.var = "LL_SODIUM")
names(dose_na) <- c("PATNUMB", "na24h", "na48h", "na72h", "na96h", "na0h", "na60d", "na7d")
dose_na <- dose_na[,c(1,6,2:5,8,7)]
rose_na <- reshape2::dcast(rose_labs, PATNUMB~FORM, value.var = "LL_SODIUM")
names(rose_na) <- c("PATNUMB", "na24h", "na48h", "na72h", "na0h","na7d")
rose_na$na96h <- NA
rose_na$na60d <- NA
rose_na <- rose_na[,c(1,5,2:4,7,6,8)]
car_na <- reshape2::dcast(car_labs, PATNUMB~FORM, value.var = "LL_SODIUM")
car_na <- car_na[,c(1:5,7,11,10)]
names(car_na) <- names(rose_na)

# potassium ------------------------------------
dose_k <- reshape2::dcast(dose_labs, PATNUMB~FORM, value.var = "LL_POTASS")
names(dose_k) <- c("PATNUMB", "k24h", "k48h", "k72h", "k96h", "k0h", "k60d", "k7d")
dose_k <- dose_k[,c(1,6,2:5,8,7)]
rose_k <- reshape2::dcast(rose_labs, PATNUMB~FORM, value.var = "LL_POTASS")
names(rose_k) <- c("PATNUMB", "k24h", "k48h", "k72h", "k0h","k7d")
rose_k$k96h <- NA
rose_k$k60d <- NA
rose_k <- rose_k[,c(1,5,2:4,7,6,8)]
car_k <- reshape2::dcast(car_labs, PATNUMB~FORM, value.var = "LL_POTASS")
car_k <- car_k[,c(1:5,7,11,10)]
names(car_k) <- names(rose_k)

# bicarb ------------------------------------
dose_co <- reshape2::dcast(dose_labs, PATNUMB~FORM, value.var = "LL_BICARB")
names(dose_co) <- c("PATNUMB", "co24h", "co48h", "co72h", "co96h", "co0h", "co60d", "co7d")
dose_co <- dose_co[,c(1,6,2:5,8,7)]
rose_co <- reshape2::dcast(rose_labs, PATNUMB~FORM, value.var = "LL_BICARB")
names(rose_co) <- c("PATNUMB", "co24h", "co48h", "co72h", "co0h","co7d")
rose_co$co96h <- NA
rose_co$co60d <- NA
rose_co <- rose_co[,c(1,5,2:4,7,6,8)]
car_co <- reshape2::dcast(car_labs, PATNUMB~FORM, value.var = "LL_BICARB")
car_co <- car_co[,c(1:5,7,11,10)]
names(car_co) <- names(rose_co)

# NTprob ------------------------------------
dose_ntpro <- reshape2::dcast(dose_labs, PATNUMB~FORM, value.var = "CL_NTPRO")
names(dose_ntpro) <- c("PATNUMB", "ntpro24h", "ntpro48h", "ntpro72h", "ntpro96h", "ntpro0h", "ntpro60d", "ntpro7d")
dose_ntpro <- dose_ntpro[,c(1,6,2:5,8,7)]
rose_ntpro <- reshape2::dcast(rose_labs, PATNUMB~FORM, value.var = "CL_NTPRO")
names(rose_ntpro) <- c("PATNUMB", "ntpro24h", "ntpro48h", "ntpro72h", "ntpro0h","ntpro7d")
rose_ntpro$ntpro96h <- NA
rose_ntpro$ntpro60d <- NA
rose_ntpro <- rose_ntpro[,c(1,5,2:4,7,6,8)]
car_ntpro <- reshape2::dcast(car_labs, PATNUMB~FORM, value.var = "CL_NTPRO")
car_ntpro <- car_ntpro[,c(1:5,7,11,10)]
names(car_ntpro) <- names(rose_ntpro)

# Other labs ------------------------------------
dose_baseline <- dose_labs[dose_labs$FORM == "BASELINE",]
rose_baseline <- rose_labs[rose_labs$FORM == "BASELINE",]
car_baseline <- car_labs[car_labs$FORM == "BASELINE",]

dose_baseline <- dose_baseline[,c(1,3,4,33,32,11:22,28)]
rose_baseline <- rose_baseline[,c(1,3,4,33,27,11:22,31)]
rose_baseline$LL_TPROT <- NA # placehold for troponin since it's missing in ROSE
car_baseline <- car_baseline[,c(1,3,4,30,37,11:22,33)]

names(rose_baseline) <- names(dose_baseline)
names(car_baseline) <- names(dose_baseline)

# Merge --------------------------------------------
creatinine <- rbind(dose_cr, rose_cr, car_cr)
bun <- rbind(dose_bun, rose_bun, car_bun)
sodium <- rbind(dose_na, rose_na, car_na)
potass <- rbind(dose_k, rose_k, car_k)
bicarb <- rbind(dose_co, rose_co, car_co)
ntpro <- rbind(dose_ntpro, rose_ntpro, car_ntpro)
ntpro[,c("ntpro24h","ntpro48h","ntpro96h")] <- NULL
other <- rbind(dose_baseline, rose_baseline, car_baseline)
other$BMI <- (other$WTLBS) / (other$HEIGHTIN^2) * 703

labs <- merge(sodium, potass, by="PATNUMB")
labs <- merge(labs, bicarb, by="PATNUMB")
labs <- merge(labs, bun, by="PATNUMB")
labs <- merge(labs, creatinine, by="PATNUMB")
labs <- merge(labs, ntpro, by="PATNUMB")
labs <- merge(labs, other, by="PATNUMB")

labs$PATNUMB <- as.factor(labs$PATNUMB)
saveRDS(labs, "datasets/labs.Rds")
rm(list = ls())

############################################################################
#################### LENGTH OF STAY AND MORTALITY ##########################
############################################################################
dose_endpts <- haven::read_sas("../data/DOSE/data/analysis/a_endpts.sas7bdat")
rose_endpts <- haven::read_sas("../data/ROSE/Data/a_endpts.sas7bdat")
car_endpts <- haven::read_sas("../data/CARRESS/main_study/data/sas/a_endpts.sas7bdat")
rose_endpts$PTDIED <- ifelse(rose_endpts$PTDIED6MO == 1 & rose_endpts$TIMETODTH6MO <= 120, 1, 0)

dose <- subset(dose_endpts, select = c(PATNUMB, DRNDIS, PTDIED))
rose <- subset(rose_endpts, select = c(PATNUMB, DRNDIS, PTDIED))
car <- subset(car_endpts, select = c(PATNUMB, DRNDIS, PTDIED))
los <- rbind(dose, rose, car)
names(los) <- c("PATNUMB", "LOS", "PTDIED")

los$PATNUMB <- as.factor(los$PATNUMB)
saveRDS(los, "datasets/los.Rds")
rm(list = ls())

############################################################################
############################### Global VAS #################################
############################################################################
dose_vas <- haven::read_sas("../data/DOSE/data/sasdata/subjsymp.sas7bdat")
rose_vas <- haven::read_sas("../data/ROSE/Data/vas.sas7bdat")
car_vas <- haven::read_sas("../data/CARRESS/main_study/data/sas/vas.sas7bdat")
dose_vas <- reshape2::dcast(dose_vas, patnumb~SELFDT, value.var = 'VASGLOBL', mean, na.rm=T)
rose_vas <- reshape2::dcast(rose_vas, patnumb~FORM, value.var = "GLOBLVAS", mean, na.rm=T)
car_vas <- reshape2::dcast(car_vas, patnumb~FORM, value.var = "GLOBLVAS", mean, na.rm=T)

dose_vas <- dose_vas[,c(1,2,5)]
rose_vas <- rose_vas[,c(1,5,4)]
car_vas <- car_vas[,c(1,2,6)]

names(dose_vas) <- names(rose_vas)
names(car_vas) <- names(rose_vas)
vas <- rbind(dose_vas, rose_vas, car_vas)
names(vas) <- c('PATNUMB', "GLVASB", 'GLVASD') 
vas$PATNUMB <- as.factor(vas$PATNUMB)
saveRDS(vas, "datasets/vas.Rds")
rm(list = ls())

############################################################################
################################## COMBINE #################################
############################################################################
base <- readRDS("datasets/base.Rds")
assessment <- readRDS("datasets/assessment.Rds")
medhx <- readRDS("datasets/medhx.Rds")
diuretics <- readRDS("datasets/diuretics.Rds")
output <- readRDS("datasets/output.Rds")
labs <- readRDS("datasets/labs.Rds")
ecg <- readRDS("datasets/ecg.Rds")
los <- readRDS("datasets/los.Rds")
vas <- readRDS("datasets/vas.Rds")

d <- merge(base, assessment, by="PATNUMB")
d <- merge(d, medhx, by="PATNUMB")
d <- merge(d, labs, by="PATNUMB")
d <- merge(d, ecg, by="PATNUMB")
d <- merge(d, diuretics, by="PATNUMB")
d <- merge(d, output, by="PATNUMB")
d <- merge(d, los, by="PATNUMB")
d <- merge(d, vas, by="PATNUMB")

saveRDS(d, "datasets/dose_rose_carress.Rds")
rm(list = ls())
