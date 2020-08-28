########################################################################
######################## Read in ESCAPE data ###########################
########################################################################
patient <- read.csv("../data/ESCAPE/data/main/analdata/patient.csv")
labs <- read.csv("../data/ESCAPE/data/main/analdata/labs.csv")
mhx <- read.csv("../data/ESCAPE/data/main/analdata/mhist.csv")
volume <- read.csv("../data/ESCAPE/data/main/sasdata/volume.csv")
meds <- read.csv("../data/ESCAPE/data/main/sasdata/inptmed3.csv")
physical <- read.csv("../data/ESCAPE/data/main/sasdata/physexam.csv")
physical <- physical[physical$FORM == 'Baseline',]
quality <- read.csv("../data/ESCAPE/data/main/analdata/quality.csv")
demog <- read.csv("../data/ESCAPE/data/pac/analdata/demog.csv")
bios <- read.csv("../data/ESCAPE/data/biosite/analdata/bios.csv")
enzymes <- read.csv("../data/ESCAPE/data/main/sasdata/enzymes.csv")
enzymes <- enzymes[enzymes$FORM == 'Baseline',]
ecg <- read.csv("../data/ESCAPE/data/main/sasdata/ecg.csv")
ecg <- ecg[ecg$FORM == 'Baseline',]
visual <- read.csv("../data/ESCAPE/data/main/sasdata/visual.csv")
visual <- reshape2::dcast(visual, DEIDNUM~FORM, value.var = "GLOBAL", mean)
names(visual)[7:8] <- c("globalb", "globald")

volume <- volume[volume$VOLDT >= 0,]
volume$VOLNET <- ifelse(volume$VOLSIGN == 1, volume$VOLNET * (-1), volume$VOLNET)
volume <- reshape2::dcast(volume, DEIDNUM~VOLDT, value.var = "VOLNET",sum)
volume$urineOutputAvg <- apply(volume[,2:4],1,mean, na.rm=T)

########################################################################
####################### Merge the data #################################
########################################################################
meds <- meds[meds$MEDNM == 1 | meds$MEDNM == 3 | meds$MEDNM == 4,] # bumet=1, furos=3, tors=4
meds[] <- lapply(meds, function(x) ifelse(is.na(x), 0, x))
meds$lasixEquiv <- ifelse(meds$MEDNM == 1, (meds$MEDTOTPO * 40) + (meds$MEDTOTIV * 40), 
  ifelse(meds$MEDNM == 3, (meds$MEDTOTPO / 2) + (meds$MEDTOTIV), (meds$MEDTOTPO * 2) + (meds$MEDTOTIV * 2)))
meds <- aggregate(meds$lasixEquiv, list(meds$DEIDNUM), max)
names(meds) <- c("DEIDNUM","lasixEquiv")
temp <- merge(subset(volume, select = c(DEIDNUM, urineOutputAvg)), meds, by="DEIDNUM")
names(temp) <- c("DEIDNUM", 'output72h', 'totalLasixIV')
temp$diureticEfficiency <- temp$output72h / temp$totalLasixIV
temp$diureticEfficiency <- ifelse(temp$diureticEfficiency == 0, NA, temp$diureticEfficiency)
is.na(temp$diureticEfficiency)<-sapply(temp$diureticEfficiency, is.infinite)
temp <- temp[!is.na(temp$diureticEfficiency),] # 424 -> 409

patient$MALE <- ifelse(patient$GENDER == 1, 1, 0)
patient$FEMALE <- ifelse(patient$MALE == 1, 0, 1)
patient$RACE <- ifelse(patient$RACE == 1, 'white', ifelse(patient$RACE == 2, 'black', 'other'))
patient$treatment <- ifelse(patient$TX_NUM == 1, 'clinical assessment', 'PAC')
physical$JVP <- ifelse(physical$JVP == 0, NA, physical$JVP)
physical$HEPMEG <- ifelse(physical$HEPMEG > 0, 1, 0)
physical$ASCITES <- ifelse(physical$ASCITES > 0, 1, 0)
quality$ORTHOPB <- quality$ORTHOPB - 1
mhx$SMOKING <- mhx$SMOKING + 1
patient$HEIGHTIN <- patient$HGHTM * 39.3701
patient$WTLBS <- patient$WTADM * 2.2
enzymes$TROPVAL <- ifelse(enzymes$TROPTYP == 1, enzymes$TROPVAL * 100, NA)
ecg$RHYTHM <- ifelse(ecg$RHYTHM == 3, 2, ifelse(ecg$RHYTHM == 2, 4, ifelse(ecg$RHYTHM == 4, NA, ifelse(ecg$RHYTHM == 5, 3, ecg$RHYTHM))))
patient$PTDIED <- ifelse(patient$DEATH == 1 & patient$T2DEATH <= 120, 1, 0)

temp <- merge(temp, patient, by='DEIDNUM')
temp <- merge(temp, labs, by='DEIDNUM')
temp <- merge(temp, mhx, by='DEIDNUM')
temp <- merge(temp, physical, by='DEIDNUM')
temp <- merge(temp, quality, by='DEIDNUM')
temp <- merge(temp, enzymes, by='DEIDNUM')
temp <- merge(temp, ecg, by='DEIDNUM')
temp <- merge(temp, visual, by='DEIDNUM')
temp <- merge(temp, subset(bios, select = c(DEIDNUM, BNPB, BNPD, BNPM1)), by="DEIDNUM", all.x = T)
temp$CL_GFR <- nephro::CKDEpi.creat(temp$CRTB, temp$MALE, temp$AGE, ifelse(temp$RACE == 'black', 1, 0))

########################################################################
############## Match formatting of other ADHF datasets #################
########################################################################
escape <- data.frame(PATNUMB=temp$DEIDNUM, AGE=temp$AGE, FEMALE=temp$FEMALE, RACE=temp$RACE,
  YIHF=NA, VLVEF=temp$EF1, treatment=temp$treatment, database='escape', HRATE=temp$HRSUP1,
  BPSYS=temp$SYSBP1, BPDIA=temp$SDIABP1, SPO2=NA, JVP=temp$JVP, RALES=temp$RALES, AUSCULTN=temp$S3,
  HEPATOM=temp$HEPMEG, ASCITES=temp$ASCITES, NYHA=temp$NYHA1, ORTHPNEA=temp$ORTHOPB, CVHSP=NA, HFHSP=temp$NUMRHSP,
  ISCHEMIC=temp$ISCHD, ANGINA=temp$ANGP, MI=temp$MI, PTCI=temp$PTCI, CABG=temp$CABG, NONISCH=NA, HYPRTESN=temp$HYPERE,
  STROKE=temp$STROKE, ATRIALFB=temp$AFIB, ICD=temp$ICD, PVD=temp$PVD, COPD=temp$COPD, DIABETES=temp$DIAB, DEPRESS=temp$DEPR,
  ALCOHOL=temp$ALCHOE, CIGARETT=temp$SMOKING, LIPIDEMA=NA, na0h=temp$SODIUM1, na24h=NA, na48h=NA, na72h=temp$SODD3,
  na96h=NA, na7d=temp$SODD7, na60d=temp$SODM2, k0h=temp$POTB, k24h=NA, k48h=NA, k72h=temp$POTD3, k96h=NA, k7d=temp$POTD7,
  k60d=temp$POTM2, co0h=NA, co24h=NA, co48h=NA, co72h=NA, co96h=NA, co7d=NA, co60d=NA, bun0h=temp$BUN1, bun24h=NA, bun48h=NA,
  bun72h=temp$BUND3, bun96h=NA, bun7d=temp$BUND7, bun60d=temp$BUNM2, cr0h=temp$CREAT1, cr24h=NA, cr48h=NA, cr72h=temp$CRTD3,
  cr96h=NA, cr7d=temp$CRTD7, cr60d=temp$CRTM2, ntpro0h=temp$BNPB, ntpro72h=NA, ntpro7d=temp$BNPD, ntpro60d=temp$BNPM1,
  HEIGHTIN=temp$HEIGHTIN, WTLBS=temp$WTLBS, CL_GFR=temp$CL_GFR, CL_TROPI=temp$TROPVAL, LL_MAGNES=NA, LL_GLUC=NA, LL_TCHOL=NA,
  LL_AST=temp$ASTB, LL_ALT=temp$ALTB, LL_ALKPHOS=NA, LL_BILIR=NA, LL_ALBUM=temp$ALBB, LL_HMG=temp$HEMB, LL_WBC=temp$WBCB,
  LL_LYMPH=NA, LL_RDW=NA, CL_CYSTC=NA, BMI=temp$BMIB, ECGRHYTH=temp$RHYTHM, ECGQRS=NA, totalLasixIV=temp$totalLasixIV,
  output72h=temp$output72h, LOS=temp$HOSPDAY, PTDIED=temp$PTDIED, GLVASB=temp$globalb, GLVASD=temp$globald)

rm(patient, labs, mhx, meds, temp, volume, bios, demog, ecg, enzymes, physical, quality, visual)
#visdat::vis_miss(escape, sort_miss = T) # visualize missingness if desired
escape$PATNUMB <- as.factor(escape$PATNUMB)

# finally, combined all the ADHF datasets
d <- readRDS("datasets/dose_rose_carress.Rds")
combined <- rbind(d, escape)

saveRDS(combined, 'datasets/combined.Rds')
rm(list = ls())



