library(rcarbon)

binsize = 50
runm = 50
amazoncurve <- mixCurves('intcal20','shcal20',p=0.5)

dataset <- read.csv("m2m_dates.csv", header=TRUE)
dataset$SiteCode <- paste("S",as.numeric(as.factor(dataset$SiteName)),sep="")


calDates <- calibrate(x = dataset$Age, errors = dataset$Error, calCurves = amazoncurve, normalised = FALSE)
bins <- binPrep(sites = dataset$SiteCode, ages = dataset$Age, h = binsize) 

# Only Polychrome Tradition

polychrome <- subset(dataset, dataset$Affiliation == "Polychrome")
polychrome$SiteCode <- paste("S",as.numeric(as.factor(polychrome$SiteName)),sep="")

calDates.poly <- calibrate(x = polychrome$Age, errors = polychrome$Error, calCurves = amazoncurve)
bins.poly <- binPrep(sites = polychrome$SiteCode, ages = polychrome$Age, h = binsize)

# Visualisations

poly.spd <- spd(x = calDates.poly, bins = bins.poly, timeRange = c(1600,100), spdnormalised = FALSE)
smooth.spd <- spd(x = calDates.poly, bins = bins.poly, timeRange = c(1600,100), runm = runm, spdnormalised = FALSE)

s = sampleDates(calDates.poly, bins=bins.poly, nsim=1000, boot=TRUE)
ckdepoly = ckde(s, timeRange=c(1600,100), bw=25, normalised=FALSE)

save(poly.spd, smooth.spd, ckdepoly,
     file = "vis_prep.RData")

# Permutation test

ap <- subset(dataset, dataset$Affiliation == "Polychrome"|dataset$Affiliation=="Arauquinoid")
bins.ap <- binPrep(sites = ap$SiteCode, ages = ap$Age, h = binsize)
calDates.ap <- calibrate(x = ap$Age, errors = ap$Error, calCurves = amazoncurve, normalised = FALSE)

bp <- subset(dataset, dataset$Affiliation == "Polychrome"|dataset$Affiliation=="Barrancoid")
bins.bp <- binPrep(sites = bp$SiteCode, ages = bp$Age, h = binsize)
calDates.bp <- calibrate(x = bp$Age, errors = bp$Error, calCurves = amazoncurve, normalised = FALSE)

save(ap, bins.ap, calDates.ap,
     bp, bins.bp, calDates.bp,
     file = "perm_prep.RData")

# Bayesian MCMC setup

calDates.poly = calDates.poly[thinDates(ages=polychrome$Age,  errors=polychrome$Error, 
                              bins=bins.poly, size=1, thresh=1,seed=420,method='splitsample')]

index = which.CalDates(calDates.poly,BP<1601&BP>149,p=0.5)
calDates.poly = calDates.poly[index]

CRA=polychrome$Age[index]
Errors=polychrome$Error[index]
SiteID = polychrome$SiteCode[index]
LabCode = polychrome$LabCode[index]

medDates = medCal(calDates.poly)

if(any(medDates>1600|medDates<150)){medDates[medDates>1600]=1600;medDates[medDates<150]=150}

obs.data = data.frame(LabCode=LabCode,CRA=CRA,Error=Errors,MedCalDate=medDates,SiteID=SiteID)
ac.dat <- as.data.frame(amazoncurve)
constants <- list(N=length(calDates.poly),calBP=ac.dat$CALBP,C14BP=ac.dat$C14BP,C14err=ac.dat$Error)
data <- list(X=obs.data$CRA,sigma=obs.data$Error)

save(calDates.poly, obs.data, constants, ac.dat, data, medDates, file = "mcmc_prep.RData")
