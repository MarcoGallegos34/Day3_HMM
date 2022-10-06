library(readr)
library(dplyr)
library(momentuHMM)
library(ggplot2)
library(lubridate)

BlacktipB <- read_delim("C:/Users/marco/Downloads/BlacktipB_original.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
BlacktipB = BlacktipB %>% 
  mutate(Time = as.POSIXct(Time,format = "%m/%d/%Y %H:%M")) %>% 
  mutate(day = day(Time),hour_to_sec =  as.integer(seconds(hm(format(Time, format = "%H:%M")))))

BlacktipB = BlacktipB %>% group_by(Time) %>% mutate(sec = row_number()) %>% ungroup() %>% 
  mutate(sec = case_when(hour_to_sec == 53160 ~ as.integer(sec + 55),
                         TRUE ~ sec),
         hour_to_sec = hour_to_sec + sec,
         group_10s = case_when(sec <= 10 ~ 1,
                               10 < sec & sec <= 20 ~ 2,
                               20 < sec & sec <= 30 ~ 3,
                               30 < sec & sec <= 40~ 4,
                               40 < sec & sec <= 50 ~ 5,
                               50 < sec & sec <= 60 ~ 6,
                               TRUE ~ -1),
         group_30s = case_when(sec <= 30 ~ 1,
                               30 < sec & sec <= 60 ~ 2,
                               TRUE ~ -1)) %>% 
  group_by(Time,group_10s) %>% mutate(ODBA_mean_10s = mean(ODBA)) %>% ungroup() %>% 
  group_by(Time,group_30s) %>% mutate(ODBA_mean_30s = mean(ODBA)) %>% ungroup() %>% 
  group_by(Time) %>% mutate(ODBA_mean_60s = mean(ODBA)) %>%
  mutate(hour_to_sec_10s = floor((hour_to_sec-1)/10+1),
         hour_to_sec_30s = floor((hour_to_sec-1)/30+1),
         hour_to_sec_60s = floor((hour_to_sec-1)/60+1)) 

BlacktipB_10s = BlacktipB %>% group_by(Time,group_10s) %>% filter(row_number() == 1) %>%
  ungroup() %>% select(Time,Depth,ODBA_mean_10s,hour_to_sec_10s)
BlacktipB_30s = BlacktipB %>% group_by(Time,group_30s) %>% filter(row_number() == 1) %>% 
  ungroup() %>% select(Time,Depth,ODBA_mean_30s,hour_to_sec_30s)
BlacktipB_60s = BlacktipB %>% group_by(Time) %>% filter(row_number() == 1) %>% 
  ungroup() %>% select(Time,Depth,ODBA_mean_60s,hour_to_sec_60s)

BlacktipB = BlacktipB %>% select(Time,Depth,ODBA,hour_to_sec)

saveRDS(BlacktipB,"BlacktipB_tidy.rds")

plot.series.ODBA = BlacktipB %>% ggplot(aes(Time,ODBA)) + geom_line()

ggsave("./plots/timeSeriesODBA.png",plot.series.ODBA,
       width = 23.05,height = 15.05,units="cm",bg="white")


BlacktipB = BlacktipB %>% filter(ODBA <= 2.0)

plot.series.ODBA.constrained = BlacktipB %>% ggplot(aes(Time,ODBA)) + geom_line()

ggsave("./plots/timeSeriesODBA_constrained.png",plot.series.ODBA.constrained,
       width = 23.05,height = 15.05,units="cm",bg="white")


BlacktipBData = prepData(BlacktipB,coordNames = NULL,
                   covNames = "hour_to_sec")

BlacktipBData_10s = prepData(BlacktipB_10s,coordNames = NULL,
                   covNames = "hour_10s")

BlacktipBData_30s = prepData(BlacktipB_30s,coordNames = NULL,
                   covNames = "hour_30s")

BlacktipBData_60s = prepData(BlacktipB_60s,coordNames = NULL,
                   covNames = "hour_60s")

hist(BlacktipBData$ODBA,breaks=80)

fit1 = fitHMM(BlacktipBData,nbStates=2,dist=list(ODBA="gamma"),Par0 = list(ODBA=c(.1,.3,1,1)))
#fit1 = readRDS("BlacktipB_m1.rds")
saveRDS(fit1,"BlacktipB_m1.rds")


fit1_10s = fitHMM(BlacktipBData_10s,nbStates=2,dist=list(ODBA_mean_10s="gamma"),
                  Par0 = list(ODBA_mean_10s=c(.1,.3,1,1)))
fit1_10s
#fit1_10s = readRDS("BlacktipB_m1_10s.rds")
saveRDS(fit1_10s,"BlacktipB_m1_10s.rds")


fit1_30s = fitHMM(BlacktipBData_30s,nbStates=2,dist=list(ODBA_mean_30s="gamma"),
                  Par0 = list(ODBA_mean_30s=c(.1,.3,1,1)))
fit1_30s
#fit1_30s=readRDS("BlacktipB_m1_30s.rds")
saveRDS(fit1_30s,"BlacktipB_m1_30s.rds")

fit1_60s = fitHMM(BlacktipBData_60s,nbStates=2,dist=list(ODBA_mean_60s="gamma"),
                  Par0 = list(ODBA_mean_60s=c(.1,.3,1,1)))
fit1_60s
#fit1_60s = readRDS("BlacktipB_m1_60s.rds")
saveRDS(fit1_60s,"BlacktipB_m1_60s.rds")

BlacktipBData %>% ggplot(aes(x=ODBA)) + geom_histogram()



set.seed(147)
fit1_s2 <- fitHMM(BlacktipBData,
                  nbState = 2,
                  dist=list(ODBA="gamma"),
                  Par0 = list(ODBA=c(.1,.3,1,1)),
                  retryFits=10)

fit1_s2 <- fitHMM(BlacktipBData,
                  nbState = 2,
                  dist=list(ODBA="gamma"),
                  Par0 = list(ODBA=c(.1,1,1,1)))



#Checking for pseudo-residuals

### IMAGE GENERATION ###
jpeg("fit1.jpg", width = 700, height = 700)
plot(fit1,breaks = 80)
dev.off()


plotPR(fit1)

AIC(fit1)


### Fiting model with covariates

formula = ~ cosinor(hour_to_sec, period = 86400)
Par0_fit2 <- getPar0(model=fit1, formula=formula)

fit2 = fitHMM(BlacktipBData,nbStates=2,
              dist=list(ODBA="gamma"),Par0 = list(ODBA=c(.1,.3,1,1)),
              formula=formula)
fit2
fit2 = readRDS("BlacktipB_m2.rds")
saveRDS(fit2,"BlacktipB_m2.rds")

plot(fit2,breaks=80)

formula_10s = ~ cosinor(hour_to_sec_10s, period = 8640)
Par0_fit2_10s <- getPar0(model=fit1_10s, formula=formula_10s)

fit2_10s = fitHMM(BlacktipBData_10s,nbStates=2,dist=list(ODBA_mean_10s="gamma"),
                  Par0 = list(ODBA_mean_10s=c(.1,.3,1,1)),
                  formula=formula_10s)
fit2_10s
saveRDS(fit2_10s,"BlacktipB_m2_10s.rds")


formula_30s = ~ cosinor(hour_to_sec_30s, period = 2880)
Par0_fit2_30s <- getPar0(model=fit1_30s, formula=formula_30s)

fit2_30s = fitHMM(BlacktipBData_30s,nbStates=2,dist=list(ODBA_mean_30s="gamma"),
                  Par0 = list(ODBA_mean_30s=c(.1,.3,1,1)),
                  formula=formula_30s)
fit2_30s
saveRDS(fit2_30s,"BlacktipB_m2_30s.rds")


formula_60s = ~ cosinor(hour_to_sec_60s, period = 1440)
Par0_fit2_60s <- getPar0(model=fit1_60s, formula=formula_60s)

fit2_60s = fitHMM(BlacktipBData,nbStates=2,dist=list(ODBA_mean_60s="gamma"),
                  Par0 = list(ODBA_mean_60s=c(.1,.3,1,1)),
                  formula=formula_60s)
fit2_60s
fit2_60s = readRDS("BlacktipB_m2_60s.rds")
saveRDS(fit2_60s,"BlacktipB_m2_60s.rds")



library(readr)
library(ggplot2)
library(dplyr)

tigerShark <- read_csv("tigershark_depthchange10min.csv")
tigerShark = tigerShark %>% mutate(change_depth = Depth - lag(Depth,default = 0))


tigerShark %>%
  filter(days != 9) %>%
  ggplot(aes(x=HST,y=Depth)) + facet_wrap(~days,scales = "free_x") + geom_line() +
  scale_x_datetime(breaks= "8 hour", date_labels = "%H:%M") + theme_minimal()


tigerShark %>% filter(days != 9) %>% 
  ggplot(aes(x=HST,y=2*abs(change_depth))) + facet_wrap(~days,scales = "free_x") + geom_line() +
  scale_x_datetime(breaks= "8 hour", date_labels = "%H:%M") + theme_minimal()


tigerShark = tigerShark %>% select(HST,days,change_depth)

tigerShark$HST[1]

#load(url(paste0("https://static-content.springer.com/esm/",
#                "art%3A10.1007%2Fs13253-017-0282-9/MediaObjects/",
#                "13253_2017_282_MOESM2_ESM.rdata")))
#data[[1]]

#rbind(data.frame(HTS=tigerShark$HTS[1],level=c("1","2i")),data.frame(tigerShark,level="2"))


#test = data.frame(tigerShark %>% filter(days == 9) %>% filter(row_number() == 1) %>% select(HST),
#           level=c("1","2i"),
#           change_depth=NA)
#test

days_range = unique(tigerShark$days)
tigerSharkData <- NULL

for(i in days_range){
  coarseInd <- data.frame(tigerShark %>% filter(days == i) %>% 
                                   filter(row_number() == 1) %>% select(HST),
                                 days = i,
                                 level=c("1","2i"),
                                 change_depth=NA)
  tmp <- rbind(coarseInd,tigerShark %>% filter(days == i) %>% mutate(level = "2") %>% 
                 select(HST,days,level,change_depth))
  tigerSharkData <- rbind(tigerSharkData,tmp)
}

head(tigerSharkData)

tigerSharkData = prepData(tigerSharkData,
                          coordNames = NULL, hierLevels = c("1","2i","2"))


# summarize prepared data
summary(tigerSharkData, dataNames = names(tigerSharkData)[-1])


tigerSharkData

library(data.tree)

library(DiagrammeR)
### define hierarchical HMM
### states 1-3 = coarse state 1 (nonforaging)
### states 4-6 = coarse state 2 (foraging)
hierStates <- data.tree::Node$new("tiger shark HHMM states")
hierStates$AddChild(name="nontravelling")
hierStates$nontravelling$AddChild(name="nt1", state=1)
hierStates$nontravelling$AddChild(name="nt2", state=2)
hierStates$nontravelling$AddChild(name="nt3", state=3)
hierStates$AddChild(name="travelling")
hierStates$travelling$AddChild(name="t1", state=4)
hierStates$travelling$AddChild(name="t2", state=5)
hierStates$travelling$AddChild(name="t3", state=6)
plot(hierStates)

#Alternative way for specifying
hierStates <- data.tree::as.Node(list(name="tiger shark HHMM states",
                                      nontravelling=list(nt1=list(state=1),
                                                nt2=list(state=2),
                                                nt3=list(state=3)),
                                      travelling=list(t1=list(state=4),
                                                    t2=list(state=5),
                                                    t3=list(state=6))))

plot(hierStates)


# data stream distributions
# level 1 = coarse scale (no data streams)
# level 2 = fine scale (dive_duration, maximum_depth, dive_wiggliness)
hierDist <- data.tree::Node$new("tiger shark HHMM dist")
hierDist$AddChild(name="level1")
hierDist$AddChild(name="level2")
hierDist$level2$AddChild(name="change_depth", dist="gamma")
plot(hierDist)


# define hierarchical t.p.m. formula(s)
hierFormula <- data.tree::Node$new("tiger shark HHMM formula")
hierFormula$AddChild(name="level1", formula=~1)
hierFormula$AddChild(name="level2", formula=~1)

# define hierarchical initial distribution formula(s)
hierFormulaDelta <- data.tree::Node$new("tiger shark HHMM formulaDelta")
hierFormulaDelta$AddChild(name="level1", formulaDelta=~1)
hierFormulaDelta$AddChild(name="level2", formulaDelta=~1)

# defining starting values
cd.mu0 = rep(c(5,50,100),hierStates$count)
cd.sigma0 = rep(c(5,15,40),hierStates$count)

Par0 = list(change_depth = c(cd.mu0,cd.sigma0))

nbStates <- length(hierStates$Get("state",filterFun=data.tree::isLeaf))

# constrain fine-scale data stream distributions to be same
cd_DM = matrix(cbind(kronecker(c(1,1,0,0),diag(3)),
                     kronecker(c(0,0,1,1),diag(3))),
               nrow=nbStates*2,
               ncol=6,
               dimnames=list(c(paste0("mean_",1:nbStates),
                               paste0("sd_",1:nbStates)),
                             paste0(rep(c("mean","sd"),each=3),
                                    c("_14:(Intercept)",
                                      "_25:(Intercept)",
                                      "_36:(Intercept)"))))
    
DM = list(change_depth = cd_DM)


# get initial parameter values for data stream probability distributions
Par <- getParDM(tigerSharkData,hierStates=hierStates,hierDist=hierDist,
                Par=Par0,DM=DM)

# check hierarchical model specification and parameters
checkPar0(tigerSharkData,hierStates=hierStates,hierDist=hierDist,Par0=Par,
          hierFormula=hierFormula,hierFormulaDelta=hierFormulaDelta,
          DM=DM)

###

plotPR(fit1)
plotPR(fit1_10s)
plotPR(fit1_30s)
plotPR(fit1_60s)


BlacktipB %>% filter(ODBA >= .4)




load(url(paste0("https://static-content.springer.com/esm/",
                "art%3A10.1007%2Fs13253-017-0282-9/MediaObjects/",
                "13253_2017_282_MOESM2_ESM.rdata")))
# convert date_time to POSIX
data <- lapply(data,function(x)
{x$date_time <- as.POSIXct(x$date_time,tz="UTC"); x})
porpoiseData <- NULL
for(i in 1:length(data)){
  coarseInd <- data.frame(date_time=as.POSIXct(format(data[[i]]$date_time[1],
                                                      format="%Y-%m-%d %H:%M"),
                                               tz="UTC"),
                          level=c("1","2i"),
                          dive_duration=NA,
                          maximum_depth=NA,
                          dive_wiggliness=NA)
  tmp <- rbind(coarseInd,data.frame(data[[i]],level="2"))
  porpoiseData <- rbind(porpoiseData,tmp)
}

table(porpoiseData$dive_wiggliness)
min(porpoiseData$maximum_depth[!is.na(porpoiseData$maximum_depth)])

hist(porpoiseData$maximum_depth,breaks=80)
hist(porpoiseData$dive_wiggliness,breaks=80)
hist(abs(tigerSharkData$change_depth),breaks=80)

min(abs(tigerShark$change_depth)[!is.na(tigerShark$change_depth)])

table(abs(tigerShark$change_depth))

242/length(tigerShark$change_depth)

518/length(porpoiseData$dive_wiggliness)

dw_DM[1:(2*nbStates),1:6]



# fit hierarchical HMM
hhmm <- fitHMM(data=tigerSharkData,hierStates=hierStates,hierDist=hierDist,
               #hierFormula=hierFormula,#hierFormulaDelta=hierFormulaDelta,
               Par0=Par,#hierBeta=hierBeta,hierDelta=hierDelta,
               DM=DM,nlmPar=list(hessian=FALSE))


### TO DO ###
# Arreglar esto
# Error in getDM(tempCovs, inputs$DM, inputs$dist, nbStates, inputs$p$parNames,  : 
# DM$change_depth should consist of 18 rows