##########################################
# Analysis of daily on and off bouts
# 2012, 2014 & 2015 data
# K. Williams
# July 2022
# for ms: Consequences of nest microclimate on incubation behavior 
##########################################
setwd("D:/Dropbox/kelly's shared files/Nest Temp research HOWA")
setwd("/Users/williak5/Dropbox/kelly's shared files/Nest Temp research HOWA")
# used offbout.loops.R to pull daily on and off bouts
# saved 3 data frames df_bouts_2012 , _2014 & _2015
# off/on  bouts defined by 1 degree change in T between readings. Missed some,
# plots saved in plots folder. Use these to determine which nests to drop
# see Data.prep.R for code and details


#### load libraries ####
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(psych)
library(emmeans)
library(nlme)
library(MuMIn)
library(viridis)
library(gridExtra)
library(lme4)
library(car)
library(jtools)
library(glmmTMB)
library(iccCounts)


#### start data exploration & analyses ####
bout_dat<-read.csv("bout_dat.csv")
length(unique(levels(factor(bout_dat$Nest.ID)))) # 59 nests
levels(factor(bout_dat$Type))
levels(factor(bout_dat$Stage))
levels(factor(bout_dat$Day.night))
levels(factor(bout_dat$aspect))
bout_dat$T.off<-bout_dat$mean.time.off.bout/60 # mean time off in minutes
bout_dat$T.on<-bout_dat$mean.time.on.bout/60 # mean time on in minutes
levels((factor(bout_dat$F.Age)))

####compass directions for aspect var name = f.aspect. dont need if read in temp_sub above vs start from scratch####
df.cut<-cut(bout_dat$aspect, breaks=c(0,90,180,270,360), labels=c("NE", "SE","SW","NW"))
#levels(factor(df.cut))
bout_dat$f.aspect<-factor(df.cut)
levels(factor(bout_dat$f.aspect))
names(bout_dat)

#### adjust on and off bouts by prop.day.r (proportion of each day the nest was observed rounded to 0.1 to adjust for time intervals of 4 min. range 0-1)
summary(bout_dat$prop.day.r)
which(bout_dat$prop.day.r == 0) # 454 455 484 485 these had 0off bouts 1 on and no real time. 
bout_dat$on.per.day<-bout_dat$number.on.bouts/bout_dat$prop.day.r
summary(bout_dat$on.per.day)
bout_dat$on.per.day[which(bout_dat$on.per.day == "Inf")]<- ""
bout_dat$on.per.day<-as.numeric(bout_dat$on.per.day)
summary(bout_dat$on.per.day) # no longer count data
hist(bout_dat$on.per.day) #normal 
which(bout_dat$on.per.day > 50)
bout_dat$off.per.day<-bout_dat$number.off.bouts/bout_dat$prop.day.r
summary(bout_dat$off.per.day) # no longer count data
hist(bout_dat$off.per.day) # normal
which(bout_dat$off.per.day > 50)

ggplot(bout_dat, aes(Day.night, off.per.day))+
       geom_boxplot()
       
#### Does TPI, aspect or amb temp affect incubation behavior ####

####subset data by Type = ambient ####
bout_amb_dat<-bout_dat %>% filter(Type == "ambient") # 391 obs
bout_amb_dat.red<-bout_amb_dat %>% filter(!Nest.ID == c("321.12")) # nest has low perc.in and weird values. weird, real but outlier.
levels(factor(bout_amb_dat.red$Nest.ID))
levels(factor(bout_amb_dat$Nest.ID))
names(bout_amb_dat)


#### order TPI ####
bout_amb_dat$tpi50_6cat<-factor(bout_amb_dat$tpi50_6cat, levels = c("ridge", "upper slope", "steep slope", "gentle slope", "lower slope", "valley") )
bout_amb_dat.red$tpi50_6cat<-factor(bout_amb_dat.red$tpi50_6cat, levels = c("ridge", "upper slope", "steep slope", "gentle slope", "lower slope", "valley") )

# explore data
hist(bout_amb_dat$number.on.bouts)
which(bout_amb_dat$number.on.bouts >= 60)
which(bout_amb_dat$number.on.bouts >= 40)
which(bout_amb_dat$number.on.bouts <= 10)
on.10<-bout_amb_dat%>%filter(number.on.bouts <= 10)
levels(factor(on.10$Nest.ID))
table(on.10$J.day, on.10$Nest.ID)
summary(bout_amb_dat$number.on.bouts)
describeBy(bout_amb_dat$number.on.bouts)
hist(bout_amb_dat$number.off.bouts)
describeBy((bout_amb_dat$number.off.bouts))
# note need to use adjusted on and off bouts instead of number
describeBy(bout_dat$on.per.day)
describeBy(bout_dat$off.per.day)
summary(bout_amb_dat$mean.time.off.bout) # skewed longest 7200s = 120min, mean = 540s = 9 min
summary(bout_amb_dat$mean.time.on.bout) #skewed longest on = 28240 = 470.67min; mean = 3146s = 52.43 min
summary(bout_amb_dat$T.on)
summary(bout_amb_dat$T.off)

#descriptive stats TPI 6 & 3
describeBy(bout_amb_dat$number.on.bouts, bout_amb_dat$tpi50_6cat, mat = TRUE)
describeBy(bout_amb_dat$number.on.bouts, bout_amb_dat$tpi50_3cat, mat = TRUE)
describeBy(bout_amb_dat$number.off.bouts, bout_amb_dat$tpi50_6cat, mat = TRUE)
describeBy(bout_amb_dat$number.off.bouts, bout_amb_dat$tpi50_3cat, mat = TRUE)
describeBy(bout_amb_dat$on.per.day, bout_amb_dat$tpi50_6cat, mat = TRUE)
describeBy(bout_amb_dat$on.per.day, bout_amb_dat$tpi50_3cat, mat = TRUE)
describeBy(bout_amb_dat$off.per.day, bout_amb_dat$tpi50_6cat, mat = TRUE)
describeBy(bout_amb_dat$off.per.day, bout_amb_dat$tpi50_3cat, mat = TRUE)
describeBy(bout_amb_dat$T.on, bout_amb_dat$tpi50_6cat, mat = TRUE)
describeBy(bout_amb_dat$T.on, bout_amb_dat$tpi50_3cat, mat = TRUE)
describeBy(bout_amb_dat$T.off, bout_amb_dat$tpi50_6cat, mat = TRUE)
describeBy(bout_amb_dat$T.off, bout_amb_dat$tpi50_3cat, mat = TRUE)
describeBy(bout_amb_dat$perc.in, bout_amb_dat$tpi50_6cat, mat = TRUE)
describeBy(bout_amb_dat$perc.in, bout_amb_dat$tpi50_3cat, mat = TRUE)

#descriptive stats f.aspect
describeBy(bout_amb_dat$number.on.bouts, bout_amb_dat$f.aspect, mat = TRUE)
describeBy(bout_amb_dat$number.off.bouts, bout_amb_dat$f.aspect, mat = TRUE)
describeBy(bout_amb_dat$on.per.day, bout_amb_dat$f.aspect, mat = TRUE)
describeBy(bout_amb_dat$off.per.day, bout_amb_dat$f.aspect, mat = TRUE)
describeBy(bout_amb_dat$T.on, bout_amb_dat$f.aspect, mat = TRUE)
describeBy(bout_amb_dat$T.off, bout_amb_dat$f.aspect, mat = TRUE)
describeBy(bout_amb_dat$perc.in, bout_amb_dat$f.aspect, mat = TRUE)

# scale J.day
range(bout_amb_dat$J.day)
bout_amb_dat$jday.s<-scale(bout_amb_dat$J.day) # scale Julian day
summary(bout_amb_dat$jday.s)
bout_amb_dat.red$jday.s<-scale(bout_amb_dat.red$J.day) # scale Julian day

# use on and off per day instead of number below.
N.on.jday.lme<-glmer(number.on.bouts ~  jday.s+ (1|Nest.ID), family = poisson, data = bout_amb_dat,  na.action=na.pass)
summary(N.on.jday.lme)
opd.lm.day<-lme(on.per.day~jday.s, random = ~1|Nest.ID, data= bout_amb_dat, na.action = na.omit)
summary(opd.lm.day) # no diff in on bouts by jday
plot(opd.lm.day)

offpd.lm.day<-lme(off.per.day~jday.s, random = ~1|Nest.ID, data = bout_amb_dat, na.action = na.omit)
summary(offpd.lm.day)

N.on.tpi.lme4<-glmer(number.on.bouts ~  F.Age + (1|Nest.ID), family = poisson, data = bout_amb_dat,  na.action=na.pass)
summary(N.on.tpi.lme4) # no diff by age
opd.lm.age<-lme(on.per.day~F.Age, random = ~1|Nest.ID, data = bout_amb_dat, na.action = na.omit)
summary(opd.lm.age)
plot(opd.lm.age)
pairs(emmeans(opd.lm.age, ~F.Age))

# since Age nor Jday sig. don't need to put in models of # on bouts

#### Does TPI3, affect # on bouts? ####
# NB number of on and off bouts not corrected. Redo all analy with on. or off.per.day
N.on.tpi.lme<-glmer(number.on.bouts ~ jday.s + tpi50_3cat + F.Age + (1|Nest.ID), family = poisson, data = bout_amb_dat,  na.action=na.pass)
N.on.tpi.lme1<-glmer(number.on.bouts ~ jday.s + tpi50_3cat + (1|Nest.ID), family = poisson, data = bout_amb_dat,  na.action=na.pass)
N.on.tpi.lme2<-glmer(number.on.bouts ~  tpi50_3cat + F.Age + (1|Nest.ID), family = poisson, data = bout_amb_dat,  na.action=na.pass)
N.on.tpi.lme3<-glmer(number.on.bouts ~  tpi50_3cat  + (1|Nest.ID), family = poisson, data = bout_amb_dat,  na.action=na.pass)

model.sel(N.on.tpi.lme, N.on.tpi.lme1, N.on.tpi.lme2, N.on.tpi.lme3, N.on.tpi.lme4)
summary(N.on.tpi.lme4)
pairs(emmeans(N.on.tpi.lme4, ~F.Age, method="tukey"))

#### TPI 3 best #### not different when drop 321.12
summary(N.on.tpi.lme3)
N.on.TPI3.emmean<- emmeans(N.on.tpi.lme3, ~tpi50_3cat, method = "tukey")
summary(N.on.TPI3.emmean, type ="response")
pairs(N.on.TPI3.emmean)
N.on.TPI3.emmean<-data.frame(N.on.TPI3.emmean)  

#### does tpi affect on or off.per.day (adj for time obs) ####
opd.tpi3.lm<- lme(on.per.day~ tpi50_3cat , random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(opd.tpi3.lm)
plot(opd.tpi3.lm)
pairs(emmeans(opd.tpi3.lm, ~tpi50_3cat))

opd.age.tpi3.lm<- lme(on.per.day~ tpi50_3cat + F.Age , random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(opd.age.tpi3.lm)
pairs(emmeans(opd.age.tpi3.lm, ~F.Age))

offpd.tpi3.lm<- lme(off.per.day~ tpi50_3cat , random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(offpd.tpi3.lm)
plot(offpd.tpi3.lm)
pairs(emmeans(offpd.tpi3.lm, ~tpi50_3cat))

#### Does aspect affect # on bouts? ####
N.on.aspect.lme<-glmer(number.on.bouts ~  f.aspect + (1|Nest.ID), family = poisson, data = bout_amb_dat,  na.action=na.pass)
summary(N.on.aspect.lme)

N.on.aspect.emmean<- emmeans(N.on.aspect.lme, ~f.aspect, method = "tukey")
summary(N.on.aspect.emmean, type ="response")
pairs(N.on.aspect.emmean)
 
# adjusted on or off bout per day
opd.aspect.lm<- lme(on.per.day~ f.aspect , random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(opd.aspect.lm)
plot(opd.aspect.lm)
pairs(emmeans(opd.aspect.lm, ~f.aspect))

offpd.aspect.lm<- lme(off.per.day~ f.aspect , random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(offpd.aspect.lm)
plot(offpd.aspect.lm)
pairs(emmeans(offpd.aspect.lm, ~f.aspect))
 
#### TPI 6 and On bouts ####
N.on.tpi6.lm<-glmer(number.on.bouts ~  tpi50_6cat + (1|Nest.ID), family = poisson,  data = bout_amb_dat, na.action=na.pass)
summary(N.on.tpi6.lm)
N.on.tpi6.emmean<- emmeans(N.on.tpi6.lm, ~tpi50_6cat, method = "tukey")
summary(N.on.tpi6.emmean, type ="response")
pairs(N.on.tpi6.emmean)

opd.tpi6.lm<- lme(on.per.day~ tpi50_6cat , random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(opd.tpi6.lm)
plot(opd.tpi6.lm)
pairs(emmeans(opd.tpi6.lm, ~tpi50_6cat))

# off per day
offpd.tpi6.lm<- lme(off.per.day~ tpi50_6cat , random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(offpd.tpi6.lm)
plot(offpd.tpi6.lm)
pairs(emmeans(offpd.tpi6.lm, ~tpi50_6cat))

#### on bout duration ####
hist(bout_amb_dat$T.on)
hist(log(bout_amb_dat$T.on)) # meets normal
bout_amb_dat$log.T.on<-log(bout_amb_dat$T.on)

# J.day & time on
T.on.day.lme<-lme(log.T.on~jday.s, random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.pass)
summary(T.on.day.lme)
plot(T.on.day.lme)
# Age
T.on.age.lme<-lme(log.T.on~F.Age, random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.pass)
summary(T.on.age.lme)
plot(T.on.age.lme)
T.on.age.emmean<- emmeans(T.on.age.lme, ~F.Age, method = "tukey")
summary(T.on.age.emmean, type ="response")
pairs(T.on.age.emmean)

#TPI3   
T.on.tpi3.lme<-lme(log.T.on~tpi50_3cat, random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.pass)
summary(T.on.tpi3.lme) # no diff if drop 321.12
#plot(T.off.tpi3.lme)
T.on.tpi3.emmean<- emmeans(T.on.tpi3.lme, ~tpi50_3cat, method = "tukey")
summary(T.on.tpi3.emmean, type ="response")
pairs(T.on.tpi3.emmean)

#tpi6  # ridge longer on than upper slope
T.on.tpi6.lme<-lme(log.T.on~tpi50_6cat, random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.pass)
summary(T.on.tpi6.lme) #no real diff if drop 321.12
plot(T.on.tpi6.lme)
T.on.tpi6.emmean<- emmeans(T.on.tpi6.lme, ~tpi50_6cat, method = "tukey")
summary(T.on.tpi6.emmean, type ="response")
pairs(T.on.tpi6.emmean)

#Aspect
T.on.aspect.lme<-lme(log.T.on~f.aspect, random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.pass)
summary(T.on.aspect.lme)
plot(T.on.aspect.lme)
T.on.aspect.emmean<- emmeans(T.on.aspect.lme, ~f.aspect, method = "tukey")
summary(T.on.aspect.emmean, type ="response")
pairs(T.on.aspect.emmean)

#### off bout duration ####
hist(bout_amb_dat$T.off)
hist(log(bout_amb_dat$T.off)) # better
summary(bout_amb_dat$T.off)
bout_amb_dat$log.T.off<-log(bout_amb_dat$T.off)
summary(bout_amb_dat$log.T.off)
bout_amb_dat.red$log.T.off<-log(bout_amb_dat.red$T.off)

# J.day 
T.off.day.lme<-lme(log.T.off~jday.s, random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(T.off.day.lme)
plot(T.off.day.lme)
# Age
T.off.age.lme<-lme(log.T.off~F.Age, random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(T.off.age.lme)
plot(T.off.age.lme)
T.off.age.emmean<- emmeans(T.off.age.lme, ~F.Age, method = "tukey")
summary(T.off.age.emmean, type ="response")
pairs(T.off.age.emmean)

#TPI3
T.off.tpi3.lme<-lme(log.T.off~tpi50_3cat, random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(T.off.tpi3.lme)
plot(T.off.tpi3.lme)
T.off.tpi3.emmean<- emmeans(T.off.tpi3.lme, ~tpi50_3cat, method = "tukey")
summary(T.off.tpi3.emmean, type ="response")
pairs(T.off.tpi3.emmean)

#tpi6  # upper slope shorter off than ridge (ridge longer off)
T.off.tpi6.lme<-lme(log.T.off~tpi50_6cat, random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(T.off.tpi6.lme) # ridge v u slop sig diff in both but post hoc p is lower using red. use full data
plot(T.off.tpi6.lme)
T.off.tpi6.emmean<- emmeans(T.off.tpi6.lme, ~tpi50_6cat, method = "tukey")
summary(T.off.tpi6.emmean, type ="response")
pairs(T.off.tpi6.emmean)

#Aspect
T.off.aspect.lme<-lme(log.T.off~f.aspect, random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(T.off.aspect.lme)
plot(T.off.aspect.lme)
T.off.aspect.emmean<- emmeans(T.off.aspect.lme, ~f.aspect, method = "tukey")
summary(T.off.aspect.emmean, type ="response")
pairs(T.off.aspect.emmean)

#### percent of time On the nest ####
hist(bout_amb_dat$perc.in)
range(bout_amb_dat$perc.in)
range(bout_amb_dat.red$perc.in)
summary(bout_amb_dat$perc.in)
which(bout_amb_dat$perc.in <= 50) # line 224 nest 321.12 first night of recording (jday 177. failed to predation jday 181. BHCO egg plus 3 others. 1 egg left in nest with a hole in it. overall, she has weird times compared to most. nest failed to predation.  made .red data set above and ran models with and without. see notes with analyses
hist(bout_amb_dat.red$perc.in)

# J.day 
T.perc.in.day.lme<-lme(perc.in~jday.s, random = ~1|Nest.ID, data = bout_amb_dat.red, na.action=na.omit)
summary(T.perc.in.day.lme)
plot(T.perc.in.day.lme)
# Age # keep 321.12 in for age
T.perc.in.age.lme<-lme(perc.in~F.Age, random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(T.perc.in.age.lme)
plot(T.perc.in.age.lme)
T.perc.in.age.emmean<- emmeans(T.perc.in.age.lme, ~F.Age, method = "tukey")
summary(T.perc.in.age.emmean, type ="response")
pairs(T.perc.in.age.emmean)

#TPI3 # no real diffs without 321 so keep
T.perc.in.tpi3.lme<-lme(perc.in~tpi50_3cat, random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(T.perc.in.tpi3.lme)
plot(T.perc.in.tpi3.lme)
T.perc.in.tpi3.emmean<- emmeans(T.perc.in.tpi3.lme, ~tpi50_3cat, method = "tukey")
summary(T.perc.in.tpi3.emmean, type ="response")
pairs(T.perc.in.tpi3.emmean)

#tpi6  #res better without 321 but no real diffs. keep
T.perc.in.tpi6.lme<-lme(perc.in~tpi50_6cat, random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(T.perc.in.tpi6.lme)
plot(T.perc.in.tpi6.lme)
T.perc.in.tpi6.emmean<- emmeans(T.perc.in.tpi6.lme, ~tpi50_6cat, method = "tukey")
summary(T.perc.in.tpi6.emmean, type ="response")
pairs(T.perc.in.tpi6.emmean)

#Aspect
T.perc.in.aspect.lme<-lme(perc.in~f.aspect, random = ~1|Nest.ID, data = bout_amb_dat, na.action=na.omit)
summary(T.perc.in.aspect.lme)
plot(T.perc.in.aspect.lme)
T.perc.in.aspect.emmean<- emmeans(T.perc.in.aspect.lme, ~f.aspect, method = "tukey")
summary(T.perc.in.aspect.emmean, type ="response")
pairs(T.perc.in.aspect.emmean)


#### does ambient temp affect incubation behavior ####
#### scale temp data ####
bout_amb_dat$T.min.s<-scale(bout_amb_dat$T.min)
bout_amb_dat$T.max.s<-scale(bout_amb_dat$T.max)
bout_amb_dat$T.mean.s<-scale(bout_amb_dat$T.mean)
bout_amb_dat$CV.s<- scale(bout_amb_dat$CV)
which(is.na(bout_amb_dat$T.min))

bout_amb_dat.red$T.min.s<-scale(bout_amb_dat.red$T.min)
bout_amb_dat.red$T.max.s<-scale(bout_amb_dat.red$T.max)
bout_amb_dat.red$T.mean.s<-scale(bout_amb_dat.red$T.mean)
bout_amb_dat.red$CV.s<- scale(bout_amb_dat.red$CV)


#### number of on bouts and temp NB: redid models with adjusted on.per.day ####
#T.min
Temp.on.lm<-lme(on.per.day~ T.min.s, random = ~ 1|Nest.ID, data = bout_amb_dat, na.action = na.omit) 
summary(Temp.on.lm)
plot(Temp.on.lm)
(5.017743^2)/(5.017743^2 + 6.750331^2) # ICC
describeBy(bout_amb_dat$on.per.day)

On.Tmin.pred.dat<-bout_amb_dat%>% filter(!is.na(T.min.s))
On.Tmin.pred.dat<-On.Tmin.pred.dat%>% filter(!is.na(on.per.day))
On.Tmin.pred.dat$pred<-predict(Temp.on.lm, type="response") #
On.Tmin.pred.dat$pred.pop<-predict(Temp.on.lm)
on.bouts.Tmin.fig<-ggplot(On.Tmin.pred.dat, aes(T.min, on.per.day))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+ #
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+
  geom_smooth(aes(y = pred.pop), method = lm, fill = NA, color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(0,65, 10), limits=c(0,65))+
  scale_x_continuous(breaks=seq(2,30, 4), limits=c(2,30))+
  labs(x = "Minimum ambient temperature (\u00B0C)", y = "On-bouts per day")+
  ggtitle("(A)")
on.bouts.Tmin.fig
ggsave(path = "Figs.OnBouts", "on.bouts.Tmin.fig.jpg", on.bouts.Tmin.fig, units = "in", height=4, width = 4, dpi = 600)

#T.max
Temp.max.on.lm<-lme(on.per.day~ T.max.s, random=~1|Nest.ID,  data = bout_amb_dat, na.action = na.omit) 
summary(Temp.max.on.lm)
plot(Temp.max.on.lm)
5.183565^2/(5.183565^2 + 7.053986^2) # icc
describeBy(bout_amb_dat$off.per.day)

On.Tmax.pred.dat<-bout_amb_dat%>% filter(!is.na(T.max.s))
On.Tmax.pred.dat<-On.Tmax.pred.dat%>% filter(!is.na(off.per.day))
On.Tmax.pred.dat$pred<-predict(Temp.max.on.lm, type="response") #
On.Tmax.pred.dat$pred.pop<-predict(Temp.max.on.lm)
on.bouts.Tmax.fig<-ggplot(On.Tmax.pred.dat, aes(T.max, on.per.day))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+ #
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+
  geom_smooth(aes(y = pred.pop), method = lm, fill = NA, color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(0,65, 10), limits=c(0,65))+
  scale_x_continuous(breaks=seq(12,42, 4), limits=c(12,42))+
  labs(x = "Maximum ambient temperature (\u00B0C)", y = "On-bouts per day")+
  ggtitle("(B)")
on.bouts.Tmax.fig
ggsave(path = "Figs.OnBouts", "on.bouts.Tmax.fig.jpg", on.bouts.Tmax.fig, units = "in", height=4, width = 4, dpi = 600)

#T.mean
Temp.mean.on.lm<-lme(on.per.day~ T.mean.s, random = ~1|Nest.ID,  data = bout_amb_dat, na.action = na.omit) 
summary(Temp.mean.on.lm)

On.Tmean.pred.dat<-bout_amb_dat%>% filter(!is.na(T.mean.s))
On.Tmean.pred.dat<-On.Tmean.pred.dat%>%filter(!is.na(on.per.day))
On.Tmean.pred.dat$pred<-predict(Temp.mean.on.lm) #
On.Tmean.pred.dat$pred.pop<-predict(Temp.mean.on.lm,  re.form = NA)
on.bouts.Tmean.fig<-ggplot(On.Tmean.pred.dat, aes(T.mean, on.per.day))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+
  geom_smooth(aes(y = pred.pop), method = "lm", fill = NA, color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(0,65, 10), limits=c(0,65))+
  #scale_x_continuous(breaks=seq(12,42, 4), limits=c(12,42))+
  labs(x = "Mean ambient temperature (\u00B0C)", y = "On-bouts per day")+
  ggtitle("(C)")
on.bouts.Tmean.fig
ggsave(path = "Figs.OnBouts", "on.bouts.Tmean.fig.jpg", on.bouts.Tmean.fig, units = "in", height=4, width = 4, dpi = 600)

# CV
Temp.CV.on.lm<-lme(on.per.day~ CV.s, random = ~1|Nest.ID, data = bout_amb_dat, na.action = na.omit) 
summary(Temp.CV.on.lm)
plot(Temp.CV.on.lm)

On.CV.pred.dat<-bout_amb_dat%>% filter(!is.na(CV.s))
On.CV.pred.dat<-On.CV.pred.dat%>%filter(!is.na(on.per.day))
On.CV.pred.dat$pred<-predict(Temp.CV.on.lm) #
On.CV.pred.dat$pred.pop<-predict(Temp.CV.on.lm)
on.bouts.CV.fig<-ggplot(On.CV.pred.dat, aes(CV, on.per.day))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+
  geom_smooth(aes(y = pred.pop), method = "lm", fill = NA, color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(0,65, 10), limits=c(0,65))+
  #scale_x_continuous(breaks=seq(12,42, 4), limits=c(12,42))+
  labs(x = "Variation (CV) in ambient temperature", y = "On-bouts per day")+
  ggtitle("(D)")
on.bouts.CV.fig
ggsave(path = "Figs.OnBouts", "on.bouts.CV.fig.jpg", on.bouts.CV.fig, units = "in", height=4, width = 4, dpi = 600)

#### number of off bouts and temp: Re did using off per day to adust for partial days ####
# note taking more on and off bouts when it is cold - can't take an on without an off. prob show just one or the other
#T.min
Temp.off.lm<-lme(off.per.day~ T.min.s, random = ~1|Nest.ID, data = bout_amb_dat, na.action = na.omit) 
summary(Temp.off.lm)
4.997687^2/( 4.997687^2 + 6.784627^2 ) # icc
describeBy(bout_amb_dat$off.per.day)
plot(Temp.off.lm)

off.Tmin.pred.dat<-bout_amb_dat%>% filter(!is.na(T.min.s))
off.Tmin.pred.dat<-off.Tmin.pred.dat%>%filter(!is.na(off.per.day))
off.Tmin.pred.dat$pred<-predict(Temp.off.lm) #
off.Tmin.pred.dat$pred.pop<-predict(Temp.off.lm)
off.bouts.Tmin.fig<-ggplot(off.Tmin.pred.dat, aes(T.min, off.per.day))+
  theme_classic()+
  geom_jitter(color= "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+
  geom_smooth(aes(y = pred.pop),method = "lm", fill = NA, color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(0,65, 10), limits=c(0,65))+
  scale_x_continuous(breaks=seq(2,30, 4), limits=c(2,30))+
  labs(x = "Minimum ambient temperature (\u00B0C)", y = "Off-bouts per day")+
  ggtitle("(A)")
off.bouts.Tmin.fig
ggsave(path = "Figs.OnBouts", "off.bouts.Tmin.fig.jpg", off.bouts.Tmin.fig, units = "in", height=4, width = 4, dpi = 600)

#T.max
Temp.max.off.lm<-lme(off.per.day~ T.max.s, random = ~1|Nest.ID, data = bout_amb_dat, na.action = na.omit) 
summary(Temp.max.off.lm)
plot(Temp.max.off.lm)

off.Tmax.pred.dat<-bout_amb_dat%>% filter(!is.na(T.max.s))
off.Tmax.pred.dat<-off.Tmax.pred.dat%>%filter(!is.na(off.per.day))
off.Tmax.pred.dat$pred<-predict(Temp.max.off.lm) #
off.Tmax.pred.dat$pred.pop<-predict(Temp.max.off.lm)
off.bouts.Tmax.fig<-ggplot(off.Tmax.pred.dat, aes(T.max, off.per.day))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+
  geom_smooth(aes(y = pred.pop),method = "lm", fill = NA, color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(0,65, 10), limits=c(0,65))+
  scale_x_continuous(breaks=seq(12,42, 4), limits=c(12,42))+
  labs(x = "Maximum ambient temperature (\u00B0C)", y = "Off-bouts per day")+
  ggtitle("(B)")
off.bouts.Tmax.fig
ggsave(path = "Figs.OnBouts", "off.bouts.Tmax.fig.jpg", off.bouts.Tmin.fig, units = "in", height=4, width = 4, dpi = 600)

#T.mean
Temp.mean.off.lm<-lme(off.per.day~ T.mean.s, random = ~1|Nest.ID, data = bout_amb_dat, na.action = na.omit) 
summary(Temp.mean.off.lm)
plot(Temp.mean.off.lm)

off.Tmean.pred.dat<-bout_amb_dat%>% filter(!is.na(T.mean.s))
off.Tmean.pred.dat<-off.Tmean.pred.dat%>%filter(!is.na(off.per.day))
off.Tmean.pred.dat$pred<-predict(Temp.mean.off.lm) #
off.Tmean.pred.dat$pred.pop<-predict(Temp.mean.off.lm)
off.bouts.Tmean.fig<-ggplot(off.Tmean.pred.dat, aes(T.mean, off.per.day))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+
  geom_smooth(aes(y = pred.pop), method = "lm", fill = NA,color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(0,65, 10), limits=c(0,65))+
  #scale_x_continuous(breaks=seq(12,42, 4), limits=c(12,42))+
  labs(x = "Mean ambient temperature (\u00B0C)", y = "Off-bouts per day")+
  ggtitle("(C)")
off.bouts.Tmean.fig
ggsave(path = "Figs.OnBouts", "off.bouts.Tmean.fig.jpg", off.bouts.Tmean.fig, units = "in", height=4, width = 4, dpi = 600)

# CV
Temp.CV.off.lm<-lme(off.per.day~ CV.s,random = ~1|Nest.ID,data = bout_amb_dat, na.action = na.omit) 
summary(Temp.CV.off.lm)
plot(Temp.CV.off.lm)

off.CV.pred.dat<-bout_amb_dat%>% filter(!is.na(CV.s))
off.CV.pred.dat<-off.CV.pred.dat%>% filter(!is.na(off.per.day))
off.CV.pred.dat$pred<-predict(Temp.CV.off.lm) #
off.CV.pred.dat$pred.pop<-predict(Temp.CV.off.lm) 
off.bouts.CV.fig<-ggplot(off.CV.pred.dat, aes(CV, off.per.day))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+
  geom_smooth(aes(y = pred.pop), method = "lm", fill = NA, color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(0,65, 10), limits=c(0,65))+
  #scale_x_continuous(breaks=seq(12,42, 4), limits=c(12,42))+
  labs(x = "Variation (CV) in ambient temperature", y = "Off-bouts per day")+
  ggtitle("(D)")
off.bouts.CV.fig
ggsave(path = "Figs.OnBouts", "off.bouts.CV.fig.jpg", off.bouts.CV.fig, units = "in", height=4, width = 4, dpi = 600)

#### Mean time on on-bouts (T.on) by ambient temps ####
# NB skip to below to analyses of reduced data
summary(bout_amb_dat$T.on)
describeBy(bout_amb_dat$T.on)
which(bout_amb_dat$T.on > 400)
#T.min
Temp.Ton.lm<-lme(log.T.on~ T.min.s, random = ~1|Nest.ID, data = bout_amb_dat, na.action = na.omit) 
summary(Temp.Ton.lm)
plot(Temp.Ton.lm)

Ton.Tmin.pred.dat<-bout_amb_dat%>% filter(!is.na(T.min.s))
Ton.Tmin.pred.dat$pred<-predict(Temp.Ton.lm, type="response") #
Ton.Tmin.pred.dat$exp.pred<- exp(Ton.Tmin.pred.dat$pred) #exponent to get back on scale
Ton.bouts.Tmin.fig<-ggplot(Ton.Tmin.pred.dat, aes(T.min,T.on))+
  theme_classic()+
  geom_point()+
  geom_line(aes(y = exp(predict(Temp.Ton.lm, level = 0))), color = "black",size = 1.1)+ # pop line
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray60", size = 0.4, linetype = "dashed" )+
  scale_y_continuous(breaks=seq(0,200, 20), limits=c(0,200))+
  scale_x_continuous(breaks=seq(5,30, 5), limits=c(5,30))+
  labs(x = "Minimum ambient temperature (\u00B0C)", y = "On-bout duration (min)")+
  ggtitle("(A)")
Ton.bouts.Tmin.fig
ggsave(path = "Figs.OnBouts", "Ton.bouts.Tmin.fig.jpg", Ton.bouts.Tmin.fig, units = "in", height=4, width = 4, dpi = 600)

#T.max
Temp.max.Ton.lm<-lme(log.T.on~ T.max.s, random = ~1|Nest.ID,  data = bout_amb_dat, na.action = na.omit) 
summary(Temp.max.Ton.lm)

Ton.Tmax.pred.dat<-bout_amb_dat%>% filter(!is.na(T.max.s))
Ton.Tmax.pred.dat$pred<-predict(Temp.max.Ton.lm, type="response") #
Ton.Tmax.pred.dat$exp.pred<- exp(Ton.Tmax.pred.dat$pred) #exponent to get back on scale
Ton.bouts.Tmax.fig<-ggplot(Ton.Tmax.pred.dat, aes(T.max, T.on))+
  theme_classic()+
  geom_point()+
  geom_line(aes(y = exp(predict(Temp.max.Ton.lm, level = 0))), color = "black",size = 1.1)+ # pop line
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray60", size = 0.4, linetype = "dashed")+ # add by Nest
   scale_y_continuous(breaks=seq(0,200, 20), limits=c(0,200))+
  scale_x_continuous(breaks=seq(12,42, 4), limits=c(12,42))+
  labs(x = "Maximum ambient temperature (\u00B0C)", y = "On-bout duration (min)")+
  ggtitle("B")
Ton.bouts.Tmax.fig
ggsave(path = "Figs.OnBouts", "Ton.bouts.Tmax.fig.jpg", Ton.bouts.Tmax.fig, units = "in", height=4, width = 4, dpi = 600)

#T.mean
Temp.mean.Ton.lm<-lme(log.T.on~ T.mean.s, random = ~1|Nest.ID,  data = bout_amb_dat, na.action = na.omit) 
summary(Temp.mean.Ton.lm)
plot(Temp.mean.Ton.lm)

Ton.Tmean.pred.dat<-bout_amb_dat%>% filter(!is.na(T.mean.s))
Ton.Tmean.pred.dat$pred<-predict(Temp.mean.Ton.lm, type="response") #
Ton.Tmean.pred.dat$exp.pred<- exp(Ton.Tmean.pred.dat$pred) #exponent to get back on scale
Ton.bouts.Tmean.fig<-ggplot(Ton.Tmean.pred.dat, aes(T.mean, T.on))+
  theme_classic()+
  geom_point()+
  geom_line(aes(y = exp(predict(Temp.mean.Ton.lm, level = 0))), color = "black",size = 1.1)+ # pop line
  geom_line(aes(y= exp.pred, group = Nest.ID, method = "lm"), color = "gray60", size = 0.4, linetype = "dashed")+ # add by Nest
  scale_y_continuous(breaks=seq(0,200, 20), limits=c(0,200))+
  scale_x_continuous(breaks=seq(10,30, 5), limits=c(10,30))+
  labs(x = "Mean ambient temperature (\u00B0C)", y = "On-bout duration (min)")+
  ggtitle("C")
Ton.bouts.Tmean.fig
ggsave(path = "Figs.OnBouts", "Ton.bouts.Tmean.fig.jpg", Ton.bouts.Tmean.fig, units = "in", height=4, width = 4, dpi = 600)

# CV
Temp.CV.Ton.lm<-lme(log.T.on~ CV.s, random = ~1|Nest.ID,  data = bout_amb_dat, na.action = na.omit) 
summary(Temp.CV.Ton.lm)

Ton.CV.pred.dat<-bout_amb_dat%>% filter(!is.na(CV.s))
Ton.CV.pred.dat$pred<-predict(Temp.CV.Ton.lm, type="response") #
Ton.CV.pred.dat$exp.pred<- exp(Ton.CV.pred.dat$pred) #exponent to get back on scale
Ton.bouts.TCV.fig<-ggplot(Ton.CV.pred.dat, aes(CV, T.on))+
  theme_classic()+
  geom_point()+
  geom_line(aes(y = exp(predict(Temp.CV.Ton.lm, level = 0))), color = "black",size = 1.1)+ # pop line
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray60", size = 0.4, linetype = "dashed" )+ # add by Nest
  scale_y_continuous(breaks=seq(0,200, 20), limits=c(0,200))+
  #scale_x_continuous(breaks=seq(12,42, 4), limits=c(12,42))+
  labs(x = "Variation (CV) in ambient temperature", y = "On-bout duration (min)")+
  ggtitle("D")
Ton.bouts.TCV.fig
ggsave(path = "Figs.OnBouts", "Ton.bouts.TCV.fig.jpg", Ton.bouts.TCV.fig, units = "in", height=4, width = 4, dpi = 600)


#### Mean time off-bouts (T.off) by ambient temps ####
summary(bout_amb_dat$T.off)
describeBy(bout_amb_dat$T.off)
which(bout_amb_dat$T.off > 100)
#T.min
Temp.Toff.lm<-lme(log.T.off~ T.min.s, random = ~1|Nest.ID, data = bout_amb_dat, na.action = na.omit) 
summary(Temp.Toff.lm)
plot(Temp.Toff.lm)

Toff.Tmin.pred.dat<-bout_amb_dat%>% filter(!is.na(T.min.s))
Toff.Tmin.pred.dat<-Toff.Tmin.pred.dat%>%filter(!is.na(log.T.off))
Toff.Tmin.pred.dat$pred<-predict(Temp.Toff.lm, type="response") #
Toff.Tmin.pred.dat$exp.pred<- exp(Toff.Tmin.pred.dat$pred) #exponent to get back on scale
Toff.bouts.Tmin.fig<-ggplot(Toff.Tmin.pred.dat, aes(T.min,T.off))+
  theme_classic()+
  geom_point()+
  geom_line(aes(y = exp(predict(Temp.Toff.lm, level = 0))), color = "black",size = 1.1)+ # pop line
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray60", size = 0.4, linetype = "dashed" )+
  scale_y_continuous(breaks=seq(0,25, 5), limits=c(0,25))+
  #scale_x_continuous(breaks=seq(5,30, 5), limits=c(5,30))+
  labs(x = "Minimum ambient temperature (\u00B0C)", y = "Off-bout duration (min)")+
  ggtitle("A")
Toff.bouts.Tmin.fig
ggsave(path = "Figs.OnBouts", "Toff.bouts.Tmin.fig.jpg", Toff.bouts.Tmin.fig, units = "in", height=4, width = 4, dpi = 600)

#T.max
Temp.max.Toff.lm<-lme(log.T.off~ T.max.s, random = ~1|Nest.ID,  data = bout_amb_dat, na.action = na.omit) 
summary(Temp.max.Toff.lm)

Toff.Tmax.pred.dat<-bout_amb_dat%>% filter(!is.na(T.max.s))
Toff.Tmax.pred.dat<-Toff.Tmax.pred.dat%>%filter(!is.na(log.T.off))
Toff.Tmax.pred.dat$pred<-predict(Temp.max.Toff.lm, type="response") #
Toff.Tmax.pred.dat$exp.pred<- exp(Toff.Tmax.pred.dat$pred) #exponent to get back on scale
Toff.bouts.Tmax.fig<-ggplot(Toff.Tmax.pred.dat, aes(T.max, T.off))+
  theme_classic()+
  geom_point()+
  geom_line(aes(y = exp(predict(Temp.max.Toff.lm, level = 0))), color = "black",size = 1.1)+ # pop line
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray60", size = 0.4, linetype = "dashed")+ # add by Nest
  scale_y_continuous(breaks=seq(0,25, 5), limits=c(0,25))+
  #scale_x_continuous(breaks=seq(12,42, 4), limits=c(12,42))+
  labs(x = "Maximum ambient temperature (\u00B0C)", y = "Off-bout duration (min)")+
  ggtitle("B")
Toff.bouts.Tmax.fig
ggsave(path = "Figs.OnBouts", "Toff.bouts.Tmax.fig.jpg", Toff.bouts.Tmax.fig, units = "in", height=4, width = 4, dpi = 600)

#T.mean
Temp.mean.Toff.lm<-lme(log.T.off~ T.mean.s, random = ~1|Nest.ID,  data = bout_amb_dat, na.action = na.omit) 
summary(Temp.mean.Toff.lm)
plot(Temp.mean.Toff.lm)

Toff.Tmean.pred.dat<-bout_amb_dat%>% filter(!is.na(T.mean.s))
Toff.Tmean.pred.dat<-Toff.Tmean.pred.dat%>%filter(!is.na(log.T.off))
Toff.Tmean.pred.dat$pred<-predict(Temp.mean.Toff.lm, type="response") #
Toff.Tmean.pred.dat$exp.pred<- exp(Toff.Tmean.pred.dat$pred) #exponent to get back on scale
Toff.bouts.Tmean.fig<-ggplot(Toff.Tmean.pred.dat, aes(T.mean, T.off))+
  theme_classic()+
  geom_point()+
  geom_line(aes(y = exp(predict(Temp.mean.Toff.lm, level = 0))), color = "black",size = 1.1)+ # pop line
  geom_line(aes(y= exp.pred, group = Nest.ID, method = "lm"), color = "gray60", size = 0.4, linetype = "dashed")+ # add by Nest
  scale_y_continuous(breaks=seq(0,25, 5), limits=c(0,25))+
  #scale_x_continuous(breaks=seq(10,30, 5), limits=c(10,30))+
  labs(x = "Mean ambient temperature (\u00B0C)", y = "Off-bout duration (min)")+
  ggtitle("C")
Toff.bouts.Tmean.fig
ggsave(path = "Figs.OnBouts", "Toff.bouts.Tmean.fig.jpg", Toff.bouts.Tmean.fig, units = "in", height=4, width = 4, dpi = 600)

# CV
Temp.CV.Toff.lm<-lme(log.T.off~ CV.s, random = ~1|Nest.ID,  data = bout_amb_dat, na.action = na.omit) 
summary(Temp.CV.Toff.lm)

Toff.CV.pred.dat<-bout_amb_dat%>% filter(!is.na(CV.s))
Toff.CV.pred.dat<-Toff.CV.pred.dat%>%filter(!is.na(log.T.off))
Toff.CV.pred.dat$pred<-predict(Temp.CV.Toff.lm, type="response") #
Toff.CV.pred.dat$exp.pred<- exp(Toff.CV.pred.dat$pred) #exponent to get back on scale
Toff.bouts.TCV.fig<-ggplot(Toff.CV.pred.dat, aes(CV, T.off))+
  theme_classic()+
  geom_point()+
  geom_line(aes(y = exp(predict(Temp.CV.Toff.lm, level = 0))), color = "black",size = 1.1)+ # pop line
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray60", size = 0.4, linetype = "dashed" )+ # add by Nest
  scale_y_continuous(breaks=seq(0,25, 5), limits=c(0,25))+
  #scale_x_continuous(breaks=seq(12,42, 4), limits=c(12,42))+
  labs(x = "Variation (CV) in ambient temperature", y = "Off-bout duration (min)")+
  ggtitle("D")
Toff.bouts.TCV.fig
ggsave(path = "Figs.OnBouts", "Toff.bouts.TCV.fig.jpg", Toff.bouts.TCV.fig, units = "in", height=4, width = 4, dpi = 600)

#### duration onbout ~ duration offbout ####
duration.on.off.lm<-lme(log.T.on~ log.T.off, random = ~1|Nest.ID,  data = bout_amb_dat, na.action = na.omit) 
summary(duration.on.off.lm)

duration.on.off.pred.dat<-bout_amb_dat%>%filter(!is.na(log.T.off))
duration.on.off.pred.dat$pred<-predict(duration.on.off.lm, type="response") #
duration.on.off.pred.dat$exp.pred<- exp(duration.on.off.pred.dat$pred) #exponent to get back on scale
duration.on.off.fig<-ggplot(duration.on.off.pred.dat, aes(T.off, T.on))+
  theme_classic()+
  geom_point()+
  geom_line(aes(y = exp(predict(duration.on.off.lm, level = 0))), color = "black",size = 1.1)+ # pop line
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray60", size = 0.4, linetype = "dashed" )+ # add by Nest
  scale_y_continuous(breaks=seq(0,200, 20), limits=c(0,200))+
  scale_x_continuous(breaks=seq(0,70, 20), limits=c(0,70))+
  labs(x = "Off-bout duration ", y = "On-bout duration (min)")
duration.on.off.fig
ggsave(path = "Figs.OnBouts", "duration.on.off.fig.jpg", duration.on.off.fig, units = "in", height=4, width = 4, dpi = 600)

#### try duration without outliers ####
dur.on.red.dat<-filter(bout_amb_dat, T.on < 120 & T.off < 30) 
summary(bout_amb_dat$T.on)
length(which(bout_amb_dat$T.on > 120))
#check which nests
rem.nests<-bout_amb_dat[c(80, 102, 120:121,129,134:135, 144,153,169, 223, 227, 306,347,361,362),1]
rem.nests
length(which(bout_amb_dat$T.on > 90))
length(which(bout_amb_dat$T.off >30))
rem.nests.off<-bout_amb_dat[c(224:228, 355, 356),1]
rem.nests.off
#desc stats
summary(bout_amb_dat$T.off)
describeBy(bout_amb_dat$T.off)
summary(dur.on.red.dat$T.on)
summary(dur.on.red.dat$T.off)
summary(dur.on.red.dat$number.on.bouts)
summary(bout_amb_dat$number.on.bouts)
summary(bout_amb_dat$number.off.bouts)
summary(dur.on.red.dat$number.off.bouts)
summary(bout_amb_dat$perc.in)
summary(dur.on.red.dat$perc.in)

#### on ~ off reduced data ####
duration.on.off.lm.red<-lme(log.T.on~ log.T.off, random = ~1|Nest.ID,  data = dur.on.red.dat, na.action = na.omit) 
summary(duration.on.off.lm.red)
plot(duration.on.off.lm.red)
0.3480861^2/(0.3480861^2 + 0.2565419^2)

duration.on.off.red.pred.dat<-dur.on.red.dat%>%filter(!is.na(log.T.off))
duration.on.off.red.pred.dat$pred<-predict(duration.on.off.lm.red, type="response") #
duration.on.off.red.pred.dat$exp.pred<- exp(duration.on.off.red.pred.dat$pred) #exponent to get back on scale
duration.on.off.red.fig<-ggplot(duration.on.off.red.pred.dat, aes(T.off, T.on))+
  theme_classic()+
  geom_jitter(color= "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest
  geom_line(aes(y = exp(predict(duration.on.off.lm.red, level = 0))), color = "black",size = 1.1)+ # pop line
 # scale_y_continuous(breaks=seq(0,200, 20), limits=c(0,200))+
  #scale_x_continuous(breaks=seq(0,70, 20), limits=c(0,70))+
  labs(x = "Off-bout duration (min)", y = "On-bout duration (min)")
duration.on.off.red.fig
ggsave(path = "Figs.OnBouts", "duration.on.off.red.fig.jpg", duration.on.off.red.fig, units = "in", height=4, width = 4, dpi = 600)

#### Mean time on on-bouts (T.on) with reduced data set by ambient temps ####

#T.min
Temp.Ton.lm.red<-lme(log.T.on~ T.min.s, random = ~1|Nest.ID, data = dur.on.red.dat, na.action = na.omit) 
summary(Temp.Ton.lm.red)
plot(Temp.Ton.lm.red)
0.471431^2/(0.471431^2 + 0.2559237^2)
Ton.Tmin.pred.dat.red<-dur.on.red.dat%>% filter(!is.na(T.min.s))
Ton.Tmin.pred.dat.red$pred<-predict(Temp.Ton.lm.red, type="response") #
Ton.Tmin.pred.dat.red$exp.pred<- exp(Ton.Tmin.pred.dat.red$pred) #exponent to get back on scale
Ton.bouts.Tmin.red.fig<-ggplot(Ton.Tmin.pred.dat.red, aes(T.min,T.on))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+
  geom_line(aes(y = exp(predict(Temp.Ton.lm.red, level = 0))), color = "black",size = 1.1)+ # pop line
  #scale_y_continuous(breaks=seq(0,200, 20), limits=c(0,200))+
  #scale_x_continuous(breaks=seq(5,30, 5), limits=c(5,30))+
  labs(x = "Minimum ambient temperature (\u00B0C)", y = "On-bout duration (min)")+
  ggtitle("(A)")
Ton.bouts.Tmin.red.fig
ggsave(path = "Figs.OnBouts", "Ton.bouts.Tmin.red.fig.jpg", Ton.bouts.Tmin.red.fig, units = "in", height=4, width = 4, dpi = 600)

#T.max
Temp.max.Ton.lm.red<-lme(log.T.on~ T.max.s, random = ~1|Nest.ID,  data = dur.on.red.dat, na.action = na.omit) 
summary(Temp.max.Ton.lm.red)
0.495089^2/(0.495089^2 + 0.2536246^2)
Ton.Tmax.pred.dat.red<-dur.on.red.dat%>% filter(!is.na(T.max.s))
Ton.Tmax.pred.dat.red$pred<-predict(Temp.max.Ton.lm.red, type="response") #
Ton.Tmax.pred.dat.red$exp.pred<- exp(Ton.Tmax.pred.dat.red$pred) #exponent to get back on scale
Ton.bouts.Tmax.red.fig<-ggplot(Ton.Tmax.pred.dat.red, aes(T.max, T.on))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed")+ # add by Nest
  geom_line(aes(y = exp(predict(Temp.max.Ton.lm.red, level = 0))), color = "black",size = 1.1)+ # pop line
  #scale_y_continuous(breaks=seq(0,200, 20), limits=c(0,200))+
  #scale_x_continuous(breaks=seq(12,42, 4), limits=c(12,42))+
  labs(x = "Maximum ambient temperature (\u00B0C)", y = "On-bout duration (min)")+
  ggtitle("(B)")
Ton.bouts.Tmax.red.fig
ggsave(path = "Figs.OnBouts", "Ton.bouts.Tmax.red.fig.jpg", Ton.bouts.Tmax.red.fig, units = "in", height=4, width = 4, dpi = 600)

#T.mean
Temp.mean.Ton.lm.red<-lme(log.T.on~ T.mean.s, random = ~1|Nest.ID,  data = dur.on.red.dat, na.action = na.omit) 
summary(Temp.mean.Ton.lm.red)
plot(Temp.mean.Ton.lm.red)
0.4962858^2/(0.4962858^2 + 0.249844^2)
Ton.Tmean.pred.dat.red<-dur.on.red.dat%>% filter(!is.na(T.mean.s))
Ton.Tmean.pred.dat.red$pred<-predict(Temp.mean.Ton.lm.red, type="response") #
Ton.Tmean.pred.dat.red$exp.pred<- exp(Ton.Tmean.pred.dat.red$pred) #exponent to get back on scale
Ton.bouts.Tmean.red.fig<-ggplot(Ton.Tmean.pred.dat.red, aes(T.mean, T.on))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed")+ # add by Nest
  geom_line(aes(y = exp(predict(Temp.mean.Ton.lm.red, level = 0))), color = "black",size = 1.1)+ # pop line
  #scale_y_continuous(breaks=seq(0,200, 20), limits=c(0,200))+
  #scale_x_continuous(breaks=seq(10,30, 5), limits=c(10,30))+
  labs(x = "Mean ambient temperature (\u00B0C)", y = "On-bout duration (min)")+
  ggtitle("(C)")
Ton.bouts.Tmean.red.fig
ggsave(path = "Figs.OnBouts", "Ton.bouts.Tmean.red.fig.jpg", Ton.bouts.Tmean.red.fig, units = "in", height=4, width = 4, dpi = 600)

# CV
Temp.CV.Ton.lm.red<-lme(log.T.on~ CV.s, random = ~1|Nest.ID,  data = dur.on.red.dat, na.action = na.omit)

summary(Temp.CV.Ton.lm.red)
0.4528945^2/(0.4528945^2 + 0.2678821^2)
Ton.CV.pred.dat.red<-dur.on.red.dat%>% filter(!is.na(CV.s))
Ton.CV.pred.dat.red$pred<-predict(Temp.CV.Ton.lm.red, type="response") #
Ton.CV.pred.dat.red$exp.pred<- exp(Ton.CV.pred.dat.red$pred) #exponent to get back on scale
Ton.bouts.TCV.red.fig<-ggplot(Ton.CV.pred.dat.red, aes(CV, T.on))+
  theme_classic()+
  geom_jitter(color="gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray60", size = 0.4, linetype = "dashed" )+ # add by Nest
  geom_line(aes(y = exp(predict(Temp.CV.Ton.lm.red, level = 0))), color = "black",size = 1.1)+ # pop line
  #scale_y_continuous(breaks=seq(0,200, 20), limits=c(0,200))+
  #scale_x_continuous(breaks=seq(12,42, 4), limits=c(12,42))+
  labs(x = "Variation (CV) in ambient temperature", y = "On-bout duration (min)")+
  ggtitle("(D)")
Ton.bouts.TCV.red.fig
ggsave(path = "Figs.OnBouts", "Ton.bouts.TCV.red.fig.jpg", Ton.bouts.TCV.red.fig, units = "in", height=4, width = 4, dpi = 600)

#### Mean time off-bouts (T.off) by ambient temps with reduced data ####
summary(bout_amb_dat$T.off)
describeBy(bout_amb_dat$T.off)
which(bout_amb_dat$T.off > 100)
#T.min
Temp.Toff.lm.red<-lme(log.T.off~ T.min.s, random = ~1|Nest.ID, data = dur.on.red.dat, na.action = na.omit) 
summary(Temp.Toff.lm.red)
plot(Temp.Toff.lm.red)
0.2529385^2/(0.2529385^2 +  0.2199547^2)
Toff.Tmin.pred.dat.red<-dur.on.red.dat%>% filter(!is.na(T.min.s))
Toff.Tmin.pred.dat.red<-Toff.Tmin.pred.dat.red%>%filter(!is.na(log.T.off))
Toff.Tmin.pred.dat.red$pred<-predict(Temp.Toff.lm.red, type="response") #
Toff.Tmin.pred.dat.red$exp.pred<- exp(Toff.Tmin.pred.dat.red$pred) #exponent to get back on scale
Toff.bouts.Tmin.red.fig<-ggplot(Toff.Tmin.pred.dat.red, aes(T.min,T.off))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+
  geom_line(aes(y = exp(predict(Temp.Toff.lm.red, level = 0))), color = "black",size = 1.1)+ # pop line
  #scale_y_continuous(breaks=seq(0,25, 5), limits=c(0,25))+
  #scale_x_continuous(breaks=seq(5,30, 5), limits=c(5,30))+
  labs(x = "Minimum ambient temperature (\u00B0C)", y = "Off-bout duration (min)")+
  ggtitle("(A)")
Toff.bouts.Tmin.red.fig
#ggsave(path = "Figs.OnBouts", "Toff.bouts.Tmin.fig.jpg", Toff.bouts.Tmin.fig, units = "in", height=4, width = 4, dpi = 600)

#T.max
Temp.max.Toff.lm.red<-lme(log.T.off~ T.max.s, random = ~1|Nest.ID,  data = dur.on.red.dat, na.action = na.omit) 
summary(Temp.max.Toff.lm.red)
0.2792172^2/(0.2792172 ^2 + 0.2129673^2) 
Toff.Tmax.pred.dat.red<-dur.on.red.dat%>% filter(!is.na(T.max.s))
Toff.Tmax.pred.dat.red<-Toff.Tmax.pred.dat.red%>%filter(!is.na(log.T.off))
Toff.Tmax.pred.dat.red$pred<-predict(Temp.max.Toff.lm.red, type="response") #
Toff.Tmax.pred.dat.red$exp.pred<- exp(Toff.Tmax.pred.dat.red$pred) #exponent to get back on scale
Toff.bouts.Tmax.red.fig<-ggplot(Toff.Tmax.pred.dat.red, aes(T.max, T.off))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed")+ # add by Nest
  geom_line(aes(y = exp(predict(Temp.max.Toff.lm.red, level = 0))), color = "black",size = 1.1)+ # pop line
  #scale_y_continuous(breaks=seq(0,25, 5), limits=c(0,25))+
  #scale_x_continuous(breaks=seq(12,42, 4), limits=c(12,42))+
  labs(x = "Maximum ambient temperature (\u00B0C)", y = "Off-bout duration (min)")+
  ggtitle("(B)")
Toff.bouts.Tmax.red.fig
ggsave(path = "Figs.OnBouts", "Toff.bouts.Tmax.red.fig.jpg", Toff.bouts.Tmax.red.fig, units = "in", height=4, width = 4, dpi = 600)

#T.mean
Temp.mean.Toff.lm.red<-lme(log.T.off~ T.mean.s, random = ~1|Nest.ID,  data = dur.on.red.dat, na.action = na.omit) 
summary(Temp.mean.Toff.lm.red)
plot(Temp.mean.Toff.lm.red)
0.2693477^2/(0.2693477^2 + 0.2152657^2)
Toff.Tmean.pred.dat.red<-dur.on.red.dat%>% filter(!is.na(T.mean.s))
Toff.Tmean.pred.dat.red<-Toff.Tmean.pred.dat.red%>%filter(!is.na(log.T.off))
Toff.Tmean.pred.dat.red$pred<-predict(Temp.mean.Toff.lm.red, type="response") #
Toff.Tmean.pred.dat.red$exp.pred<- exp(Toff.Tmean.pred.dat.red$pred) #exponent to get back on scale
Toff.bouts.Tmean.red.fig<-ggplot(Toff.Tmean.pred.dat.red, aes(T.mean, T.off))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed")+ # add by Nest
  geom_line(aes(y = exp(predict(Temp.mean.Toff.lm.red, level = 0))), color = "black",size = 1.1)+ # pop line
  #scale_y_continuous(breaks=seq(0,25, 5), limits=c(0,25))+
  #scale_x_continuous(breaks=seq(10,30, 5), limits=c(10,30))+
  labs(x = "Mean ambient temperature (\u00B0C)", y = "Off-bout duration (min)")+
  ggtitle("(C)")
Toff.bouts.Tmean.red.fig
ggsave(path = "Figs.OnBouts", "Toff.bouts.Tmean.fig.jpg", Toff.bouts.Tmean.fig, units = "in", height=4, width = 4, dpi = 600)

# CV
Temp.CV.Toff.lm.red<-lme(log.T.off~ CV.s, random = ~1|Nest.ID,  data = dur.on.red.dat, na.action = na.omit) 
summary(Temp.CV.Toff.lm.red)
0.2517153^2/(0.2517153^2 + 0.2204819^2)
Toff.CV.pred.dat.red<-dur.on.red.dat%>% filter(!is.na(CV.s))
Toff.CV.pred.dat.red<-Toff.CV.pred.dat.red%>%filter(!is.na(log.T.off))
Toff.CV.pred.dat.red$pred<-predict(Temp.CV.Toff.lm.red, type="response") #
Toff.CV.pred.dat.red$exp.pred<- exp(Toff.CV.pred.dat.red$pred) #exponent to get back on scale
Toff.bouts.TCV.red.fig<-ggplot(Toff.CV.pred.dat.red, aes(CV, T.off))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= exp.pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest
  geom_line(aes(y = exp(predict(Temp.CV.Toff.lm.red, level = 0))), color = "black",size = 1.1)+ # pop line
  #scale_y_continuous(breaks=seq(0,25, 5), limits=c(0,25))+
  #scale_x_continuous(breaks=seq(12,42, 4), limits=c(12,42))+
  labs(x = "Variation (CV) in ambient temperature", y = "Off-bout duration (min)")+
  ggtitle("(D)")
Toff.bouts.TCV.red.fig
#ggsave(path = "Figs.OnBouts", "Toff.bouts.TCV.fig.jpg", Toff.bouts.TCV.fig, units = "in", height=4, width = 4, dpi = 600)

#### Percent time on nest and ambient T ####

hist(bout_amb_dat.red$perc.in) # range 
summary(bout_amb_dat.red$perc.in)
describeBy(bout_amb_dat.red$perc.in)
#Tmin
#deleted code with % in since wrong distribution. use proportion and beta dist

#try with beta dist - make proportion but some 1 so - 0.01
bout_amb_dat.red$p.in<-bout_amb_dat.red$perc.in/100 - 0.01
summary(bout_amb_dat.red$p.in)
describeBy(bout_amb_dat.red$p.in)
which(bout_amb_dat.red$p.in ==1)
bout_amb_age_dat<-bout_amb_dat.red %>% filter(!F.Age == c("AHY"))
levels(factor(bout_amb_age_dat$F.Age))

#age prop in & age
p.in.agelm<-glmmTMB(p.in~F.Age + (1|Nest.ID), family=beta_family(), data=bout_amb_dat.red)
summary(p.in.agelm)
pairs(emmeans(p.in.agelm, ~F.Age))
# use bout_amb_age_dat Jul 2024
p.in.agelm2<-glmmTMB(p.in~F.Age + (1|Nest.ID), family=beta_family(), data=bout_amb_age_dat)
summary(p.in.agelm2) # 0.068
describeBy(bout_amb_age_dat$p.in, bout_amb_age_dat$F.Age)


# prop in and temps
prop.in.Tmin.beta<-glmmTMB(p.in~T.min.s + (1|Nest.ID), family=beta_family(), data=bout_amb_dat.red)
summary(prop.in.Tmin.beta)

Tmin.prop.pred.dat<-bout_amb_dat.red%>% filter(!is.na(T.min.s))
Tmin.prop.pred.dat$pred<-predict(prop.in.Tmin.beta, type="response") #
Tmin.prop.pred.dat$pred.pop<-predict(prop.in.Tmin.beta, type = "response", re.form = NA)
Tmin.prop.fig<-ggplot(Tmin.prop.pred.dat, aes(T.min, p.in))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest
  geom_line(aes(y = pred.pop), color = "black",size = 1.1)+ # pop line
  #scale_y_continuous(breaks=seq(0,200, 10), limits=c(0,200))+
  #scale_x_continuous(breaks=seq(10,30, 5), limits=c(10,30))+
  labs(x = "Minimum ambient temperature (\u00B0C)", y = "Proportion of time on the nest")+
  ggtitle("(A)")
Tmin.prop.fig
ggsave(path = "Figs.OnBouts", "Tmin.prop.fig.jpg", Tmin.prop.fig, units = "in", height=4, width = 4, dpi = 600)

#Tmax
#deleted code with % in since wrong distribution. use proportion and beta dist

prop.in.Tmax.beta<-glmmTMB(p.in~T.max.s + (1|Nest.ID), family=beta_family(), data=bout_amb_dat.red)
summary(prop.in.Tmax.beta)

#Tmean
#deleted code with % in since wrong distribution. use proportion and beta dist
prop.in.Tmean.beta<-glmmTMB(p.in~T.mean.s + (1|Nest.ID), family=beta_family(), data=bout_amb_dat.red)
summary(prop.in.Tmean.beta)

#CV
#deleted code with % in since wrong distribution. use proportion and beta dist

prop.in.TCV.beta<-glmmTMB(p.in~CV.s + (1|Nest.ID), family=beta_family(), data=bout_amb_dat.red)
summary(prop.in.TCV.beta)
TCV.prop.pred.dat<-bout_amb_dat.red%>% filter(!is.na(CV.s))
TCV.prop.pred.dat$pred<-predict(prop.in.TCV.beta, type="response") #
TCV.prop.pred.dat$pred.pop<-predict(prop.in.TCV.beta, type = "response", re.form = NA)
TCV.prop.fig<-ggplot(TCV.prop.pred.dat, aes(CV, p.in))+
  theme_classic()+
  geom_jitter(color="gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest
  geom_line(aes(y = pred.pop), color = "black",size = 1.1)+ # pop line
  #scale_y_continuous(breaks=seq(0,200, 10), limits=c(0,200))+
  #scale_x_continuous(breaks=seq(10,30, 5), limits=c(10,30))+
  labs(x = "Variation (CV) in ambient temperature", y = "Proportion of time on the nest")+
  ggtitle("(B)")
TCV.prop.fig
ggsave(path = "Figs.OnBouts", "TCV.prop.fig.jpg", TCV.prop.fig, units = "in", height=4, width = 4, dpi = 600)


#### Are there diffs in incubation behavior by Age ####
#subset out AHY
length(which(bout_amb_dat$F.Age == "AHY"))
bout_amb_dat[c(170:179, 357:361), 1]
bout_age_dat<-bout_amb_dat %>% filter(!F.Age == c("AHY"))
levels(factor(bout_age_dat$F.Age))
## order age
bout_age_dat$F.Age<-factor(bout_age_dat$F.Age, levels = c("SY", "ASY"))

# check temps by age
describeBy(bout_age_dat$T.min, bout_age_dat$F.Age)
min.age.lm<-lme(T.min.s~F.Age,random = ~1|Nest.ID, data = bout_age_dat, na.action = na.omit)
summary(min.age.lm)


#N on. redid with on.per.day
n.on.age.lm<-lme(on.per.day~F.Age,random = ~1|Nest.ID, data = bout_age_dat, na.action = na.omit)
summary(n.on.age.lm)
plot(n.on.age.lm)                           

#N off
n.off.age.lm<-lme(off.per.day~F.Age,random = ~1|Nest.ID, data = bout_age_dat, na.action = na.omit)
summary(n.off.age.lm)

# T.on
T.on.age.lm<-lme(log.T.on~ F.Age, random = ~1|Nest.ID, data= bout_age_dat, na.action = na.omit)
summary(T.on.age.lm)
plot(T.on.age.lm)

T.on.age.tpi.lm<-lme(log.T.on~ F.Age + tpi50_3cat, random = ~1|Nest.ID, data= bout_age_dat, na.action = na.omit)
summary(T.on.age.tpi.lm)

table(bout_amb_dat$F.Age, bout_amb_dat$tpi50_6cat)

  #log scale so use describeBy instead of emmean since no other factors 
T.on.desc<-describeBy(bout_age_dat$T.on, bout_age_dat$F.Age, mat = TRUE)
T.on.desc
T.on.desc$group1<-factor(T.on.desc$group1, levels = c("SY", "ASY"))

library(ggbeeswarm)
Age.T.on.fig<-ggplot(T.on.desc, aes(group1, mean))+
  theme_classic()+
  geom_beeswarm(bout_age_dat, mapping= aes(F.Age, T.on), color = "gray", alpha = 0.5, size = 0.8)+
  geom_point(aes(y=mean), size = 2)+
  #geom_beeswarm(bout_age_dat, mapping = aes(F.Age, T.on), color = "grey")+
  geom_errorbar(aes(ymin = mean - se, ymax=mean+se))+
  labs(x = "Female age", y = "Duration of on-bout (min)")+
  scale_y_continuous(breaks=seq(0,200, 20), limits=c(0,200))
Age.T.on.fig
ggsave(path = "Figs.OnBouts", "Age.T.on.fig.jpg", Age.T.on.fig, units = "in", height=4, width = 4, dpi = 600)

# T.off
T.off.age.lm<-lme(log.T.off~ F.Age, random = ~1|Nest.ID, data= bout_age_dat, na.action = na.omit)
summary(T.off.age.lm)

#### does incubation behavior affect nest survival? ####
levels(factor(bout_amb_dat$survive))
#### scale incubation behaviors ####
bout_amb_dat$number.on.bouts.s<-scale(bout_amb_dat$number.on.bouts)
bout_amb_dat$number.off.bouts.s<-scale(bout_amb_dat$number.off.bouts)
bout_amb_dat$on.per.day.s<-scale(bout_amb_dat$on.per.day)
bout_amb_dat$off.per.day.s<-scale(bout_amb_dat$off.per.day)
bout_amb_dat$T.on.s<-scale(bout_amb_dat$T.on)
bout_amb_dat$T.off.s<-scale(bout_amb_dat$T.off)
bout_amb_dat$perc.in.s<-scale(bout_amb_dat$perc.in)

# start with on bouts  # re did with on per day
surv.age.lm<-glmer(survive ~ on.per.day.s + T.on.s  + perc.in.s +F.Age+ expos +(1|Nest.ID), family = binomial, data = bout_amb_dat, na.action = na.omit) # wouldn't converge w/T.off - too related?

surv.lm<-glmer(survive ~ on.per.day.s + T.on.s  + perc.in.s +expos +(1|Nest.ID), family = binomial, data = bout_amb_dat, na.action = na.omit) #
summary(surv.lm) # none sig
vif(surv.lm)
surv.lm1<-glmer(survive ~ on.per.day.s + T.on.s  + expos + (1|Nest.ID), family = binomial, data = bout_amb_dat, na.action = na.omit) 
summary(surv.lm1)
surv.lm2<-glmer(survive ~ on.per.day.s  + perc.in.s + expos + (1|Nest.ID), family = binomial, data = bout_amb_dat, na.action = na.omit) 
summary(surv.lm2)
surv.lm3<-glmer(survive ~  T.on.s  + perc.in.s + expos + (1|Nest.ID), family = binomial, data = bout_amb_dat, na.action = na.omit) 
summary(surv.lm3)
surv.lm4<-glmer(survive ~ on.per.day.s  + expos+ (1|Nest.ID), family = binomial, data = bout_amb_dat, na.action = na.omit) 
surv.lm5<-glmer(survive ~  T.on.s + expos +  (1|Nest.ID), family = binomial, data = bout_amb_dat, na.action = na.omit)
summary(surv.lm5)
surv.lm6<-glmer(survive ~  perc.in.s + expos + (1|Nest.ID), family = binomial, data = bout_amb_dat, na.action = na.omit) 

model.sel(surv.lm, surv.lm1, surv.lm2, surv.lm3, surv.lm4, surv.lm5, surv.lm6)
summary(surv.lm6) 


# off bouts? # redid with off.per.day.s
surv.off.lm<-glmer(survive ~  T.off.s + off.per.day.s +expos + (1|Nest.ID), family = binomial, data = bout_amb_dat, na.action = na.omit) 
summary(surv.off.lm)
vif(surv.off.lm)
surv.t.off.lm<-glmer(survive ~  T.off.s + expos + (1|Nest.ID), family = binomial, data = bout_amb_dat, na.action = na.omit) 
summary(surv.t.off.lm)
surv.n.off.bouts.lm<-glmer(survive ~ off.per.day.s +  expos  + (1|Nest.ID), family = binomial, data = bout_amb_dat, na.action = na.omit) 
summary(surv.n.off.bouts.lm) # without expos converges but not sig.
summary(bout_amb_dat$number.off.bouts.s)
hist(bout_amb_dat$number.off.bouts.s)
vif(surv.n.off.bouts.lm)

#### fledged per nest and female behave ####
fledged.lm <- glmer(fledged ~ on.per.day.s + T.on.s  + perc.in.s + expos +(1|Nest.ID), family = poisson, data = bout_amb_dat)
vif(fledged.lm)
summary(fledged.lm) 

#### Young per successful nest ####
suc.nests.dat<-bout_amb_dat[bout_amb_dat$survive ==0,]
summary(suc.nests.dat$fledged)
describeBy(suc.nests.dat$fledged) # but repeats need unique nests
nests.uniq.dat<- distinct(suc.nests.dat, Nest.ID, .keep_all = TRUE)
head(nests.uniq.dat)
describeBy(nests.uniq.dat$fledged)


fledged.inc.lm <- glmer(fledged ~ on.per.day.s + T.on.s  + perc.in.s + expos +(1|Nest.ID), family = poisson, data = suc.nests.dat)
vif(fledged.inc.lm)
summary(fledged.inc.lm) # failed to converge
fledged.inc.lm1 <- glmer(fledged ~ on.per.day.s + T.on.s  +  expos +(1|Nest.ID), family = poisson, data = suc.nests.dat)
summary(fledged.inc.lm1) # reported this since last with perc.in alone did not converge 
fledged.inc.lm2 <- glmer(fledged ~ on.per.day.s +  perc.in.s + expos +(1|Nest.ID), family = poisson, data = suc.nests.dat)
fledged.inc.lm3 <- glmer(fledged ~ T.on.s  + perc.in.s + expos +(1|Nest.ID), family = poisson, data = suc.nests.dat)
fledged.inc.lm4 <- glmer(fledged ~ on.per.day.s  + expos +(1|Nest.ID), family = poisson, data = suc.nests.dat)
fledged.inc.lm5 <- glmer(fledged ~  T.on.s   + expos +(1|Nest.ID), family = poisson, data = suc.nests.dat)
fledged.inc.lm6 <- glmer(fledged ~  perc.in.s + expos +(1|Nest.ID), family = poisson, data = suc.nests.dat)
summary(fledged.inc.lm6)
model.sel(fledged.inc.lm, fledged.inc.lm1, fledged.inc.lm2, fledged.inc.lm3, fledged.inc.lm4, fledged.inc.lm5, fledged.inc.lm6)


fledged.inc.off.lm<- glmer(fledged ~ off.per.day.s + T.off.s  + expos +(1|Nest.ID), family = poisson, data = suc.nests.dat)
vif(fledged.inc.off.lm)
summary(fledged.inc.off.lm) # failed to converge do separate
fledged.inc.off.lm1<- glmer(fledged ~ off.per.day.s   + expos +(1|Nest.ID), family = poisson, data = suc.nests.dat)
summary(fledged.inc.off.lm1) # failed to converge
fledged.inc.off.lm2<- glmer(fledged ~  T.off.s  + expos +(1|Nest.ID), family = poisson, data = suc.nests.dat)
summary(fledged.inc.off.lm2)

#### Does incubation behavior affect Nest Temps ####
bout_nest_dat<-bout_dat %>% filter(Type == "nest") # n = 418
summary(bout_nest_dat$T.min)
summary(bout_nest_dat$T.max)
summary(bout_nest_dat$T.mean)
#### scale incubation behaviors ####
bout_nest_dat$number.on.bouts.s<-scale(bout_nest_dat$number.on.bouts)
bout_nest_dat$number.off.bouts.s<-scale(bout_nest_dat$number.off.bouts)
bout_nest_dat$T.on.s<-scale(bout_nest_dat$T.on)
bout_nest_dat$T.off.s<-scale(bout_nest_dat$T.off)
bout_nest_dat$perc.in.s<-scale(bout_nest_dat$perc.in)
bout_nest_dat$on.per.day.s<-scale(bout_nest_dat$on.per.day)
bout_nest_dat$off.per.day.s<-scale(bout_nest_dat$off.per.day)
summary(bout_nest_dat$on.per.day.s)
summary(bout_nest_dat$off.per.day.s)
which(bout_nest_dat$T.on >120)
which(bout_age_dat$T.off >30)
#### reduced data set for duration of on and off bout analyses of nest temps ####
dur.nest.red.dat<-filter(bout_nest_dat, T.on < 120 & T.off < 30) 
write.csv(dur.nest.red.dat, "dur.nest.red.dat.csv",  row.names = FALSE)
write.csv(bout_amb_dat, "bout_amb_dat.csv", row.names = FALSE)
#### T.min and inc behav ####
describeBy(bout_nest_dat$T.min)
#redid on and off bouts with on.per.day.s or off.per.day.s
T.min.nest.lm<- lme(T.min~ on.per.day.s + T.on.s + perc.in.s  + T.off.s,  random = ~1|Nest.ID, method= "ML", data = bout_nest_dat, na.action = na.omit)
vif(T.min.nest.lm) # can't use N off bouts vif = 1219 
#Top model when all combos. BUT behavs are not independent. Run separate
summary(T.min.nest.lm)

# number on
T.min.n.on.lm<- lme(T.min~ on.per.day.s,  random = ~1|Nest.ID, method= "REML", data = bout_nest_dat, na.action = na.omit)
summary(T.min.n.on.lm)
#Jul 2024 nest.dat.age
T.min.n.on.age.lm<- lme(T.min~ on.per.day.s + F.Age,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(T.min.n.on.age.lm)

Tmin.n.on.pred.dat<-bout_nest_dat%>% filter(!is.na(T.min))
Tmin.n.on.pred.dat<-Tmin.n.on.pred.dat%>%filter(!is.na(on.per.day.s))
Tmin.n.on.pred.dat$pred<-predict(T.min.n.on.lm, type="response") #
Tmin.n.on.fig<-ggplot(Tmin.n.on.pred.dat, aes(on.per.day, T.min))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest     
  geom_line(aes(y = predict(T.min.n.on.lm, level = 0)), color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(12,38, 4), limits=c(12,38))+
  scale_x_continuous(breaks=seq(0,60, 10), limits=c(0,60))+
  labs(y = "Minimum nest temperature (\u00B0C)", x = "On-bouts per day")+
  ggtitle("(A)")
Tmin.n.on.fig
ggsave(path = "Figs.OnBouts", "Tmin.n.on.fig.jpg", Tmin.n.on.fig, units = "in", height=4, width = 4, dpi = 600)

# number off
T.min.n.off.lm<- lme(T.min~ off.per.day.s,  random = ~1|Nest.ID, method= "REML", data = bout_nest_dat, na.action = na.omit)
summary(T.min.n.off.lm)
describeBy(bout_nest_dat$off.per.day)
Tmin.n.off.pred.dat<-bout_nest_dat%>% filter(!is.na(T.min))
Tmin.n.off.pred.dat<-Tmin.n.off.pred.dat%>%filter(!is.na(off.per.day.s))
Tmin.n.off.pred.dat$pred<-predict(T.min.n.off.lm, type="response") #
Tmin.n.off.fig<-ggplot(Tmin.n.off.pred.dat, aes(off.per.day, T.min))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest     
  geom_smooth(aes(y = predict(T.min.n.on.lm, level = 0)), method = "lm", fill = NA, color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(12,38, 4), limits=c(12,38))+
  scale_x_continuous(breaks=seq(0,60, 10), limits=c(0,60))+
  labs(y = "Minimum nest temperature (\u00B0C)", x = "Off-bouts per day")+
  ggtitle("(A)")
Tmin.n.off.fig
ggsave(path = "Figs.OnBouts", "Tmin.n.off.fig.jpg", Tmin.n.off.fig, units = "in", height=4, width = 4, dpi = 600)

# on bout duration NB: reran with dur.nest.red.dat 
T.min.T.on.lm<- lme(T.min~ T.on.s,  random = ~1|Nest.ID, method= "REML", data = dur.nest.red.dat, na.action = na.omit)
summary(T.min.T.on.lm) # 
plot(T.min.T.on.lm)

Tmin.dur.on.pred.dat<-dur.nest.red.dat%>% filter(!is.na(T.min))
Tmin.dur.on.pred.dat$pred<-predict(T.min.T.on.lm, type="response") #
Tmin.dur.on.fig<-ggplot(Tmin.dur.on.pred.dat, aes(T.on, T.min))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest 
  geom_line(aes(y = predict(T.min.T.on.lm, level = 0)), color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(12,38, 4), limits=c(12,38))+
  #scale_x_continuous(breaks=seq(0,60, 10), limits=c(0,60))+
  labs(y = "Minimum nest temperature (\u00B0C)", x = "On-bout duration (min)")+
  ggtitle("(B)")
Tmin.dur.on.fig
ggsave(path = "Figs.OnBouts", "Tmin.dur.on.fig.jpg", Tmin.dur.on.fig, units = "in", height=4, width = 4, dpi = 600)

# off bout duration
T.min.T.off.lm<- lme(T.min~ T.off.s,  random = ~1|Nest.ID, method= "REML", data = dur.nest.red.dat, na.action = na.omit)
summary(T.min.T.off.lm)
Tmin.T.off.pred.dat<-dur.nest.red.dat%>% filter(!is.na(T.min))
Tmin.T.off.pred.dat<-Tmin.T.off.pred.dat%>% filter(!is.na(T.off.s))
Tmin.T.off.pred.dat$pred<-predict(T.min.T.off.lm, type="response") #
#Tmin.T.off.pred.dat$pred.pop<-predict(T.min.T.off.lm, type = "response", re.form = NA)
Tmin.T.off.fig<-ggplot(Tmin.T.off.pred.dat, aes(T.off, T.min))+
  theme_classic()+
  geom_jitter(color="gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.5, linetype = "dashed" )+ # add by Nest   
  geom_line(aes(y = predict(T.min.T.off.lm, level = 0)), color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(12,38, 4), limits=c(12,38))+
  #scale_x_continuous(breaks=seq(0,60, 10), limits=c(0,60))+
  labs(y = "Minimum nest Temperature (\u00B0C)", x = "Off-bout duration (min)")+
  ggtitle("(C)")
Tmin.T.off.fig
ggsave(path = "Figs.OnBouts", "Tmin.T.off.fig.jpg", Tmin.T.off.fig, units = "in", height=4, width = 6, dpi = 300)
 
bout_nest_dat.red<-bout_nest_dat %>% filter(!Nest.ID == c("321.12"))   
T.min.P.on.lm<- lme(T.min~ perc.in.s,  random = ~1|Nest.ID, method= "REML", data = bout_nest_dat.red, na.action = na.omit)
summary(T.min.P.on.lm)
Tmin.P.on.pred.dat<-bout_nest_dat.red%>% filter(!is.na(T.min))
Tmin.P.on.pred.dat$pred<-predict(T.min.P.on.lm, type="response")
P.on.T.min.fig<-ggplot(Tmin.P.on.pred.dat, aes(perc.in, T.min))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest
  geom_line(aes(y = predict(T.min.P.on.lm, level = 0)), color = "black",size = 1.1)+ # pop line
  labs(y = "Minimum nest temperature (\u00B0C)", x = "Percent time on the nest")+
  scale_y_continuous(breaks=seq(12,38, 4), limits=c(12,38))+
  ggtitle("(D)")
P.on.T.min.fig
ggsave(path = "Figs.OnBouts", "P.on.T.min.fig.jpg", P.on.T.min.fig, units = "in", height=4, width = 4, dpi = 600)

#### T.max and inc behav ####
describeBy(bout_nest_dat$T.max)
hist(bout_nest_dat$T.max)
which(bout_nest_dat$T.max >45)
T.max.dat<-filter(bout_nest_dat, T.max < 45) # remove bad data. 
describeBy(T.max.dat$T.max)
T.max.dur.red.dat<-filter(T.max.dat, T.on < 120 & T.off < 30)
describeBy(T.max.dur.red.dat$T.max)
# number on
T.max.n.on.lm<- lme(T.max~ on.per.day.s,  random = ~1|Nest.ID, method= "REML", data = T.max.dat, na.action = na.omit)
summary(T.max.n.on.lm)

# off
T.max.n.off.lm<- lme(T.max~ off.per.day.s,  random = ~1|Nest.ID, method= "REML", data = T.max.dat, na.action = na.omit)
summary(T.max.n.off.lm)

# on bout duration
T.max.T.on.lm<- lme(T.max~ T.on.s,  random = ~1|Nest.ID, method= "REML", data = T.max.dur.red.dat , na.action = na.omit)
summary(T.max.T.on.lm) # 

# off bout duration
T.max.T.off.lm<- lme(T.max~ T.off.s,  random = ~1|Nest.ID, method= "REML", data = T.max.dur.red.dat, na.action = na.omit)
summary(T.max.T.off.lm)
Tmax.T.off.pred.dat<-dur.nest.red.dat%>% filter(!is.na(T.max))
Tmax.T.off.pred.dat<-Tmax.T.off.pred.dat%>% filter(!is.na(T.off.s))
Tmax.T.off.pred.dat$pred<-predict(T.max.T.off.lm, type="response") #

Tmax.n.on.fig<-ggplot(Tmax.T.off.pred.dat, aes(T.off, T.max))+
  theme_classic()+
  geom_point()+
  geom_line(aes(y = predict(T.max.T.off.lm, level = 0)), color = "black",size = 1.1)+ # pop line
  geom_line(aes(y= pred, group = Nest.ID), color = "gray60", size = 0.5, linetype = "dashed" )+ # add by Nest   
  #scale_y_continuous(breaks=seq(25,40, 5), limits=c(25,40))+
  #scale_x_continuous(breaks=seq(0,60, 10), limits=c(0,60))+
  labs(y = "Maximum nest temperature (\u00B0C)", x = "Duration of off-bout (min)")
Tmax.n.on.fig
ggsave(path = "Figs.OnBouts", "Tmin.n.on.fig.jpg", Tmin.n.on.fig, units = "in", height=4, width = 6, dpi = 300)

 
T.max.P.on.lm<- lme(T.max~ perc.in.s,  random = ~1|Nest.ID, method= "REML", data = T.max.dat, na.action = na.omit)
summary(T.max.P.on.lm) 

#### T.mean and inc behavs ####
summary(bout_nest_dat$T.mean)
describeBy(bout_nest_dat$T.mean)

# number on
T.mean.n.on.lm<- lme(T.mean~ on.per.day.s,  random = ~1|Nest.ID, method= "REML", data = bout_nest_dat, na.action = na.omit)
summary(T.mean.n.on.lm)

Tmean.n.on.pred.dat<-bout_nest_dat%>% filter(!is.na(T.mean))
Tmean.n.on.pred.dat<-Tmean.n.on.pred.dat%>%filter(!is.na(on.per.day))
Tmean.n.on.pred.dat$pred<-predict(T.mean.n.on.lm, type="response") #
Tmean.n.on.fig<-ggplot(Tmean.n.on.pred.dat, aes(on.per.day, T.mean))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest
  geom_line(aes(y = predict(T.mean.n.on.lm, level = 0)), color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(20,40, 5), limits=c(22,40))+
  #scale_x_continuous(breaks=seq(0,60, 10), limits=c(0,60))+
  labs(y = "Mean nest temperature (\u00B0C)",x = "On-bouts per day")+
  ggtitle("(A)")
Tmean.n.on.fig
ggsave(path = "Figs.OnBouts", "Tmean.n.on.fig.jpg", Tmean.n.on.fig, units = "in", height=4, width = 4, dpi = 600)

# number off
T.mean.n.off.lm<- lme(T.mean~ off.per.day.s,  random = ~1|Nest.ID, method= "REML", data = bout_nest_dat, na.action = na.omit)
summary(T.mean.n.off.lm)

Tmean.n.off.pred.dat<-bout_nest_dat%>% filter(!is.na(T.mean))
Tmean.n.off.pred.dat<-Tmean.n.off.pred.dat%>%filter(!is.na(off.per.day))
Tmean.n.off.pred.dat$pred<-predict(T.mean.n.off.lm, type="response") #
Tmean.n.off.fig<-ggplot(Tmean.n.off.pred.dat, aes(off.per.day, T.mean))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest
  geom_line(aes(y = predict(T.mean.n.off.lm, level = 0)), color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(20,40, 5), limits=c(22,40))+
  #scale_x_continuous(breaks=seq(0,60, 10), limits=c(0,60))+
  labs(y = "Mean nest temperature (\u00B0C)",x = "Off-bouts per day")+
  ggtitle("(A)")
Tmean.n.off.fig
ggsave(path = "Figs.OnBouts", "Tmean.n.off.fig.jpg", Tmean.n.off.fig, units = "in", height=4, width = 4, dpi = 600)

# on bout duration
T.mean.T.on.lm<- lme(T.mean~ T.on.s,  random = ~1|Nest.ID, method= "REML", data = dur.nest.red.dat, na.action = na.omit)
summary(T.mean.T.on.lm) # 

Tmean.T.on.pred.dat<-dur.nest.red.dat%>% filter(!is.na(T.mean))
Tmean.T.on.pred.dat$pred<-predict(T.mean.T.on.lm, type="response") #
Tmean.T.on.fig<-ggplot(Tmean.T.on.pred.dat, aes(T.on, T.mean))+
  theme_classic()+
  geom_jitter(color= "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest 
  geom_line(aes(y = predict(T.mean.T.on.lm, level = 0)), color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(20,40, 5), limits=c(22,40))+
  #scale_x_continuous(breaks=seq(0,60, 10), limits=c(0,60))+
  labs(y = "Mean nest temperature (\u00B0C)",x = "On-bout duration (min)")+
  ggtitle("(B)")
Tmean.T.on.fig

# off bout duration
T.mean.T.off.lm<- lme(T.mean~ T.off.s,  random = ~1|Nest.ID, method= "REML", data = dur.nest.red.dat, na.action = na.omit)
summary(T.mean.T.off.lm)
Tmean.T.off.pred.dat<-dur.nest.red.dat%>% filter(!is.na(T.mean))
Tmean.T.off.pred.dat<-Tmean.T.off.pred.dat%>% filter(!is.na(T.off.s))
Tmean.T.off.pred.dat$pred<-predict(T.mean.T.off.lm, type="response") #
nest.MeanT.T.off.fig<-ggplot(Tmean.T.off.pred.dat, aes(T.off, T.mean))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest
  geom_line(aes(y = predict(T.mean.T.off.lm, level = 0)), color = "black",size = 1.1)+ # pop line
  scale_y_continuous(breaks=seq(20,40, 5), limits=c(22,40))+
  #scale_x_continuous(breaks=seq(0,60, 10), limits=c(0,60))+
  labs(y = "Mean nest temperature (\u00B0C)", x = "Off-bout duration (min)")+
  ggtitle("(C)")
nest.MeanT.T.off.fig


   
T.mean.P.on.lm<- lme(T.mean~ perc.in.s,  random = ~1|Nest.ID, method= "REML", data = bout_nest_dat.red, na.action = na.omit)
summary(T.mean.P.on.lm)
Tmean.P.on.pred.dat<-bout_nest_dat.red%>% filter(!is.na(T.mean))
Tmean.P.on.pred.dat$pred<-predict(T.mean.P.on.lm, type="response",re.form=NA )
P.on.T.mean.fig<-ggplot(Tmean.P.on.pred.dat, aes(perc.in, T.mean))+
  theme_classic()+
  geom_jitter(color="gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest
  geom_line(aes(y = predict(T.mean.P.on.lm, level = 0)), color = "black",size = 1.1)+ # pop line
  labs(y = "Mean nest temperature (\u00B0C)", x = "Percent time on the nest")+
  scale_y_continuous(breaks=seq(20,40, 5), limits=c(22,40))+
  ggtitle("(D)")
P.on.T.mean.fig
ggsave(path = "Figs.OnBouts", "P.on.T.mean.fig.jpg", P.on.T.mean.fig, units = "in", height=4, width = 6, dpi = 300)

#### CV and inc behavs ####
summary(bout_nest_dat$CV)
describeBy(bout_nest_dat$CV)

# number on
CV.n.on.lm<- lme(CV~ on.per.day.s,  random = ~1|Nest.ID, method= "REML", data = bout_nest_dat, na.action = na.omit)
summary(CV.n.on.lm)

CV.n.on.pred.dat<-bout_nest_dat%>% filter(!is.na(CV))
CV.n.on.pred.dat<-CV.n.on.pred.dat%>% filter(!is.na(on.per.day.s))
CV.n.on.pred.dat$pred<-predict(CV.n.on.lm, type="response") #
CV.n.on.fig<-ggplot(CV.n.on.pred.dat, aes(on.per.day, CV))+
  theme_classic()+
  geom_jitter(color="gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest
  geom_line(aes(y = predict(CV.n.on.lm, level = 0)), color = "black",size = 1.1)+ # pop line
  #scale_y_continuous(breaks=seq(10,40, 5), limits=c(10,40))+
  #scale_x_continuous(breaks=seq(0,60, 10), limits=c(0,60))+
  labs(y = "Variation (CV) in nest temperature",x = "Number of on-bouts")+
  ggtitle("(A)")
CV.n.on.fig
#ggsave(path = "Figs.OnBouts", "CV.n.on.fig.jpg", CV.n.on.fig, units = "in", height=4, width = 4, dpi = 600)

# number off
CV.n.off.lm<- lme(CV~ off.per.day.s,  random = ~1|Nest.ID, method= "REML", data = bout_nest_dat, na.action = na.omit)
summary(CV.n.off.lm)

CV.n.off.pred.dat<-bout_nest_dat%>% filter(!is.na(CV))
CV.n.off.pred.dat<-CV.n.off.pred.dat%>% filter(!is.na(off.per.day.s))
CV.n.off.pred.dat$pred<-predict(CV.n.off.lm, type="response") #
CV.n.off.fig<-ggplot(CV.n.off.pred.dat, aes(off.per.day, CV))+
  theme_classic()+
  geom_jitter(color="gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest
  geom_line(aes(y = predict(CV.n.off.lm, level = 0)), color = "black",size = 1.1)+ # pop line
  #scale_y_continuous(breaks=seq(10,40, 5), limits=c(10,40))+
  #scale_x_continuous(breaks=seq(0,60, 10), limits=c(0,60))+
  labs(y = "Variation (CV) in nest temperature",x = "Off-bouts per day")+
  ggtitle("(A)")
CV.n.off.fig
ggsave(path = "Figs.OnBouts", "CV.n.off.fig.jpg", CV.n.off.fig, units = "in", height=4, width = 4, dpi = 600)

# on bout duration
CV.T.on.lm<- lme(CV~ T.on.s,  random = ~1|Nest.ID, method= "REML", data = dur.nest.red.dat, na.action = na.omit)
summary(CV.T.on.lm) # 

CV.T.on.pred.dat<-dur.nest.red.dat%>% filter(!is.na(CV))
CV.T.on.pred.dat<-CV.T.on.pred.dat%>% filter(!is.na(T.on.s))
CV.T.on.pred.dat$pred<-predict(CV.T.on.lm, type="response") #
CV.T.on.fig<-ggplot(CV.T.on.pred.dat, aes(T.on, CV))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest 
  geom_line(aes(y = predict(CV.T.on.lm, level = 0)), color = "black", size = 1.1)+ # pop line
  #scale_y_continuous(breaks=seq(10,40, 5), limits=c(10,40))+
  #scale_x_continuous(breaks=seq(0,60, 10), limits=c(0,60))+
  labs(y = "Variation (CV) in nest temperature", x = "On-bout duration (min)")+
  ggtitle("(B)")
CV.T.on.fig

# off bout duration
CV.T.off.lm<- lme(CV~ T.off.s,  random = ~1|Nest.ID, method= "REML", data = dur.nest.red.dat, na.action = na.omit)
summary(CV.T.off.lm)
CV.T.off.pred.dat<-dur.nest.red.dat%>% filter(!is.na(CV))
CV.T.off.pred.dat<-CV.T.off.pred.dat%>% filter(!is.na(T.off.s))
CV.T.off.pred.dat$pred<-predict(CV.T.off.lm, type="response") #
CV.T.off.fig<-ggplot(CV.T.off.pred.dat, aes(T.off, CV))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = 0.4, linetype = "dashed" )+ # add by Nest   
  geom_line(aes(y = predict(CV.T.off.lm, level = 0)), color = "black", size = 1.1)+ # pop line
#scale_y_continuous(breaks=seq(10,40, 5), limits=c(10,40))+
#scale_x_continuous(breaks=seq(0,60, 10), limits=c(0,60))+
 labs(y = "Variation (CV) in nest temperature", x = "Off-bout duration (min)")+
  ggtitle("(C)")
CV.T.off.fig
ggsave(path = "Figs.OnBouts", "CV.T.off.fig.jpg", CV.T.off.fig, units = "in", height=4, width = 6, dpi = 300)


CV.P.on.lm<- lme(CV~ perc.in.s,  random = ~1|Nest.ID, method= "REML", data = bout_nest_dat.red, na.action = na.omit)
summary(CV.P.on.lm)
CV.P.on.pred.dat<-bout_nest_dat.red%>% filter(!is.na(CV))
CV.P.on.pred.dat$pred<-predict(CV.P.on.lm, type="response")
P.on.CV.fig<-ggplot(CV.P.on.pred.dat, aes(perc.in, CV))+
  theme_classic()+
  geom_jitter(color = "gray40", alpha = 0.2, size = 0.8)+
  geom_line(aes(y= pred, group = Nest.ID), color = "gray29", size = .4, linetype = "dashed" )+ # add by Nest
  geom_line(aes(y = predict(CV.P.on.lm, level = 0)), color = "black",size = 1.1)+ # pop line
  #geom_smooth(aes(x = perc.in, y= pred), method = "lm", color = "black")+
  labs(y = "Variation (CV) in nest temperature", x = "Percent time on the nest")+
  ggtitle("(D)")
P.on.CV.fig
ggsave(path = "Figs.OnBouts", "P.on.CV.fig.jpg", P.on.CV.fig, units = "in", height=4, width = 4, dpi = 600)

#### incubation temps and female age ####
nest.dat.age<-bout_nest_dat %>% filter(!F.Age == c("AHY"))
names(nest.dat.age)
describeBy(nest.dat.age$T.min, nest.dat.age$F.Age)
minT.age.lm<- lme(T.min~ F.Age,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(minT.age.lm)
describeBy(nest.dat.age$T.max, nest.dat.age$F.Age)
maxT.age.lm<- lme(T.max~ F.Age,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(maxT.age.lm)
meanT.age.lm<- lme(T.mean~ F.Age,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(meanT.age.lm)
maxT.age.lm<- lme(T.max~ F.Age,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(maxT.age.lm)
cv.age.lm<- lme(CV~ F.Age,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(cv.age.lm)

#### July 2024. check models of inc behav temps by age ####
CV.T.on.age.lm<- lme(CV~ T.on.s + F.Age,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(CV.T.on.age.lm) # 

minT.dur.age.lm<- lme(T.min~ F.Age + T.on.s,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(minT.dur.age.lm)

Tmean.dur.age.lm<- lme(T.mean~ F.Age + T.on.s,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(Tmean.dur.age.lm)

Tmax.dur.age.lm<- lme(T.max~ F.Age + T.on.s,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(Tmax.dur.age.lm)

#### Apr 2024 first and last off bouts ####
describeBy(bout_nest_dat$first.off, list(bout_nest_dat$tpi50_3cat, bout_nest_dat$F.Age))
str(bout_nest_dat$first.off)
describeBy(bout_nest_dat$last.on,list(bout_nest_dat$tpi50_3cat, bout_nest_dat$F.Age)) 

describeBy(nest.dat.age$first.off, list(nest.dat.age$tpi50_3cat, nest.dat.age$F.Age), mat = TRUE)
describeBy(nest.dat.age$last.on, list(nest.dat.age$tpi50_3cat, nest.dat.age$F.Age), mat = TRUE)


first.off.age.tpi<-lme(first.off~ tpi50_3cat + F.Age + J.day,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(first.off.age.tpi)
emmeans(first.off.age.tpi, ~tpi50_3cat, method = "tukey")
pairs(emmeans(first.off.age.tpi, ~tpi50_3cat, method = "tukey"))
first.off.age<-lme(first.off~ F.Age + J.day,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(first.off.age)


last.on.age.tpi<-lme(last.on~ tpi50_3cat + F.Age + J.day,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(last.on.age.tpi)
emmeans(last.on.age.tpi, ~tpi50_3cat, method = "tukey")
pairs(emmeans(last.on.age.tpi, ~tpi50_3cat, method = "tukey"))

last.on.age<-lme(last.on~  F.Age + J.day,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(last.on.age)

last.on.tpi<-lme(last.on~ tpi50_3cat + J.day,  random = ~1|Nest.ID, method= "REML", data = nest.dat.age, na.action = na.omit)
summary(last.on.tpi)

ggplot(nest.dat.age, aes(J.day, last.on))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~tpi50_3cat)

ggplot(nest.dat.age, aes(J.day, first.off), color = F.Age)+
  geom_point()+
  geom_smooth()+
  facet_wrap(~tpi50_3cat)


#### figure grids ####
#### function for mylegend ####
library(gridExtra)

#### no.on.bouts.amb.T.fig####
no.on.bouts.amb.T.fig<-grid.arrange(arrangeGrob(on.bouts.Tmin.fig +theme(legend.position = "none"),
                                            on.bouts.Tmax.fig+theme(legend.position = "none"), 
                                            on.bouts.Tmean.fig+theme(legend.position = "none"), 
                                            on.bouts.CV.fig+theme(legend.position = "none"),
                                            nrow=2),
                                             nrow = 2, heights = c(10,1))

ggsave(path = "Figs.OnBouts","no.on.bouts.amb.T.fig.jpg", no.on.bouts.amb.T.fig, units= "in", height = 8, width = 7, dpi = 600)

#### no.off.bouts.amb.T.fig####
no.off.bouts.amb.T.fig<-grid.arrange(arrangeGrob(off.bouts.Tmin.fig +theme(legend.position = "none"),
                                                off.bouts.Tmax.fig+theme(legend.position = "none"), 
                                                off.bouts.Tmean.fig+theme(legend.position = "none"), 
                                                off.bouts.CV.fig+theme(legend.position = "none"),
                                                nrow=2),
                                    nrow = 2, heights = c(10,1))

ggsave(path = "Figs.OnBouts","no.off.bouts.amb.T.fig.jpg", no.off.bouts.amb.T.fig, units= "in", height = 8, width = 7, dpi = 600)

#### duration on.bouts.amb.T.fig####
dur.on.bouts.amb.T.fig<-grid.arrange(arrangeGrob(Ton.bouts.Tmin.fig +theme(legend.position = "none"),
                                                 Ton.bouts.Tmax.fig+theme(legend.position = "none"), 
                                                 Ton.bouts.Tmean.fig+theme(legend.position = "none"), 
                                                 Ton.bouts.TCV.fig+theme(legend.position = "none"),
                                                 nrow=2),
                                     nrow = 2, heights = c(10,1))

ggsave(path = "Figs.OnBouts","dur.on.bouts.amb.T.fig.jpg", dur.on.bouts.amb.T.fig, units= "in", height = 8, width = 7, dpi = 600)

#### duration off.bouts.amb.T.fig####
dur.off.bouts.amb.T.fig<-grid.arrange(arrangeGrob(Toff.bouts.Tmin.fig +theme(legend.position = "none"),
                                                 Toff.bouts.Tmax.fig+theme(legend.position = "none"), 
                                                 Toff.bouts.Tmean.fig+theme(legend.position = "none"), 
                                                 Toff.bouts.TCV.fig+theme(legend.position = "none"),
                                                 nrow=2),
                                     nrow = 2, heights = c(10,1))

ggsave(path = "Figs.OnBouts","dur.off.bouts.amb.T.fig.jpg", dur.off.bouts.amb.T.fig, units= "in", height = 8, width = 7, dpi = 600)

#### duration on.bouts.amb.T.fig with reduced data ####
dur.on.bouts.amb.T.red.fig<-grid.arrange(arrangeGrob(Ton.bouts.Tmin.red.fig +theme(legend.position = "none"),
                                                  Ton.bouts.Tmax.red.fig+theme(legend.position = "none"), 
                                                  Ton.bouts.Tmean.red.fig+theme(legend.position = "none"), 
                                                  Ton.bouts.TCV.red.fig+theme(legend.position = "none"),
                                                  nrow=2),
                                      nrow = 2, heights = c(10,1))

ggsave(path = "Figs.OnBouts","dur.on.bouts.amb.T.red.fig.jpg", dur.on.bouts.amb.T.red.fig, units= "in", height = 8, width = 7, dpi = 600)

#### duration off.bouts.amb.T.fig with reduced data####
dur.off.bouts.amb.T.red.fig<-grid.arrange(arrangeGrob(Toff.bouts.Tmin.red.fig +theme(legend.position = "none"),
                                                  Toff.bouts.Tmax.red.fig+theme(legend.position = "none"), 
                                                  Toff.bouts.Tmean.red.fig+theme(legend.position = "none"), 
                                                  Toff.bouts.TCV.red.fig+theme(legend.position = "none"),
                                                  nrow=2),
                                      nrow = 2, heights = c(10,1))

ggsave(path = "Figs.OnBouts","dur.off.bouts.amb.T.red.fig.jpg", dur.off.bouts.amb.T.red.fig, units= "in", height = 8, width = 7, dpi = 600)

#### proportion time on fig####
#note Tmax and Tmean no sig diff so not in figure
prop.on.amb.T.fig<-grid.arrange(arrangeGrob(Tmin.prop.fig +theme(legend.position = "none"),
                                            TCV.prop.fig+theme(legend.position = "none"), 
                                                  nrow=2),
                                                  nrow = 2, heights = c(10,1))
ggsave(path = "Figs.OnBouts","prop.on.amb.T.fig.jpg", prop.on.amb.T.fig, units= "in", height = 8, width = 4, dpi = 600)

prop.on.amb.T.row.fig<-grid.arrange(arrangeGrob(Tmin.prop.fig +theme(legend.position = "none"),
                                            TCV.prop.fig+theme(legend.position = "none"), 
                                            ncol=2))
                               
ggsave(path = "Figs.OnBouts","prop.on.amb.T.row.fig.jpg", prop.on.amb.T.row.fig, units= "in", height = 4, width = 7, dpi = 600)

#### Min nest temperatures and inc #### # nb redid with daily off bouts
Min.nest.T.figs<-grid.arrange(arrangeGrob(Tmin.n.off.fig+theme(legend.position = "none"),
                                          Tmin.dur.on.fig+theme(legend.position = "none"), 
                                          Tmin.T.off.fig+theme(legend.position = "none"), 
                                          P.on.T.min.fig+theme(legend.position = "none"),
                                                 nrow=2),
                                     nrow = 2, heights = c(10,1))

ggsave(path = "Figs.OnBouts","Min.nest.T.figs.jpg", Min.nest.T.figs, units= "in", height = 8, width = 7, dpi = 600)

#### Mean nest temps in inc #### # replaced with off bouts
Mean.nest.T.figs<-grid.arrange(arrangeGrob(Tmean.n.off.fig +theme(legend.position = "none"),
                                           Tmean.T.on.fig+theme(legend.position = "none"), 
                                           nest.MeanT.T.off.fig+theme(legend.position = "none"), 
                                          P.on.T.mean.fig+theme(legend.position = "none"),
                                          nrow=2),
                              nrow = 2, heights = c(10,1))

ggsave(path = "Figs.OnBouts","Mean.nest.T.figs.jpg", Mean.nest.T.figs, units= "in", height = 8, width = 7, dpi = 600)


#### CV nest temps in inc #### # replace on with off
CV.nest.T.figs<-grid.arrange(arrangeGrob(CV.n.off.fig +theme(legend.position = "none"),
                                           CV.T.on.fig+theme(legend.position = "none"), 
                                           CV.T.off.fig+theme(legend.position = "none"), 
                                           P.on.CV.fig+theme(legend.position = "none"),
                                           nrow=2),
                               nrow = 2, heights = c(10,1))

ggsave(path = "Figs.OnBouts","CV.nest.T.figs.jpg", CV.nest.T.figs, units= "in", height = 8, width = 7, dpi = 600)
