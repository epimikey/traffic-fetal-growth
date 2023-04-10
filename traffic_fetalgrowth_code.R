###***************************************************************************
###
### R code for analysis of 16-23 week ultrasounds in:
###
### "Leung M, Modest AM, Hacker MR, Wylie BJ, Wei Y, Schwartz J, Iyer HS,
###  Hart JE, Coull BA, Laden F, Weisskopf MG, Papatheodorou S. 
###  Traffic-related air pollution and ultrasound parameters of fetal growth
###  in Eastern Massachusetts, USA. American Journal of Epidemiology. 2023."
###
###***************************************************************************

###*************
### N: Notes ###
###*************

#' @param ga_cat gestational age category (1: 16-23 week scan, 2: 24-31 week scan, 3: 32+ week scan)
#' @param week gestational week
#' @param zscore ultrasound parameter z-score based on INTERGROWTH-21st standards
#' @param doc.numeric day of conception within the study period (to control for seasonality & long-term trends)
#' @param mage maternal age at conception
#' @param meduc maternal educational attainment
#' @param parity_first parity
#' @param sex fetal sex
#' @param ins insurance type
#' @param adi_natrank Area Deprivation Index national percentiles
#' @param smoking maternal smoking during pregnancy
#' @param ndlag NO2 lag at a given lag (e.g., "ndlag01" is the NO2 1 week prior to the ultrasound measurement)
#' @param templag temperature lag at a given lag 
#' @param pmlag PM2.5 lag at a given lag 


###****************************
### 0: Read in data & setup ###
###****************************

# 0a. install required packages
requiredPackages = c('tidyverse','gamm4','dlnm','splines')
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE))
    install.packages(p)
  library(p, character.only = TRUE)
}

# 0b. read in analytic datasets
df.bpd <- readRDS("DIRECTORY/df_bpd.rds") # biparietal diameter ultrasounds
df.hc <- readRDS("DIRECTORY/df_hc.rds") # head circumference ultrasounds
df.fl <- readRDS("DIRECTORY/df_fl.rds") # femur length ultrasounds
df.ac <- readRDS("DIRECTORY/df_ac.rds") # abdominal circumference ultrasounds

# 0c. set up ultrasound datasets
df.bpd1 <- df.bpd %>% filter(ga_cat == 1) # 16-23 Week Scan
df.bpd2 <- df.bpd %>% filter(ga_cat == 2) # 24-31 Week Scan
df.bpd3 <- df.bpd %>% filter(ga_cat == 3) # 32+ Week Scan

df.hc1 <- df.hc %>% filter(ga_cat == 1)
df.hc2 <- df.hc %>% filter(ga_cat == 2)
df.hc3 <- df.hc %>% filter(ga_cat == 3)

df.fl1 <- df.fl %>% filter(ga_cat == 1)
df.fl2 <- df.fl %>% filter(ga_cat == 2)
df.fl3 <- df.fl %>% filter(ga_cat == 3)

df.ac1 <- df.ac %>% filter(ga_cat == 1)
df.ac2 <- df.ac %>% filter(ga_cat == 2)
df.ac3 <- df.ac %>% filter(ga_cat == 3)

# 0d. define function to test AIC
testAIC <- function(nd_hist, temp_hist, n_df1, data){
  cb.nd <<- crossbasis(nd_hist, argvar = list("lin"), arglag = list(df=n_df1))
  cb.temp <<- crossbasis(temp_hist, argvar = list(df=3), arglag = list(df=4))
  
  model <- glm(zscore ~ cb.nd + cb.temp + 
                 ns(doc.numeric,35) +
                 poly(mage,2) + as.factor(meduc) + as.factor(mrace) + as.factor(parity_first) +
                 as.factor(sex) + as.factor(ins) + poly(adi_natrank,2) + as.factor(smoking), 
               family = gaussian(), 
               data = data)
  AIC(model)
}

###**********************************
### 1. Ultrasound 16-23 Week Scan ###
###**********************************

### 1a. Biparietal Diameter

# 1a.i. Join with weekly exposure data
df.bpd1 <- df.bpd1 %>% 
  left_join(nd.week16, by="baby_mrn") %>% 
  left_join(temp.week16, by="baby_mrn") %>%
  left_join(pm.week16, by="baby_mrn") %>%
  na.omit()

nd.bpd1 <- df.bpd1 %>% select(ndlag01:ndlag16) # NO2 exposure history from week 1-16  for each person (wide format)
temp.bpd1 <- df.bpd1 %>% select(templag01:templag16) # temperature exposure history from week 1-16 for each person (wide format)
pm.bpd1 <- df.bpd1 %>% select(pmlag01:pmlag16) # PM exposure history from week 1-16 for each person (wide format)

# 1a.ii. test AIC
bpd1.aic <- data.frame(aic = numeric())
for(i in 2:5){bpd1.aic <- rbind(bpd1.aic, testAIC(nd.bpd1, temp.bpd1, i, df.bpd1))} # loop through 2-5 dfs for the lag constraint 
bpd1.aic <- bpd1.aic %>% mutate(n_df1 = 2:5) %>% rename(aic = 1) %>% filter(aic==min(aic)) %>% as.matrix() # select the one that minimizes the AIC

# 1a.iii. create cross basis
cb.nd.bpd1 <- crossbasis(nd.bpd1, argvar = list("lin"), arglag = list(df=bpd1.aic[2]))
cb.temp.bpd1 <- crossbasis(temp.bpd1, argvar = list(df=3), arglag = list(df=4))

# 1a.iv. models
model.bpd1 <- gamm4(zscore ~ cb.nd.bpd1 + cb.temp.bpd1 + 
                      ns(doc.numeric,35) + # control for seasonality & long-term trends
                      poly(mage,2) + as.factor(meduc) + as.factor(mrace) + as.factor(parity_first) +
                      as.factor(sex) + as.factor(ins) + poly(adi_natrank,2) + as.factor(smoking), 
                    family = gaussian(), 
                    random = ~(1|baby_mrn), # random intercept for each fetus
                    data = df.bpd1)$gam

# 1a.v. predictions
pred.bpd1 <- crosspred(cb.nd.bpd1, model.bpd1, at=30, cen=20) # for a 10-ppb contrast

# 1a.vi. extract week-specific estimate & 95% CIs
bpd1.nd <- cbind(data.frame(pred.bpd1$matfit[1,]), data.frame(pred.bpd1$matlow[1,]), data.frame(pred.bpd1$mathigh[1,])) %>%
  mutate(week = 1:16, outcome = "BPD", time="16-23 Week Scan") 
colnames(bpd1.nd) <- c("est","lci","uci","week","outcome","time") 

# 1a.vii. extract the cumulative estimate & 95% CI
bpd1.nd.cum <- cbind(data.frame(pred.bpd1$allfit[1]), data.frame(pred.bpd1$alllow[1]), data.frame(pred.bpd1$allhigh[1])) 
colnames(bpd1.nd.cum) <- c("est","lci","uci")
bpd1.nd.cum <- bpd1.nd.cum %>%
  mutate(est = paste(formatC(est,digits=2,format="f"), " (",formatC(lci,digits=2,format="f"),", ",formatC(uci,digits=2,format="f"),")",sep="")) %>%
  select(est) %>% mutate(outcome = "BPD", time = "16-23 Week Scan")

### 1b. Head Circumference

# 1b.i. join with weekly exposure data
df.hc1 <- df.hc1 %>% 
  left_join(nd.week16, by="baby_mrn") %>% 
  left_join(temp.week16, by="baby_mrn") %>%
  na.omit()

nd.hc1 <- df.hc1 %>% select(ndlag01:ndlag16)
temp.hc1 <- df.hc1 %>% select(templag01:templag16)

# 1b.ii. test AIC
hc1.aic <- data.frame(aic = numeric())
for(i in 2:5){hc1.aic <- rbind(hc1.aic, testAIC(nd.hc1, temp.hc1, i, df.hc1))}
hc1.aic <- hc1.aic %>% mutate(n_df1 = 2:5) %>% rename(aic = 1) %>% filter(aic==min(aic)) %>% as.matrix()

# 1b.iii. create cross basis
cb.nd.hc1 <- crossbasis(nd.hc1, argvar = list("lin"), arglag = list(df=hc1.aic[2]))
cb.temp.hc1 <- crossbasis(temp.hc1, argvar = list(df=3), arglag = list(df=4))

# 1b.iv. models
model.hc1 <- gamm4(zscore ~ cb.nd.hc1 + cb.temp.hc1 +  
                     ns(doc.numeric,35) +
                     poly(mage,2) + as.factor(meduc) + as.factor(mrace) + as.factor(parity_first) +
                     as.factor(sex) + as.factor(ins) + poly(adi_natrank,2) + as.factor(smoking), 
                   random = ~(1|baby_mrn),
                   data = df.hc1)$gam

# 1b.v. predictions
pred.hc1 <- crosspred(cb.nd.hc1, model.hc1, at=30, cen=20)

# 1b.vi. extract week-specific estimate & 95% CIs
hc1.nd <- cbind(data.frame(pred.hc1$matfit[1,]), data.frame(pred.hc1$matlow[1,]), data.frame(pred.hc1$mathigh[1,])) %>%
  mutate(week = 1:16, outcome = "HC", time="16-23 Week Scan")
colnames(hc1.nd) <- c("est","lci","uci","week","outcome","time")

# 1b.vii. extract the cumulative estimate & 95% CIs
hc1.nd.cum <- cbind(data.frame(pred.hc1$allfit[1]), data.frame(pred.hc1$alllow[1]), data.frame(pred.hc1$allhigh[1]))
colnames(hc1.nd.cum) <- c("est","lci","uci")
hc1.nd.cum <- hc1.nd.cum %>%
  mutate(est = paste(formatC(est,digits=2,format="f"), " (",formatC(lci,digits=2,format="f"),", ",formatC(uci,digits=2,format="f"),")",sep="")) %>%
  select(est) %>% mutate(outcome = "HC", time = "16-23 Week Scan")


### 1c. Femur Length

# 1c.i. join with weekly exposure data
df.fl1 <- df.fl1 %>% 
  left_join(nd.week16, by="baby_mrn") %>% 
  left_join(temp.week16, by="baby_mrn") %>%
  na.omit()

nd.fl1 <- df.fl1 %>% select(ndlag01:ndlag16)
temp.fl1 <- df.fl1 %>% select(templag01:templag16)

# 1c.ii. test AIC
fl1.aic <- data.frame(aic = numeric())
for(i in 2:5){fl1.aic <- rbind(fl1.aic, testAIC(nd.fl1, temp.fl1, i, df.fl1))}
fl1.aic <- fl1.aic %>% mutate(n_df1 = 2:5) %>% rename(aic = 1) %>% filter(aic==min(aic)) %>% as.matrix()

# 1c.iii. create cross basis
cb.nd.fl1 <- crossbasis(nd.fl1, argvar = list("lin"), arglag = list(df=fl1.aic[2]))
cb.temp.fl1 <- crossbasis(temp.fl1, argvar = list(df=3), arglag = list(df=4))

# 1c.iv. models
model.fl1 <- gamm4(zscore ~ cb.nd.fl1 + cb.temp.fl1 + 
                     ns(doc.numeric,35) +
                     poly(mage,2) + as.factor(meduc) + as.factor(mrace) + as.factor(parity_first) +
                     as.factor(sex) + as.factor(ins) + poly(adi_natrank,2) + as.factor(smoking), 
                   family = gaussian(), 
                   random = ~(1|baby_mrn),
                   data = df.fl1)$gam

# 1c.v. predictions
pred.fl1 <- crosspred(cb.nd.fl1, model.fl1, at=30, cen=20)

# 1c.vi. extract week-specific estimate & 95% CIs
fl1.nd <- cbind(data.frame(pred.fl1$matfit[1,]), data.frame(pred.fl1$matlow[1,]), data.frame(pred.fl1$mathigh[1,])) %>%
  mutate(week = 1:16, outcome = "FL", time="16-23 Week Scan")
colnames(fl1.nd) <- c("est","lci","uci","week","outcome","time")

# 1c.vii. extract the cumulative estimate & 95% CIs
fl1.nd.cum <- cbind(data.frame(pred.fl1$allfit[1]), data.frame(pred.fl1$alllow[1]), data.frame(pred.fl1$allhigh[1]))
colnames(fl1.nd.cum) <- c("est","lci","uci")
fl1.nd.cum <- fl1.nd.cum %>%
  mutate(est = paste(formatC(est,digits=2,format="f"), " (",formatC(lci,digits=2,format="f"),", ",formatC(uci,digits=2,format="f"),")",sep="")) %>%
  select(est) %>% mutate(outcome = "FL", time = "16-23 Week Scan")


### 1d. Abdominal Circumference

# 1d.i. join with weekly exposure data
df.ac1 <- df.ac1 %>% 
  left_join(nd.week16, by="baby_mrn") %>% 
  left_join(temp.week16, by="baby_mrn") %>%
  na.omit()

nd.ac1 <- df.ac1 %>% select(ndlag01:ndlag16)
temp.ac1 <- df.ac1 %>% select(templag01:templag16)

# 1d.ii. test AIC
ac1.aic <- data.frame(aic = numeric())
for(i in 2:5){ac1.aic <- rbind(ac1.aic, testAIC(nd.ac1, temp.ac1, i, df.ac1))}
ac1.aic <- ac1.aic %>% mutate(n_df1 = 2:5) %>% rename(aic = 1) %>% filter(aic==min(aic)) %>% as.matrix()

# 1d.iii. create cross basis
cb.nd.ac1 <- crossbasis(nd.ac1, argvar = list("lin"), arglag = list(df=ac1.aic[2]))
cb.temp.ac1 <- crossbasis(temp.ac1, argvar = list(df=3), arglag = list(df=4))

# 1d.iv. models
model.ac1 <- gamm4(zscore ~ cb.nd.ac1 + cb.temp.ac1 + 
                     ns(doc.numeric,35) +
                     poly(mage,2) + as.factor(meduc) + as.factor(mrace) + as.factor(parity_first) +
                     as.factor(sex) + as.factor(ins) + poly(adi_natrank,2) + as.factor(smoking), 
                   family = gaussian(), 
                   random = ~(1|baby_mrn),
                   data = df.ac1)$gam

# 1d.v. predictions
pred.ac1 <- crosspred(cb.nd.ac1, model.ac1, at=30, cen=20)

# 1d.vi. extract week-specific estimate & 95% CIs
ac1.nd <- cbind(data.frame(pred.ac1$matfit[1,]), data.frame(pred.ac1$matlow[1,]), data.frame(pred.ac1$mathigh[1,])) %>%
  mutate(week = 1:16, outcome = "AC", time="16-23 Week Scan")
colnames(ac1.nd) <- c("est","lci","uci","week","outcome","time")

# 1d.vii. extract the cumulative estimate & 95% CIs
ac1.nd.cum <- cbind(data.frame(pred.ac1$allfit[1]), data.frame(pred.ac1$alllow[1]), data.frame(pred.ac1$allhigh[1]))
colnames(ac1.nd.cum) <- c("est","lci","uci")
ac1.nd.cum <- ac1.nd.cum %>%
  mutate(est = paste(formatC(est,digits=2,format="f"), " (",formatC(lci,digits=2,format="f"),", ",formatC(uci,digits=2,format="f"),")",sep="")) %>%
  select(est) %>% mutate(outcome = "AC", time = "16-23 Week Scan")
