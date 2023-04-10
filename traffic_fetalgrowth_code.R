###***************************************************************************
###
### R code for analysis of 16-23 week biparietal diameter ultrasounds in:
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
#' @param ndlag NO2 exposure at a given gestational week (e.g., "ndlag01" is the NO2 value in the first gestational week)
#' @param templag temperature exposure at a given gestational week
#' @param pmlag PM2.5 exposure at a given gestational week


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

# 0c. divide ultrasound datasets based on time of measurement
df.bpd1 <- df.bpd %>% filter(ga_cat == 1) # 16-23 Week Scan
df.bpd2 <- df.bpd %>% filter(ga_cat == 2) # 24-31 Week Scan
df.bpd3 <- df.bpd %>% filter(ga_cat == 3) # 32+ Week Scan

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

# 1a.i. join with weekly exposure data
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

# 1a.v. generate DLM predictions
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
