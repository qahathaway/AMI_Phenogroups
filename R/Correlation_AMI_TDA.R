#initialize and load the libraries
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(caTools)
library(ROCR)
library(stringr)
library(tidyverse)
library(linelist)
library(broom)
library(survival)
library(survminer)
library(riskRegression)
library(survcomp)


#reading the csv file with the data
df<- read.csv("/Users/quincy/Documents/Research/Rutgers/TDA/paper/Prism_Stats.csv" ,header=TRUE)
#which(is.na(df))

# remove certain unwanted columns
drop1 <- c(colnames(df)[0:1])
df = df[,!(names(df) %in% drop1)]
dff1 <- df[colMeans(is.na(df)) <= 0.1 & colMeans((df == 0), na.rm = T) <= 0.1]
dff1<-df
dff1[is.na(dff1)] <- 0

#Univariate Analysis
explanatory_vars <- c(colnames(dff1)[1:74])

explanatory_vars %>% str_c("Death_1yr ~ ", .)

models <- explanatory_vars %>%       # begin with variables of interest
  str_c("Death_1yr ~ ", .) %>%         # combine each variable into formula ("outcome ~ variable of interest")
  
  # iterate through each univariate formula
  map(                               
    .f = ~glm(                       # pass the formulas one-by-one to glm()
      formula = as.formula(.x),      # within glm(), the string formula is .x
      family = "binomial",           # specify type of glm (logistic)
      data = dff1)) %>%          # dataset
  
  # tidy up each of the glm regression outputs from above
  map(
    .f = ~tidy(
      .x, 
      exponentiate = TRUE,           # exponentiate 
      conf.int = TRUE)) %>%          # return confidence intervals
  
  # collapse the list of regression outputs in to one data frame
  bind_rows() %>% 
  
  # round all numeric columns
  mutate(across(where(is.numeric), round, digits = 4))

#Summary
models
summary(models)
write.csv(models, file="/Users/quincy/Documents/Research/Rutgers/TDA/paper/Univariate.csv")


#Multivariate Analysis
logistic_model <- glm(Death_1yr ~Groups+E.e..septal..112.,
                      data = dff1,
                      family=binomial)

# Summary
logistic_model
summary(logistic_model)

#Odds Ratio
require(MASS)
exp(cbind(coef(logistic_model), confint(logistic_model)))

#write.csv(models, file="/Users/quincy/Documents/Research/Rutgers/TDA/paper/Univariate.csv")



######HF Risk Prediction#####
##Training/Testing Split##
#set.seed(100)
#index <- createDataPartition(final$HF, p = 0.6, list = FALSE)
#train <- final[index, ]
#test <- final[-index, ]

##COXPH Univariate Analysis##
##Univariate##
covariates <- c(colnames(dff1)[5:89])
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(Death_1yr_Time, Death_1yr)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = dff1)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=4)
                         wald.test<-signif(x$wald["test"], digits=4)
                         beta<-signif(x$coef[1], digits=4);#coeficient beta
                         HR <-signif(x$coef[2], digits=4);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 4)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],4)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res1 <- as.data.frame(res)
print(res1[order(res1$p.value, decreasing = FALSE), ]   )

write.csv(res1[order(res1$p.value, decreasing = FALSE), ], file = "/Users/quincy/Documents/Research/Rutgers/TDA/paper/Cox_Univariate.csv")


# detach mlr package because of 'plotCalibration' conflict
detach("package:mlr", unload = TRUE)

##Survial Models##
coxDem <- coxph(Surv(Death_1yr_Time, Death_1yr)~Grace.score.at.6.months,data=dff1,x=TRUE,y=TRUE)
coxDemFunc <- coxph(Surv(Death_1yr_Time, Death_1yr)~Grace.score.at.6.months+LV.GLS..166.,data=dff1,x=TRUE,y=TRUE)
coxDemFuncOmics <- coxph(Surv(Death_1yr_Time, Death_1yr)~Grace.score.at.6.months+LV.GLS..166.+Group.C.Prob_NY,data=dff1,x=TRUE,y=TRUE)
coxDemOmics <- coxph(Surv(Death_1yr_Time, Death_1yr)~Grace.score.at.6.months+Group.C.Prob_NY,data=dff1,x=TRUE,y=TRUE)


AUC=Score(list(coxDem,coxDemFunc, coxDemFuncOmics), formula=Surv(Death_1yr_Time, Death_1yr)~1,
                  data=dff1, metrics="auc", plots=c("calibration","ROC"), times=seq(0:365))

aucgraphDem <- plotAUC(AUC, models = "coxph", which = "score", xlim = c(0,365), ylim = c(0,1), col = "bisque", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
aucgraphFunc <- plotAUC(AUC, models = "coxph.1", which = "score", xlim = c(0,365), ylim = c(0,1), col = "burlywood", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
aucgraphBoth <- plotAUC(AUC, models = "coxph.2", which = "score", xlim = c(0,365), ylim = c(0,1), col = "coral", lwd = 5, conf.int =TRUE, legend = FALSE)

Brier=Score(list(coxDem,coxDemFunc,coxDemFuncOmics), formula=Surv(Death_1yr_Time, Death_1yr)~1,
                    data=dff1, metrics="brier", plots=c("calibration","ROC"), times=seq(0:365))

briergraphDem <- plotBrier(Brier, models = "coxph", which = "score", xlim = c(0,365), ylim = c(0,1), col = "bisque", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
briergraphFunc <- plotBrier(Brier, models = "coxph.1", which = "score", xlim = c(0,365), ylim = c(0,1), col = "burlywood", lwd = 5, conf.int =TRUE, legend = FALSE)
par(new=TRUE)
briergraphBoth <- plotBrier(Brier, models = "coxph.2", which = "score", xlim = c(0,365), ylim = c(0,1), col = "coral", lwd = 5, conf.int =TRUE, legend = FALSE)


######Predict Risk#####
risk_score_Dem <- predictRisk(coxDem,times=c(365),newdata=dff1)
risk_score_DemFunc <- predictRisk(coxDemFunc,times=c(365),newdata=dff1)
risk_score_DemFuncOmics <- predictRisk(coxDemFuncOmics,times=c(365),newdata=dff1)
risk_score_DemOmics <- predictRisk(coxDemOmics,times=c(365),newdata=dff1)


######Concordance#####

ConDem <- concordance.index(x=risk_score_Dem, surv.time = dff1$Death_1yr_Time, surv.event = dff1$Death_1yr)
ConDemFunc <- concordance.index(x=risk_score_DemFunc, surv.time = dff1$Death_1yr_Time, surv.event = dff1$Death_1yr)
ConDemFuncOmics <- concordance.index(x=risk_score_DemFuncOmics, surv.time = dff1$Death_1yr_Time, surv.event = dff1$Death_1yr)
ConDemOmics <- concordance.index(x=risk_score_DemOmics, surv.time = dff1$Death_1yr_Time, surv.event = dff1$Death_1yr)

ApparrentCindex <- pec::cindex(list("COXPH Dem"=coxDem,
                                      "COXPH Function"=coxDemFunc,
                                      "COXPH Omics"=coxDemOmics,
                                      "COXPH Both"=coxDemFuncOmics),
                                 formula=Surv(Death_1yr_Time, Death_1yr)~1,data=dff1,
                                 eval.times=seq(0,365,1), pred.times=seq(0,365,1))

col = c("bisque", "coral", "burlywood", "chocolate")

plot(ApparrentCindex, legend = FALSE, xlim=c(0,400), ylim=c(0.5,1.0), lwd = 5, col = col)
write.csv(ApparrentCindex$AppCindex, file = "/Users/quincy/Documents/Research/Rutgers/TDA/paper/Concordance.csv")


#####Hazard Ratio#####
##Plot the baseline survival function##
fit <- surv_fit(Surv(Death_1yr_Time, Death_1yr) ~Group.A+Group.B+Group.C,
                data = dff1)

#ggsurvplot(fit, data = final, pval = TRUE, break.time.by = 500)

ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = FALSE,         # show confidence intervals for 
#  fun="event",
#  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in days",
  ylab = "Survival",# customize X axis label.
  fontsize = 2,
  break.time.by = 73,     # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = TRUE,
  risk.table.font = 3,
  risk.table.title = "Number at Risk (%)",
  risk.table.pos = "out",
  cumevents.title = "Cumulative Events",
  cumevents = FALSE,
  cumevents.font = 3,
  font.x = 12,
  font.y = 12,
  ncensor.plot = FALSE,      # plot the number of censored subjects at time t
  #  surv.median.line = "hv",  # add the median survival pointer.
  ylim = c(0.7, 1),
  xlim = c(0, 365),
  legend.labs = 
    c("Cluster C", "Cluster B", "Cluster A"),    # change legend labels.
  palette = 
    c("springgreen4", "darkorange2", "dodgerblue3") # custom color palettes.
)
