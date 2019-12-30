library(Matching)
library(MASS)
library(readxl)
data <- read_xlsx("C:\\Users\\test\\Desktop\\test.xlsx")
names(data)
View(data)

male<-as.numeric(data$gender=='男')
College<-as.numeric(data$degree=='专科')
bachelor<-as.numeric(data$degree=='本科')
master<-as.numeric(data$degree=='硕士')
docter<-as.numeric(data$degree=='博士')
independent<-as.numeric(data$liveCount)
careercount<-as.numeric(data$careercount)
isplatform<-as.numeric(data$isplatform)
Officialcertification<-as.numeric(data$Officialcertification)
treatment<-as.numeric(data$treatment=='2')


mydata<-cbind(male,College,bachelor,master,docter,independent,careercount,isplatform,Officialcertification,treatment)
mydata<-data.frame(mydata)
View(mydata)
psmodel<-glm(treatment~male+College+bachelor+master+docter+careercount+isplatform+Officialcertification,family=binomial(),data=mydata)
summary(psmodel)
pscore<-psmodel$fitted.values


#-------------------------------do greedy matching on logit(PS) using Match with a caliper-----------------------------

logit <- function(p) {log(p)-log(1-p)}
psmatch<-Match(Tr=mydata$treatment,M=1,X=logit(pscore),replace=FALSE,caliper=.2)
matched<-mydata[unlist(psmatch[c("index.treated","index.control")]), ]
xvars<-c("male","College","bachelor","master","docter","careercount","isplatform","Officialcertification")

#get standardized differences
matchedtab1<-CreateTableOne(vars=xvars, strata ="treatment", 
                            data=matched, test = FALSE)
print(matchedtab1, smd = TRUE)

#outcome analysis
y_trt<-matched$independent[matched$treatment==1]
y_con<-matched$independent[matched$treatment==0]

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)








#------------------------------------mnps Method for Multi-level treatment PSM-----------------------------------

library(twang)
library(survey)

male<-as.numeric(data$gender=='男')
College<-as.numeric(data$degree=='专科')
bachelor<-as.numeric(data$degree=='本科')
master<-as.numeric(data$degree=='硕士')
docter<-as.numeric(data$degree=='博士')
independent<-as.numeric(data$liveCount)
careercount<-as.numeric(data$careercount)
isplatform<-as.numeric(data$isplatform)
Officialcertification<-as.numeric(data$Officialcertification)
treatment<-as.character(data$treatment)
mydata<-cbind(male,College,bachelor,master,docter,independent,careercount,isplatform,Officialcertification,treatment)
mydata<-data.frame(mydata)
View(mydata)

##Estimate propensity scores for each treatment by running mnps() function in twang
mnps.mydata <- mnps(treatment ~ male+College+bachelor+master+docter+careercount+isplatform+Officialcertification,
               data = mydata,
               estimand = "ATE",
               verbose = FALSE,
               stop.method = c("es.mean", "ks.mean"),
               n.trees = 3000)

##Assess balance using desired balance metric(s)
plot(mnps.mydata, plots = 1)[[1]]
plot(mnps.mydata, plots = 1)[[2]]
plot(mnps.mydata, plots = 1)[[3]]
##Assess overlap among treatment samples
plot(mnps.mydata, plots = 2, subset = "es.mean")[[1]]
plot(mnps.mydata, plots = 2, subset = "es.mean")[[2]]
plot(mnps.mydata, plots = 2, subset = "es.mean")[[3]]
plot(mnps.mydata, plots = 3, pairwiseMax = FALSE, figureRows = 3)
plot(mnps.mydata, plots = 4)

bal.table(mnps.mydata, digits = 2)
summary(mnps.mydata)


# Estimate the weighted mean for each treatment
mydata$w <- get.weights(mnps.mydata, stop.method = "es.mean")



#--------------------------------------------------------
## Estimate causal treatment effects
design.mnps <- svydesign(ids=~1, weights=~w,data = mydata)
glm1 <- svyglm(as.numeric(independent) ~ treatment, design = design.mnps)
summary(glm1)

#the mean difference estimater μ2-μ1 of coefficient β1,μ3-μ1 of coefficient β2,μ2-μ3 of coefficient β1-β2
mydata_T1_T2 <- data.frame(svycontrast(glm1, quote(treatment2 - treatment3)))
names(mydata_T1_T2) <- c("Estimate", "Std. Error")
mydata_T1_T2$"t value" <- mydata_T1_T2$Estimate/mydata_T1_T2$"Std. Error"
mydata_T1_T2$"Pr(>|t|)" <- 2*(1-pt(abs(mydata_T1_T2$"t value"),
                                  df=summary(glm1)$df.residual))
mydata_all <- rbind(summary(glm1)$coef[2:3,], mydata_T1_T2)
rownames(mydata_all) <- c("T2 vs. T1", "T3 vs. T1", "T2 vs. T3")
mydata_all

#--------------------------------------------------------
## Combine all three mean
## each treatment mean μ1，μ2，μ3
#The intercept estimates the mean for T3
glm2 <- svyglm(as.numeric(independent) ~ treatment, 
               design = design.mnps,contrast=list(treatment=contr.sum))
summary(glm2)
T3_mean$"Estimate" <- -sum(coef(glm2)[-1])
T3_mean$"Std. Error" <- sqrt(c(-1,-1) %*% summary(glm2)$cov.scaled[-1,-1] %*% c(-1,-1))
T3_mean$"t value" <- T3_mean$"Estimate"/T3_mean$"Std. Error"
T3_mean$"Pr(>|t|)" <- 2*(1-pt(abs(T3_mean$"t value"),
                                   df=summary(glm2)$df.residual))
T1_mean <- summary(glm2)$coef[2,]
T2_mean <- summary(glm2)$coef[3,]
T3_mean <- cbind(T3_mean$"Estimate", T3_mean$"Std. Error", T3_mean$"t value", T3_mean$"Pr(>|t|)")
means_all <- rbind(T1_mean, T2_mean, T3_mean)
rownames(means_all) <- c("T1", "T2", "T3")
means_all




###independent value is count number
glm1 <- svyglm(as.numeric(independent) ~ treatment, design = design.mnps, family=quasipoisson(link = 'log'))
summary(glm1)
glm2 <- svyglm(as.numeric(independent) ~ treatment, design = design.mnps, family=quasipoisson(link = 'log'), contrast=list(treatment=contr.sum))
summary(glm2)






#---------------------------------------------IPTW method for PSM----------------------------------
#RHC Example

#install packages (if needed)
install.packages("tableone")
install.packages("ipw")
install.packages("sandwich")
install.packages("survey")

#load packages
library(tableone)
library(ipw)
library(sandwich) #for robust variance estimation
library(survey)

expit <- function(x) {1/(1+exp(-x)) }
logit <- function(p) {log(p)-log(1-p)}

#read in data
load(url("http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/rhc.sav"))
#view data
View(rhc)

#treatment variables is swang1
#x variables that we will use
#cat1: primary disease category
#age
#sex
#meanbp1: mean blood pressure

#create a data set with just these variables, for simplicity
ARF<-as.numeric(rhc$cat1=='ARF')
CHF<-as.numeric(rhc$cat1=='CHF')
Cirr<-as.numeric(rhc$cat1=='Cirrhosis')
colcan<-as.numeric(rhc$cat1=='Colon Cancer')
Coma<-as.numeric(rhc$cat1=='Coma')
COPD<-as.numeric(rhc$cat1=='COPD')
lungcan<-as.numeric(rhc$cat1=='Lung Cancer')
MOSF<-as.numeric(rhc$cat1=='MOSF w/Malignancy')
sepsis<-as.numeric(rhc$cat1=='MOSF w/Sepsis')
female<-as.numeric(rhc$sex=='Female')
died<-as.integer(rhc$death=='Yes')
age<-rhc$age
treatment<-as.numeric(rhc$swang1=='RHC')
meanbp1<-rhc$meanbp1

#new dataset
mydata<-cbind(ARF,CHF,Cirr,colcan,Coma,lungcan,MOSF,sepsis,
              age,female,meanbp1,treatment,died)
mydata<-data.frame(mydata)

#covariates we will use (shorter list than you would use in practice)
xvars<-c("age","female","meanbp1","ARF","CHF","Cirr","colcan",
         "Coma","lungcan","MOSF","sepsis")

#look at a table 1
table1<- CreateTableOne(vars=xvars,strata="treatment", data=mydata, test=FALSE)
## include standardized mean difference (SMD)
print(table1,smd=TRUE)

#propensity score model
psmodel <- glm(treatment ~ age + female + meanbp1+ARF+CHF+Cirr+colcan+
                 Coma+lungcan+MOSF+sepsis,
               family  = binomial(link ="logit"))

## value of propensity score for each subject
ps <-predict(psmodel, type = "response")

#create weights
weight<-ifelse(treatment==1,1/(ps),1/(1-ps))

#apply weights to data
weighteddata<-svydesign(ids = ~ 1, data =mydata, weights = ~ weight)

#weighted table 1
weightedtable <-svyCreateTableOne(vars = xvars, strata = "treatment", 
                                  data = weighteddata, test = FALSE)
## Show table with SMD
print(weightedtable, smd = TRUE)

#to get a weighted mean for a single covariate directly:
mean(weight[treatment==1]*age[treatment==1])/(mean(weight[treatment==1]))

#get causal risk difference
glm.obj<-glm(died~treatment,weights=weight,family=quasibinomial(link="identity"))
#summary(glm.obj)
betaiptw<-coef(glm.obj)
SE<-sqrt(diag(vcovHC(glm.obj, type="HC0")))

causalrd<-(betaiptw[2])
lcl<-(betaiptw[2]-1.96*SE[2])
ucl<-(betaiptw[2]+1.96*SE[2])
c(lcl,causalrd,ucl)

#get causal relative risk. Weighted GLM
glm.obj<-glm(died~treatment,weights=weight,family=quasibinomial(link=log))
#summary(glm.obj)
betaiptw<-coef(glm.obj)
#to properly account for weighting, use asymptotic (sandwich) variance
SE<-sqrt(diag(vcovHC(glm.obj, type="HC0")))

#get point estimate and CI for relative risk (need to exponentiate)
causalrr<-exp(betaiptw[2])
lcl<-exp(betaiptw[2]-1.96*SE[2])
ucl<-exp(betaiptw[2]+1.96*SE[2])
c(lcl,causalrr,ucl)

#truncate weights at 10

truncweight<-replace(weight,weight>10,10)
#get causal risk difference
glm.obj<-glm(died~treatment,weights=truncweight,family=quasibinomial(link="identity"))
#summary(glm.obj)
betaiptw<-coef(glm.obj)
SE<-sqrt(diag(vcovHC(glm.obj, type="HC0")))

causalrd<-(betaiptw[2])
lcl<-(betaiptw[2]-1.96*SE[2])
ucl<-(betaiptw[2]+1.96*SE[2])
c(lcl,causalrd,ucl)

#############################
#alternative: use ipw package
#############################

#first fit propensity score model to get weights
weightmodel<-ipwpoint(exposure= treatment, family = "binomial", link ="logit",
                      denominator= ~ age + female + meanbp1+ARF+CHF+Cirr+colcan+
                        Coma+lungcan+MOSF+sepsis, data=mydata)
#numeric summary of weights
summary(weightmodel$ipw.weights)
#plot of weights
ipwplot(weights = weightmodel$ipw.weights, logscale = FALSE,
        main = "weights", xlim = c(0, 22))
mydata$wt<-weightmodel$ipw.weights

#fit a marginal structural model (risk difference)
msm <- (svyglm(died ~ treatment, design = svydesign(~ 1, weights = ~wt,
                                                    data =mydata)))
coef(msm)
confint(msm)


# fit propensity score model to get weights, but truncated
weightmodel<-ipwpoint(exposure= treatment, family = "binomial", link ="logit",
                      denominator= ~ age + female + meanbp1+ARF+CHF+Cirr+colcan+
                        Coma+lungcan+MOSF+sepsis, data=mydata,trunc=.01)

#numeric summary of weights
summary(weightmodel$weights.trun)
#plot of weights
ipwplot(weights = weightmodel$weights.trun, logscale = FALSE,
        main = "weights", xlim = c(0, 22))
mydata$wt<-weightmodel$weights.trun
#fit a marginal structural model (risk difference)
msm <- (svyglm(died ~ treatment, design = svydesign(~ 1, weights = ~wt,
                                                    data =mydata)))
coef(msm)
confint(msm)
