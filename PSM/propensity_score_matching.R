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


#do greedy matching on logit(PS) using Match with a caliper

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








#------------------GBM Method for Multi-level PSM-----------------------------------

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

## Estimate causal treatment effects
design.mnps <- svydesign(ids=~1, weights=~w,data = mydata)
glm1 <- svyglm(as.numeric(independent) ~ treatment, design = design.mnps, family=quasipoisson(link = log))
summary(glm1)





#------------------mean estimater-----------------------------------

treatment1<-as.numeric(data$treatment==1)
treatment2<-as.numeric(data$treatment==2)
treatment3<-as.numeric(data$treatment==3)
treatment<-as.character(data$treatment)
mydata1<-cbind(male,College,bachelor,master,
              docter,independent,careercount,
              isplatform,Officialcertification,treatment1,treatment2,treatment3,treatment)
mydata1<-data.frame(mydata1)

mnps.mydata <- mnps(treatment ~ male+College+bachelor+master+docter+careercount+isplatform+Officialcertification,
                    data = mydata1,
                    estimand = "ATE",
                    verbose = FALSE,
                    stop.method = c("es.mean", "ks.mean"),
                    n.trees = 3000)

mydata1$w <- get.weights(mnps.mydata, stop.method = "es.mean")
design.mnps <- svydesign(ids=~1, weights=~w,data = mydata1)
glm2 <- svyglm(as.numeric(independent) ~ treatment1+treatment2+male+College+bachelor+master+docter+careercount+isplatform+Officialcertification, design = design.mnps, family=quasipoisson(link = log))
summary(glm2)

#--------------------------------------------------------
#the estimater μ1-μ3 of coefficient β1,μ2-μ3 of coefficient β2,μ1-μ2 of coefficient β1-β2
mydata_T1_T2 <- data.frame(svycontrast(glm2, quote(treatment1 - treatment2)))
names(mydata_T1_T2) <- c("Estimate", "Std. Error")
mydata_T1_T2$"t value" <- mydata_T1_T2$Estimate/mydata_T1_T2$"Std. Error"
mydata_T1_T2$"Pr(>|t|)" <- 2*(1-pt(abs(mydata_T1_T2$"t value"),
                                  df=summary(glm2)$df.residual))
mydata_all <- rbind(summary(glm2)$coef[2:3,], mydata_T1_T2)
rownames(mydata_all) <- c("T1 vs. T3",
                        "T2 vs. T3", "T1 vs. T2")
mydata_all

#--------------------------------------------------------
## Combine all three mean
## each treatment mean μ1，μ2，μ3
#The intercept estimates the mean for T3
com_mean <- summary(glm2)$coef[1,]

pop_means <- data.frame(svycontrast(glm2, list(T1=c(1,1,0,0,0,0,0,0,0,0,0),
                                               T2=c(1,0,1,0,0,0,0,0,0,0,0))))
names(pop_means) <- c("Estimate", "Std. Error")
pop_means$"t value" <- pop_means$Estimate/pop_means$"Std. Error"
pop_means$"Pr(>|t|)" <- 2*(1-pt(abs(pop_means$"t value"), df=summary(glm2)$df.residual))
pop_means <- rbind(T3=com_mean, pop_means)
pop_means

