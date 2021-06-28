library(glmmTMB)
library(performance)
library(DHARMa)
library(car)
library(AICcmodavg)
library(MuMIn)
library(emmeans)
library(multcomp)
library(mgcv)
library(gamm4)
library(itsadug)
library(readxl)

source(system.file("other_methods","lsmeans_methods.R",package="glmmTMB"))

AGO <- read_excel("Supplementary Dataset 1.xlsx", 
                  sheet = "AGOcov", col_types = c("numeric", 
                                                  "date", "numeric", "text", "numeric", 
                                                  "text", "text", "text", "text", "text", 
                                                  "text", "text", "text", "text", "text", 
                                                  "text", "numeric", "numeric", "text", 
                                                  "text", "numeric", "numeric", "text", 
                                                  "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", 
                                                  "text", "numeric", "numeric", "numeric", 
                                                  "numeric", "numeric", "numeric", 
                                                  "numeric","numeric"))


AGO$IntPhase <- relevel(as.factor(AGO$IntPhase), ref = "Pre")

#Evaluating the best distribution for the dataset with count data
AGOPs <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                 + Year*IntPhase + PRCP + TAVG + (1|Community/HouseID) 
                 + (1|Week), family = "poisson", data = AGO)
AGONB1 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                  + Year*IntPhase + PRCP + TAVG + (1|Community/HouseID) 
                  + (1|Week), family = "nbinom1", data = AGO)
AGONB2 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                  + Year*IntPhase + PRCP + TAVG + (1|Community/HouseID) 
                  + (1|Week), family = "nbinom2", data = AGO)

#Checking the AIC values for each model
AIC(AGOPs)
AIC(AGONB1)
AIC(AGONB2)

#Evaluation of the overdispersion of the dataset 
check_overdispersion(AGOPs)

#Visualization of the data
AGOPs.V1 <- simulateResiduals(fittedModel = AGOPs, n = 250, plot = TRUE)
AGONB1.V1 <- simulateResiduals(fittedModel = AGONB1, n = 250, plot = TRUE)
AGONB2.V1 <- simulateResiduals(fittedModel = AGONB2, n = 250, plot = TRUE)

#explore interaction contrasts of best fit
a1 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
              + Year*IntPhase + PRCP + TAVG + (1|Community/HouseID) + (1|Week), 
              family = "nbinom2", data = AGO)
summary(a1)
a2 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
              + Year*IntPhase + TAVG + (1|Community/HouseID) + (1|Week), 
              family = "nbinom2", data = AGO)
summary(a2)
a3 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
              + Year*IntPhase + (1|Community/HouseID) + (1|Week), 
              family = "nbinom2", data = AGO)
summary(a3)
a4 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Placement + Income
              + Year*IntPhase + (1|Community/HouseID) + (1|Week), 
              family = "nbinom2", data = AGO)
summary(a4)
a5 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Placement + Income + 
                Year + IntPhase + (1|Community/HouseID) + (1|Week), 
              family = "nbinom2", data = AGO)
summary(a5)

anova(a1,a2,a3,a4,a5)

AICtab(a1, a2, a3, a4, a5, base = T, weights = T, sort = FALSE)

r.squaredGLMM(a3, pj2014 = FALSE)

#exploring the models a2 and a3
anova(a2, a3)

#This next code takes a lot of time, skip to the joint_test and plot effects if already tested
confint(a3, method = "uniroot")

car::Anova(a3, type = "III")

#Checking for multicollinearity
results <- check_collinearity(a3)
results
plot(results)

#Evaluating if the intervention had an effect on the abundance of female Ae. aegypti.
a3.emm.s <- emmeans(a3, ~ IntPhase*Year)
plot(allEffects(a3))

pairs(a3.emm.s)
multcomp::cld(a3.emm.s)

#following [mean(Cont)-mean(Int)]+\-1.96{sqrt[SE(Cont)^2+SE(Int)^2]} to calculate z-confidence interval at (alpha=0.05)
#for comparison between Control and Intervention for the intervention period of 2017 indoor abundance
(-0.453-0.474)-(1.96*sqrt((0.363^2)+(0.359^2))) #2017
(-0.453-0.474)+(1.96*sqrt((0.363^2)+(0.359^2))) 

write.csv(a3.emm.s, file = "a3emm.csv")


#GLMM for transiet effects, checking Short time frames

#Evaluating the best distribution for the dataset with count data
AGOweekPs <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                     + Year*IntPhase + PRCP + TAVG + Week + (1|Community/HouseID) 
                     , family = "poisson", data = AGO)
AGOweekNB1 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                      + Year*IntPhase + PRCP + TAVG + Week + (1|Community/HouseID) 
                      , family = "nbinom1", data = AGO)
AGOweekNB2 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                      + Year*IntPhase + PRCP + TAVG + Week + (1|Community/HouseID) 
                      , family = "nbinom2", data = AGO)

#Checking the AIC values for each model
AIC(AGOweekPs)
AIC(AGOweekNB1)
AIC(AGOweekNB2)

#Evaluation of the overdispersion of the dataset 
check_overdispersion(AGOweekPs)

#Visualization of the data
AGOweekPs.V1 <- simulateResiduals(fittedModel = AGOweekPs, n = 250, plot = TRUE)
AGOweekNB1.V1 <- simulateResiduals(fittedModel = AGOweekNB1, n = 250, plot = TRUE)
AGweekONB2.V1 <- simulateResiduals(fittedModel = AGOweekNB2, n = 250, plot = TRUE)

#Model reduction
a6 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
              + Year*IntPhase + PRCP + TAVG + Week + (1|Community/HouseID) 
              , family = "nbinom2", data = AGO)
summary(a6)
a7 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
              + Year*IntPhase + TAVG + Week + (1|Community/HouseID) 
              , family = "nbinom2", data = AGO)
summary(a7)
a8 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
              + Year*IntPhase + Week + (1|Community/HouseID) 
              , family = "nbinom2", data = AGO)
summary(a8)

AICtab(a6,a7,a8, base = T, weights = T, sort = FALSE)

r.squaredGLMM(a6, pj2014 = FALSE)

#This next code takes a lot of time, skip to the joint_test and plot effects if already tested
confint(a6, method = "uniroot")

car::Anova(a6, type = "III")

#Checking for multicollinearity
results <- check_collinearity(a6)
results
plot(results)

#Evaluating if the intervention on a weekly scale had an effect on the abundance of female Ae. aegypti.
a6.emm.s <- emmeans(a6, ~ IntPhase*Year)
plot(allEffects(a6))

pairs(a6.emm.s)
multcomp::cld(a6.emm.s)

#following [mean(Cont)-mean(Int)]+\-1.96{sqrt[SE(Cont)^2+SE(Int)^2]} to calculate z-confidence interval at (alpha=0.05)
#for comparison between Control and Intervention for the intervention period of 2017
(-0.972+0.039)-(1.96*sqrt((0.268^2)+(0.264^2))) #2017
(-0.972+0.039)+(1.96*sqrt((0.268^2)+(0.264^2)))

(-0.600+1.187)-(1.96*sqrt((0.263^2)+(0.267^2))) #2018
(-0.600+1.187)+(1.96*sqrt((0.263^2)+(0.267^2)))


write.csv(a6.emm.s, file = "a6emm.csv")


#GLMM for transiet effects, checking with shorter time frames Monthly

#Evaluating the best distribution for the dataset with count data
AGOmonthPs <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                      + Year*IntPhase + PRCP + TAVG + MonthN + (1|Community/HouseID) 
                      , family = "poisson", data = AGO)
AGOmonthNB1 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                       + Year*IntPhase + PRCP + TAVG + MonthN + (1|Community/HouseID) 
                       , family = "nbinom1", data = AGO)
AGOmonthNB2 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                       + Year*IntPhase + PRCP + TAVG + MonthN + (1|Community/HouseID) 
                       , family = "nbinom2", data = AGO)

#Checking the AIC values for each model
AIC(AGOmonthPs)
AIC(AGOmonthNB1)
AIC(AGOmonthNB2)

#Evaluation of the overdispersion of the dataset 
check_overdispersion(AGOmonthPs)

#Visualization of the data
AGOmonthPs.V1 <- simulateResiduals(fittedModel = AGOmonthPs, n = 250, plot = TRUE)
AGOmonthNB1.V1 <- simulateResiduals(fittedModel = AGOmonthNB1, n = 250, plot = TRUE)
AGOmonthNB2.V1 <- simulateResiduals(fittedModel = AGOmonthNB2, n = 250, plot = TRUE)

#Model reduction
a9 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
              + Year*IntPhase + PRCP + TAVG + MonthN + (1|Community/HouseID) 
              , family = "nbinom2", data = AGO)
summary(a9)
a10 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
               + Year*IntPhase + TAVG + MonthN + (1|Community/HouseID) 
               , family = "nbinom2", data = AGO)
summary(a10)
a11 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
               + Year*IntPhase + MonthN + (1|Community/HouseID) 
               , family = "nbinom2", data = AGO)
summary(a11)

AICtab(a9,a10,a11, base = T, weights = T, sort = FALSE)

r.squaredGLMM(a9, pj2014 = FALSE)

#This next code takes a lot of time, skip to the joint_test and plot effects if already tested
confint(a9, method = "uniroot")

car::Anova(a9, type = "III")

#Checking for multicollinearity
results <- check_collinearity(a9)
results
plot(results)

#Evaluating if the intervention on a weekly scale had an effect on the abundance of female Ae. aegypti.
a9.emm.s <- emmeans(a9, ~ IntPhase*Year)
plot(allEffects(a9))

pairs(a9.emm.s)
multcomp::cld(a9.emm.s)

#following [mean(Cont)-mean(Int)]+\-1.96{sqrt[SE(Cont)^2+SE(Int)^2]} to calculate z-confidence interval at (alpha=0.05)
#for comparison between Control and Intervention for the intervention period of 2017
(-0.937-0.088)-(1.96*sqrt((0.269^2)+(0.264^2))) #2017
(-0.937-0.088)+(1.96*sqrt((0.269^2)+(0.264^2)))

(-0.559+1.157)-(1.96*sqrt((0.264^2)+(0.268^2))) #2018
(-0.559+1.157)+(1.96*sqrt((0.264^2)+(0.268^2)))

write.csv(a9.emm.s, file = "a9emm.csv")

#GLMM for transiet effects delay 1 week

AGO$Delay1 <- relevel(as.factor(AGO$Delay1), ref = "Pre")

#Evaluating the best distribution for the dataset with count data
AGOlag1Ps <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                     + Year*Delay1 + PRCP + TAVG + (1|Community/HouseID) + (1|Week)
                     , family = "poisson", data = AGO)
AGOlag1NB1 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                      + Year*Delay1 + PRCP + TAVG + (1|Community/HouseID) + (1|Week)
                      , family = "nbinom1", data = AGO)
AGOlag1NB2 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                      + Year*Delay1 + PRCP + TAVG + (1|Community/HouseID) + (1|Week) 
                      , family = "nbinom2", data = AGO)

#Checking the AIC values for each model
AIC(AGOlag1Ps)
AIC(AGOlag1NB1)
AIC(AGOlag1NB2)

#Evaluation of the overdispersion of the dataset 
check_overdispersion(AGOlag1Ps)

#Visualization of the data
AGOlag1Ps.V1 <- simulateResiduals(fittedModel = AGOlag1Ps, n = 250, plot = TRUE)
AGOlag1NB1.V1 <- simulateResiduals(fittedModel = AGOlag1NB1, n = 250, plot = TRUE)
AGOlag1NB2.V1 <- simulateResiduals(fittedModel = AGOlag1NB2, n = 250, plot = TRUE)

#Model reduction
a12 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
               + Year*Delay1 + PRCP + TAVG + (1|Community/HouseID) + (1|Week)
               , family = "nbinom2", data = AGO)
summary(a12)
a13 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
               + Year*Delay1 + TAVG + (1|Community/HouseID) + (1|Week)
               , family = "nbinom2", data = AGO)
summary(a13)
a14 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
               + Year*Delay1 + (1|Community/HouseID) + (1|Week)
               , family = "nbinom2", data = AGO)
summary(a14)

AICtab(a12,a13,a14, base = T, weights = T, sort = FALSE)

r.squaredGLMM(a13, pj2014 = FALSE)

#This next code takes a lot of time, skip to the joint_test and plot effects if already tested
confint(a13, method = "uniroot")

car::Anova(a13, type = "III")

#Checking for multicollinearity
results <- check_collinearity(a13)
results
plot(results)

#Evaluating if the intervention on a weekly scale had an effect on the abundance of female Ae. aegypti.
a13.emm.s <- emmeans(a13, ~ Delay1*Year)
plot(allEffects(a13))

pairs(a13.emm.s)
multcomp::cld(a13.emm.s)

#following [mean(Cont)-mean(Int)]+\-1.96{sqrt[SE(Cont)^2+SE(Int)^2]} to calculate z-confidence interval at (alpha=0.05)
#for comparison between Control and Intervention for the intervention period of 2017
(-0.839-0.175)-(1.96*sqrt((0.342^2)+(0.342^2))) #2017
(-0.839-0.175)+(1.96*sqrt((0.342^2)+(0.342^2)))

(-0.471+1.157)-(1.96*sqrt((0.339^2)+(0.344^2))) #2018
(-0.471+1.157)+(1.96*sqrt((0.339^2)+(0.344^2)))

write.csv(a13.emm.s, file = "a13emm.csv")

#GLMM for transiet effects delay 2 weeks

AGO$Delay2 <- relevel(as.factor(AGO$Delay2), ref = "Pre")

#Evaluating the best distribution for the dataset with count data
AGOlag2Ps <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                     + Year*Delay2 + PRCP + TAVG + (1|Community/HouseID) + (1|Week)
                     , family = "poisson", data = AGO)
AGOlag2NB1 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                      + Year*Delay2 + PRCP + TAVG + (1|Community/HouseID) + (1|Week)
                      , family = "nbinom1", data = AGO)
AGOlag2NB2 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                      + Year*Delay2 + PRCP + TAVG + (1|Community/HouseID) + (1|Week) 
                      , family = "nbinom2", data = AGO)

#Checking the AIC values for each model
AIC(AGOlag2Ps)
AIC(AGOlag2NB1)
AIC(AGOlag2NB2)

#Evaluation of the overdispersion of the dataset 
check_overdispersion(AGOlag2Ps)

#Visualization of the data
AGOlag2Ps.V1 <- simulateResiduals(fittedModel = AGOlag2Ps, n = 250, plot = TRUE)
AGOlag2NB1.V1 <- simulateResiduals(fittedModel = AGOlag2NB1, n = 250, plot = TRUE)
AGOlag2NB2.V1 <- simulateResiduals(fittedModel = AGOlag2NB2, n = 250, plot = TRUE)

#Model reduction
a15 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
               + Year*Delay2 + PRCP + TAVG + (1|Community/HouseID) + (1|Week)
               , family = "nbinom2", data = AGO)
summary(a15)
a16 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
               + Year*Delay2 + TAVG + (1|Community/HouseID) + (1|Week)
               , family = "nbinom2", data = AGO)
summary(a16)
a17 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
               + Year*Delay2 + (1|Community/HouseID) + (1|Week)
               , family = "nbinom2", data = AGO)
summary(a17)

AICtab(a15,a16,a17, base = T, weights = T, sort = FALSE)

r.squaredGLMM(a16, pj2014 = FALSE)

#This next code takes a lot of time, skip to the joint_test and plot effects if already tested
confint(a16, method = "uniroot")

car::Anova(a16, type = "III")

#Checking for multicollinearity
results <- check_collinearity(a16)
results
plot(results)

#Evaluating if the intervention on a weekly scale had an effect on the abundance of female Ae. aegypti.
a16.emm.s <- emmeans(a16, ~ Delay2*Year)
plot(allEffects(a16))

plot_model(a16, type = "int")

pairs(a16.emm.s)
multcomp::cld(a16.emm.s)

#following [mean(Cont)-mean(Int)]+\-1.96{sqrt[SE(Cont)^2+SE(Int)^2]} to calculate z-confidence interval at (alpha=0.05)
#for comparison between Control and Intervention for the intervention period of 2017
(-0.9508-0.0596)-(1.96*sqrt((0.338^2)+(0.340 ^2))) #2017
(-0.9508-0.0596)+(1.96*sqrt((0.338^2)+(0.340 ^2)))

(-0.5915+1.3019)-(1.96*sqrt((0.336^2)+(0.342^2))) #2018
(-0.5915+1.3019)+(1.96*sqrt((0.336^2)+(0.342^2)))

write.csv(a16.emm.s, file = "a16emm.csv")

#Checking the weight of all the AIC models

AICtab(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17, base = T, weights = T, sort = FALSE)

# GLMM for coverage

AGOcovPs <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income + Placement + Income*Placement
                    + PRCP + TAVG + CovRate200 + Participation + (1|Week)
                    + (1 + CovRate200|Community), family = "poisson", data = AGO)
AGOcovNB1 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income + Placement + Income*Placement
                     + PRCP + TAVG + CovRate200 + Participation + (1|Week)
                     + (1 + CovRate200|Community), family = "nbinom1", data = AGO)
AGOcovNB2 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income + Placement + Income*Placement
                     + PRCP + TAVG + CovRate200 + Participation + (1|Week)
                     + (1 + CovRate200|Community), family = "nbinom2", data = AGO)

#Checking the AIC values for each model
AIC(AGOcovPs)
AIC(AGOcovNB1)
AIC(AGOcovNB2)

#Evaluation of the overdispersion of the dataset 
check_overdispersion(AGOcovPs)

#Visualization of the data
AGOcovPs.V1 <- simulateResiduals(fittedModel = AGOcovPs, n = 250, plot = TRUE)
AGOcovNB1.V1 <- simulateResiduals(fittedModel = AGOcovNB1, n = 250, plot = TRUE)
AGOcovNB2.V1 <- simulateResiduals(fittedModel = AGOcovNB2, n = 250, plot = TRUE)

#Coverage rate test for random slopes
a18 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income + Placement + Income*Placement
               + PRCP + TAVG + CovRate200 + Participation + (1|Week)
               + (1 + CovRate200|Community), family = "nbinom2", data = AGO)
summary(a18)
a19 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income + Placement + Income*Placement
               + PRCP + TAVG + CovRate200 + Participation + (1|Week)
               + (1|Community), family = "nbinom2", data = AGO)
summary(a19)

AICtab(a18,a19,base = T, weights = T, sort = FALSE)

r.squaredGLMM(a18, pj2014 = FALSE)

Anova3 <- anova(a18, a19)
Anova3

#This next code takes a lot of time, skip to the joint_test and plot effects if already tested
confint(a18, method = "uniroot")

car::Anova(a18, type = "III")

#Checking for multicollinearity
results <- check_collinearity(a18)
results
plot(results)

# GAMM with the variable of Coverage

#Testing the distribution of the dataset for the GAMM models

GamPsAll <- gamm4(AEAfemale ~ offset(log(daystrapping)) + Year + Income*Placement 
                  + s(Week) + s(CovRate200), random =~ (1|HouseID) + (1|Community), 
                  data = AGO, family=poisson)
summary(GamPsAll$mer)
summary(GamPsAll$gam)

GamNBAll <- gamm4(AEAfemale ~ offset(log(daystrapping)) + Year + Income*Placement
                  + s(Week) + s(CovRate200), random =~ (1|HouseID) + (1|Community), 
                  data = AGO, family=neg.bin(1))
summary(GamNBAll$mer)
summary(GamNBAll$gam)

anova(GamNBAll$gam)

GamNB2All <- gamm4(AEAfemale ~ offset(log(daystrapping)) + Year + Income*Placement 
                   + s(Week) + CovRate200, random =~ (1|HouseID) + (1|Community), 
                   data = AGO, family=neg.bin(1))

summary(GamNB2All$mer)
summary(GamNB2All$gam)

#Testing only trap density without houses
GamNB3All <- gamm4(AEAfemale ~ offset(log(daystrapping)) + Year + Income*Placement 
                   + s(Week) + s(Cov200m), random =~ (1|HouseID) + (1|Community), 
                   data = AGO, family=neg.bin(1))
summary(GamNB3All$mer)
summary(GamNB3All$gam)

plot(GamNB3All$gam, select = 1, rug=TRUE, residuals=TRUE, se = TRUE, shade=TRUE, ylim=c(0,0.20), trans = plogis, shift = coef(GamNB3All$gam)[1], shade.col="lightblue", xlab="Week", ylab=expression(paste(italic("f "),"(Female ", paste(italic("Aedes aegypti"), " captures)"))))
plot(GamNB3All$gam, select = 2, rug=TRUE, residuals=TRUE, se = TRUE, shade=TRUE, ylim=c(0,0.20), trans = plogis, shift = coef(GamNB3All$gam)[1], shade.col="lightblue", xlab="CovRate 200m", ylab=expression(paste(italic("f "),"(Female ", paste(italic("Aedes aegypti"), " captures)"))))

#AIC comparison for the two types of distributions
AICtab(GamNBAll$mer, GamPsAll$mer)

fvisgam(GamNBAll$gam, view = c("Week", "CovRate200"))
fvisgam(GamNBAll$gam, view = c("Week", "CovRate200"), zlim = c(-6,1.3), cond = list(IntPhase="Int", Year="2017"), main="Intervention 2017")
fvisgam(GamNBAll$gam, view = c("Week", "CovRate200"), zlim = c(-6,1.3), cond = list(IntPhase="Int", Year="2018"), main="Intervention 2018")

#Figure transformation for the plogis models
tiff("GAMM_plot.tiff",width=6.2, height=4.4, units="in", res=300)

par(mfrow=c(2,2), mai = c(1,1,0.2,0.1))

plot(GamNBAll$gam, select = 1, rug=TRUE, residuals=TRUE, se = TRUE, shade=TRUE, ylim=c(-3,3), pch = 20, cex = 0.1, shade.col="#F79680", xlab="Epidemiological Week", ylab=expression(paste(italic("f "),"(Epidemiological Week, 8.67)")))

par(xpd = TRUE)
legend(-19, 4.8, legend="A", bty="n", text.font = 2, cex=3, pt.cex = 2, xpd=T)

par(xpd = FALSE)
plot(GamNBAll$gam, select = 2, rug=TRUE, residuals=TRUE, se = TRUE, shade=TRUE, ylim=c(-2,2), pch = 20, cex = 0.1, shade.col="#F79680", xlab="CovRate 200m", ylab=expression(paste(italic("f "),"(CovRate 200m, 2.96)")))

par(xpd = TRUE)
legend(-1.1, 3.2, legend="B", bty="n", text.font = 2, cex=3, pt.cex = 2, xpd=T)

dev.off()

M3Pred <- predict(GamNBAll$gam, se=TRUE, type = 'response')

tiff("Fitedvalues_plot.tiff",width=6.2, height=2.1, units="in", res=300)
par(mfrow=c(1,2), mai = c(1,1,0.2,0.1))

coeff_bigger <- 1
plot(M3Pred$fit, AGO$AEAfemale,  cex=AGO$CovRate200*coeff_bigger, xlim=c(0,45), ylim=c(0,45), pch= 16, col=rgb(0,0,1,0.5), ylab=expression(paste("Female ", paste(italic("Aedes aegypti"), " captures"))), xlab=expression(paste("Female ", paste(italic("Aedes aegypti"), " fitted values"))))

plot(AGO$CovRate200, AGO$AEAfemale, cex=0.7, ylim=c(0,45), xlim=c(1,2.6), pch =20, col=rgb(red =0, green= 0, blue= 0, alpha =0.7),xlab= 'CovRate 200m', ylab="")

I <- order(AGO$CovRate200)
lines(AGO$CovRate200[I], M3Pred$fit[I], lwd=1, col=rgb(red =0, green= 0, blue= 1, alpha =0.8))

dev.off()

lines(AGO$CovRate200[I], M3Pred$fit[I] + 2*M3Pred$se.fit[I], lty=2, lwd=2, col=rgb(red =1, green= 0, blue= 0, alpha =0.4))
lines(AGO$CovRate200[I], M3Pred$fit[I] - 2*M3Pred$se.fit[I], lty=2, lwd=2, col=rgb(red =1, green= 0, blue= 0, alpha =0.4))

head(predict(GamNBAll$gam, type='response'))

confint(GamNBAll$gam, parm = 'PlacementOUT', level=0.5)

p <- predict(GamNBAll$gam, se.fit = TRUE, type="link")

predict(GamNBAll$gam)

upr <- p$fit + (2 * p$se.fit)
lwr <- p$fit - (2 * p$se.fit)

upr1 <- GamNBAll$gam$family$linkinv(upr)
lwr1 <- GamNBAll$gam$family$linkinv(lwr)

#extract predicted values by case

DataAGO1 <- AGO[AGO$CovRate200 > 0.8,]
AGO1.1 <- DataAGO1[,"Pi"]
AGO1.1

GammTest <- gamm4(AEAfemale ~ offset(log(daystrapping)) + Year + Income*Placement
                  + s(Week) + s(CovRate200), random =~ (1|HouseID) + (1|Community), 
                  data = AGO1, family=neg.bin(1))
summary(GammTest$mer)
summary(GammTest$gam)

M4Pred <- predict(GammTest$gam, se=TRUE, type = 'response')

plot(GammTest$gam, select = 2, rug=TRUE, residuals=TRUE, se = TRUE,ylim=c(-2,2),shade=TRUE, pch = 20, cex = 0.1, shade.col="#F79680", xlab="CovRate 200m", ylab=expression(paste(italic("f "),"(Female ", paste(italic("Aedes aegypti"), " captures)"))))

AGO1 <- dplyr::filter(AGO, IntPhase=="Int")

plot(AGO1$CovRate200, AGO1$AEAfemale, cex=0.5, ylim=c(0,70), pch =20, main= 'Negative binomial GAMM', xlab= 'CovRate 200m', ylab= 'Aedes aegypti females')

I <- order(AGO1$CovRate200)
lines(AGO1$CovRate200[I], M4Pred$fit[I], lwd=1, col= "brown")
lines(AGO1$CovRate200[I], M4Pred$fit[I] + 2*M4Pred$se.fit[I], lty=2, lwd=2, col=gray(0.8))
lines(AGO1$CovRate200[I], M4Pred$fit[I] - 2*M4Pred$se.fit[I], lty=2, lwd=2, col=gray(0.8))
for(i in 1:52) {
  y <- rnbinom(100, size =11.8, mu = M4Pred$fit[i])
  points(rep(AGO1$CovRate200[i], 100), y, cex =0.5)
}

AGO2 <- dplyr::filter(AGO1, Placement=="OUT")
plot(AGO2$CovRate200, AGO2$AEAfemale, cex=0.5, pch =20, main= 'Negative binomial GAMM', xlab= 'CovRate 200m', ylab= 'Aedes aegypti females')

I <- order(AGO2$CovRate200)
lines(AGO2$CovRate200[I], M3Pred$fit[I], lwd=1, col= "brown")
lines(AGO2$CovRate200[I], M3Pred$fit[I] + 2*M3Pred$se.fit[I], lty=2, lwd=2, col=gray(0.8))
lines(AGO2$CovRate200[I], M3Pred$fit[I] - 2*M3Pred$se.fit[I], lty=2, lwd=2, col=gray(0.8))
for(i in 1:52) {
  y <- rnbinom(100, size =11.8, mu = M4Pred$fit[i])
  points(rep(AGO2$CovRate200[i], 100), y, cex =0.5)
}

#IT methodology based on factors associated with mosquito interventions
#Immediate effects
IT1 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
                  + Year*IntPhase + (1|Community/HouseID) 
                  + (1|Week), family = "nbinom2", data = AGO)

summary(IT1)

IT2 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income + Placement 
               + Year*IntPhase + (1|Community/HouseID) 
               + (1|Week), family = "nbinom2", data = AGO)
summary(IT2)

IT3 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Placement 
               + Year + IntPhase + (1|Community/HouseID) 
               + (1|Week), family = "nbinom2", data = AGO)
summary(IT3)

AICtab(IT1,IT2,IT3, base = T, weights = T, sort = FALSE)
#Reduced weekly
IT4 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
               + Year*IntPhase + Week + (1|Community/HouseID) 
               , family = "nbinom2", data = AGO)

summary(IT4)

IT5 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income + Placement 
               + Year*IntPhase + Week + (1|Community/HouseID) 
               , family = "nbinom2", data = AGO)
summary(IT5)

IT6 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Placement 
               + Year + IntPhase + Week + (1|Community/HouseID) 
               , family = "nbinom2", data = AGO)
summary(IT6)

AICtab(IT4,IT5,IT6, base = T, weights = T, sort = FALSE)
#Reduced monthly
IT7 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
               + Year*IntPhase + MonthN + (1|Community/HouseID) 
               , family = "nbinom2", data = AGO)

summary(IT7)

IT8 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income + Placement 
               + Year*IntPhase + MonthN + (1|Community/HouseID) 
               , family = "nbinom2", data = AGO)
summary(IT8)

IT9 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Placement 
               + Year + IntPhase + MonthN + (1|Community/HouseID) 
               , family = "nbinom2", data = AGO)
summary(IT9)

AICtab(IT7,IT8,IT9, base = T, weights = T, sort = FALSE)
#Delayed 1 week
IT10 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
               + Year*Delay1 + (1|Community/HouseID) 
               + (1|Week), family = "nbinom2", data = AGO)

summary(IT10)

IT11 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income + Placement 
               + Year*Delay1 + (1|Community/HouseID) 
               + (1|Week), family = "nbinom2", data = AGO)
summary(IT11)

IT12 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Placement 
               + Year + Delay1 + (1|Community/HouseID) 
               + (1|Week), family = "nbinom2", data = AGO)
summary(IT12)

AICtab(IT10,IT11,IT12, base = T, weights = T, sort = FALSE)
#Delayed 2 weeks
IT13 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
               + Year*Delay2 + (1|Community/HouseID) 
               + (1|Week), family = "nbinom2", data = AGO)

summary(IT13)

IT14 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income + Placement 
               + Year*Delay2 + (1|Community/HouseID) 
               + (1|Week), family = "nbinom2", data = AGO)
summary(IT14)

IT15 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Placement 
               + Year + Delay2 + (1|Community/HouseID) 
               + (1|Week), family = "nbinom2", data = AGO)
summary(IT15)

AICtab(IT13,IT14,IT15, base = T, weights = T, sort = FALSE)

#Comparing the best fit models from the IT selection approach
AICtab(IT1,IT4,IT7,IT10,IT13, base = T, weights = T, sort = FALSE)


#Best fit models from the stepwise approach
a3 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
              + Year*IntPhase + (1|Community/HouseID) + (1|Week), 
              family = "nbinom2", data = AGO)

a6 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
              + Year*IntPhase + PRCP + TAVG + Week + (1|Community/HouseID) 
              , family = "nbinom2", data = AGO)

a9 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
              + Year*IntPhase + PRCP + TAVG + MonthN + (1|Community/HouseID) 
              , family = "nbinom2", data = AGO)

a13 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
               + Year*Delay1 + TAVG + (1|Community/HouseID) + (1|Week)
               , family = "nbinom2", data = AGO)

a16 <- glmmTMB(AEAfemale ~ offset(log(daystrapping)) + Income*Placement 
               + Year*Delay2 + TAVG + (1|Community/HouseID) + (1|Week)
               , family = "nbinom2", data = AGO)

AICtab(a3,a6,a9,a13,a16, base = T, weights = T, sort = FALSE)

summary(a16)
anova(a16, IT10)




#Figure transformation for the plogis models
tiff("GAMM_plotS.tiff",width=9.2, height=6.4, units="in", res=300)

par(mfrow=c(2,2), mai = c(0.7,0.7,0.1,0.1), oma = c(0,0,0,0))

plot(GamNBAll$gam, select = 1, rug=TRUE, residuals=TRUE, se = TRUE, shade=TRUE, ylim=c(-3,3), pch = 20, cex = 0.1, shade.col="#F79680", xlab="Epidemiological Week", ylab=expression(paste(italic("f "),"(Epidemiological Week, 8.67)")))

plot(GamNBAll$gam, select = 2, rug=TRUE, residuals=TRUE, se = TRUE, shade=TRUE, ylim=c(-2,2), pch = 20, cex = 0.1, shade.col="#F79680", xlab="CovRate 200m", ylab=expression(paste(italic("f "),"(CovRate 200m, 2.96)")))

coeff_bigger <- 1
plot(M3Pred$fit, AGO$AEAfemale,  cex=AGO$CovRate200*coeff_bigger, xlim=c(0,45), ylim=c(0,45), pch= 16, col=rgb(0,0,1,0.5), ylab=expression(paste("Female ", paste(italic("Aedes aegypti"), " captures"))), xlab=expression(paste("Female ", paste(italic("Aedes aegypti"), " fitted values"))))

plot(AGO$CovRate200, AGO$AEAfemale, cex=0.7, ylim=c(0,45), xlim=c(1,2.6), pch =20, col=rgb(red =0, green= 0, blue= 0, alpha =0.7),xlab= 'CovRate 200m', ylab="")

I <- order(AGO$CovRate200)
lines(AGO$CovRate200[I], M3Pred$fit[I], lwd=1, col=rgb(red =0, green= 0, blue= 1, alpha =0.8))

dev.off()