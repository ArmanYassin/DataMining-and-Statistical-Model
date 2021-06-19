rm(list = ls())

install.packages("car")
install.packages("rrcov")
install.packages("robustbase")
install.packages("rospca")
install.packages("sp")
install.packages("cluster")
install.packages("rgl")
install.packages("DAAG")
install.packages("leaps")
install.packages("SMPracticals")

library(car)
library(rrcov)
library(robustbase)
library(rospca) 
library(sp)
library(cluster)
library(rgl)
library(MASS)
library(DAAG)
library(leaps)
library(pls)
library(SMPracticals)

set.seed(0470345)
data <- read.table("diabetes.txt", header = TRUE)
selVec <- c(sample(1:dim(data)[1],300))
Xstudent <- data[selVec,]

#This project consists of a regression analysis of the dia-
#  betes dataset. The data is based on 442 individuals. Ten baseline variables were obtained
#for each of these 442 diabetes patients as well as the response of interest (y), a quantitative
#measure of disease progression one year after baseline. The ten baseline variables are age,
#sex (male= 1, female= 0), body mass index, average blood pressure, and six blood serum
#measurements.

dim(Xstudent)
head(Xstudent)
summary(Xstudent)
str(Xstudent)

Xstudent$SEX <- as.factor(Xstudent$SEX)
str(Xstudent)
plot(Xstudent)

boxplot(Xstudent[,11])
for (n in 3 : 11) {
  boxplot(Xstudent[,n], main = paste0("bx (n = ", n, ")"))
}
for (n in 3 : 11) {
  plot(density(Xstudent[,n]), main = paste0("dens (n = ", n, ")"))
}
plot(Xstudent[,5],Xstudent[,6])
robAn <- cbind(Xstudent$S1,Xstudent$Y)
plot(robAn, xlim = c(-10, 350 ), ylim = c(-20, 350))
rd <- sqrt(qchisq(0.95, df =2))
mcd1X <- covMcd(robAn, alpha = 0.5)
mcd2X <- covMcd(robAn, alpha = 0.75)
car::ellipse(center = mcd1X$center,            shape = mcd1X$cov,          radius = rd, col = "blue")
car::ellipse(center = mcd2X$center,            shape = mcd2X$cov,          radius = rd, col = "magenta")
plot(mcd1X)
plot(mcd2X)



boxplot((Xstudent[,11])^0.25)
outliers = which(Xstudent$S2 > boxplot(Xstudent$S2)$stats[5, ])
outliers[3]
Xstudent2 = Xstudent[-outliers[3], ]
boxplot(Xstudent2[,6])
plot(Xstudent2)
plot(Xstudent)

#multicoll

Xstudent3 <-subset(Xstudent2, select = -S2) 
Rr <- cor(subset(Xstudent3, select = -SEX))
Rd <- diag(solve(Rr))
max(Rd)
mean(Rd)

ev <- eigen(cor(subset(Xstudent2, select = -SEX)))
ev$values / sum(ev$values)              # some small values: indication of multicollinearity
sqrt( max(ev$values) / min(ev$values) )
plot(density(Xstudent2$Y^0.4))



boxcox(Y ~ ., data = Xstudent2)
lambdaBC <- powerTransform(Xstudent2$AGE ~ 1)
lambdaBC
qqnorm((Xstudent2$Y))
qqline((Xstudent2$Y))
shapiro.test(log(Xstudent2$S5))

#lukt niet 
bcBMI <- bcPower(Xstudent2$BMI, -0.5)
qqnorm(bcBMI)
qqline(bcBMI)
shapiro.test(Xstudent$BMI^-0.5)

#
bcBP <- bcPower(Xstudent2$BP, -0.75)
qqnorm(bcBP)
qqline(bcBP)
shapiro.test(bcBP)

bc2X2 <- boxCox(Xstudent$Y ~ 1)
lambdaX2_2 <- bc2X2$x[which.max(bc2X2$y)]
lambdaX2_2

#reg
Ylm <- lm(Y ~ ., data = Xstudent3)
#model.matrix(Ylm)
summary(Ylm)

#Model variabel selection 
acal <- Xstudent3[1:200, ]     # calibration data
aval <- Xstudent3[-c(1:200), ]

Yall <- regsubsets(Y ~ . , data = acal, nbest = 1, nvmax = 9)#9 predictor - 1 categorical + 1 dummy var

str(summary(Yall))
summary(Yall, matrix.logical = TRUE)
plot(Yall, scale = "adjr2")

#r2
summary(Yall)$which
numvar <- as.numeric(row.names(summary(Yall)$which))
numvar
rsq <- summary(Yall)$rsq
plot(numvar, rsq, pch = 16)
chosenR2 <- which.max(rsq)
chosenR2

#adjr2
adjrsq <- summary(Yall)$adjr2
plot(numvar, adjrsq, pch = 16) 
chosenadjR2 <- which.max(adjrsq)
chosenadjR2
summary(Yall)$which[chosenadjR2, ]
max(adjrsq)
plot(Yall, scale = "adjr2")

#AIC
aic <- summary(Yall)$bic + (2 - log(nrow(acal))) * numvar # BIC => AIC conversion
plot(numvar, aic, pch = 16)
chosenaic <- which.min(aic)
chosenaic
summary(Yall)$which[chosenaic, ]

#Cp
cp <- summary(Yall)$cp
plot(numvar, cp, pch = 16)
chosencp <- which.min(cp)
chosencp
summary(Yall)$which[chosencp, ]

#PRESS
designcal <- model.matrix(lm(Y ~ . , data = acal))
designval <- model.matrix(lm(Y ~ . , data = aval))

#r2
summary(Yall)$which[chosenR2, ]
selR2 <- names(which( summary(Yall)$which[chosenR2, !colnames(summary(Yall)$which) %in% c("(Intercept)")] == 1))
acal2 <- cbind(designcal[, selR2], acal[, "Y", drop = FALSE])
aval2 <- cbind(designval[, selR2], aval[, "Y", drop = FALSE])

mod2 <- lm(Y ~ . , data = acal2)

cvpress2 <- press(mod2)
mcvpress2 <- cvpress2/nrow(acal2) 
mcvpress2

#MSEP on valset
fittedval2 <- predict(mod2, newdata = aval2, interval = "prediction")[, 1]
sum((fittedval2-aval2$Y)^2)/nrow(aval2) #2625.737

#adjr2
summary(Yall)$which[chosenadjR2, ]
seladj <- names(which( summary(Yall)$which[chosenadjR2, !colnames(summary(Yall)$which) %in% c("(Intercept)")] == 1))
acal3 <- cbind(designcal[, seladj], acal[, "Y", drop = FALSE])
aval3 <- cbind(designval[, seladj], aval[, "Y", drop = FALSE])

mod3 <- lm(Y ~ . , data = acal3)

cvpress3 <- press(mod3)
mcvpress3 <- cvpress3/nrow(acal3) 
mcvpress3
#MSEP on valset
fittedval3 <- predict(mod3, newdata = aval3, interval = "prediction")[, 1]
sum((fittedval3-aval3$Y)^2)/nrow(aval3) #2560.675

#AIC

summary(Yall)$which[chosenaic, ]
selaic <- names(which( summary(Yall)$which[chosenaic, !colnames(summary(Yall)$which) %in% c("(Intercept)")] == 1))
acal4 <- cbind(designcal[, selaic], acal[, "Y", drop = FALSE])
aval4 <- cbind(designval[, selaic], aval[, "Y", drop = FALSE])

mod4 <- lm(Y ~ . , data = acal4)
summary(mod4)

cvpress4 <- press(mod4)
mcvpress4 <- cvpress4/nrow(acal4) 
mcvpress4
#MSEP on valset
fittedval4 <- predict(mod4, newdata = aval4, interval = "prediction")[, 1]
sum((fittedval4-aval4$Y)^2)/nrow(aval4) #same as adj??



#stepwise regression
mod7 <- lm(Y ~ 1, data = acal2)


## CREATE FORMULA
listofvars <- colnames(designcal)
listofvars <- listofvars[2:length(listofvars)]  # remove intercept
fullformula <- as.formula(paste("~", paste(listofvars, collapse = "+")))


## FORWARDS, goed motiveren waarom!
modEF <- stepAIC(mod7, list(lower = ~ 1, upper = fullformula), direction = "forward", k = 2)
summary(modEF)$coeff[, 1]#7 variabele met 

#stepwize
modES <- stepAIC(mod7, list(lower = ~ 1, upper = fullformula), direction = "both", k = 2)
summary(modES)$coeff[, 1]
# stepwize en adr, aic en cp geven dezelfde gegevens

#gauss markov
pBMI<- Xstudent3$BMI^-0.5
#Y ~ SEX + BMI + BP + S1  +S4 +S5
fit_Y <- lm(Y ~ SEX + pBMI+ BP + S1 +S4  +log(S5) , data = Xstudent3)
e <- fit_Y$residuals # = acal4$Y - fitted(mod4)
plot(e, xlab = "index", ylab = "Residual")                      #residual vs index
plot(fitted(fit_Y), e, xlab = "fitted", ylab = "Residual")     #residual vs fitted, fitted() extract fitted values
plot(Xstudent3$BP, e, ylab = "Residuals", xlab = "BP")  #residual vs independent var
plot(fit_Y)
abline(h = 0, lty = 2)   
qqnorm(e)
qqline(e)
shapiro.test(e)

es <- stdres(fit_Y)
pairs(fitted(fit_Y), es, xlab = "fitted", ylab = "Residual",)
plot(es, xlab = "index", ylab = "Standardized Residual") #uitleg(m4749O6)!!
plot(fitted(fit_Y), es, xlab = "fitted", ylab = "Std Residual") 
abline(h=-2.5, lty=2)
abline(h=2.5, lty=2)
qqnorm(es)
qqline(es)
shapiro.test(es)
plot(Xstudent3$BP,abs(e))
cor(Xstudent3$BP,e)

#weighted regression
stdev <- lm(abs(e)~XsS5)
summary(stdev)


wts <- 1/fitted(lm(abs(fit_Y2$residuals) ~ Xstudent3$S5))^2
fit_Y2 <- lm(Y^0.2 ~ SEX + BMI + BP +S4 +S5 , data = Xstudent3, weights = wts)
plot(bcBMI, residuals(fit_Y2)*sqrt(wts))
plot(fitted(fit_Y2), residuals(fit_Y2))
summary(lm(abs(e) ~ fitted(fit_Y)))
plot(fit_Y2)


summary(fit_Y)
summary(fit_Y)$r.squared
summary(fit_Y)$adj.r.squared
anova(fit_Y)
summary(aov(fit_Y))
NYinf <- lm.influence(fit_Y)
h <- NYinf$hat #diagonal of the 'hat' matrix
es1 <- es - e/(sqrt(sum(e^2)/(299-6))*(1-h)^.5)


robFit <- ltsReg(Y ~SEX +  BMI + BP + S1  +S4 +(S5) , data = Xstudent)
summary(robFit)
plot(robFit, id.n = 6)
plot(robFit, which = "rdiag")
rRes<-robFit$resid
robDis <- robFit$RD
vO <- identify(robDis, rRes , row.names(Xstudent3))
Xstudent4 = Xstudent3[-vO,]
fit_Y4 <- lm(Y ~ SEX + BMI + BP  +S5 , data = Xstudent4)
summary(fit_Y4)
plot(fitted(fit_Y4), residuals(fit_Y4))
summary(fit_Y)$adj.r.squared
ano <- anova(fit_Y4)
ano[nrow (ano), "Mean Sq"]