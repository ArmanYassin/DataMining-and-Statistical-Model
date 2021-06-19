rm(list = ls())

library(car)
library(rrcov)
library(robustbase)
library(rospca) 
library(sp)
library(cluster)
library(rgl)

#analyse reele data  wijn karakteristieken

#exploratory phase

wijn <- read.csv("wine567(1).csv", header = TRUE)
dim(wijn)
head(wijn)
summary(wijn)
sample1 <- sample(which(wijn[,12]==5),200)
sample2 <- sample(which(wijn[,12]==6),200)
sample3 <- sample(which(wijn[,12]==7),200)
mydata <- wijn[c(sample1,sample2,sample3),]


dim(mydata)
head(mydata)
summary(mydata)

datazq <- subset( mydata, select = -quality  )
#datazq <- subset( datazqa, select = -pH  )
dim(datazq)
plot(datazq)
matplot(t(datazq),lty=1,type="l")

datazqSt <- scale(datazq, center = TRUE, scale = TRUE)

boxplot(datazq[c(1:3)])
for (n in 1 : 11) {
  boxplot(datazq[,n],main = paste0("boxplot (n = ", n, ")"))
}
for (n in 1 : 11) {
  qqnorm(datazq[,n], main = paste0("Normal Q-Q Plot (n = ", n, ")"))
  qqline(datazq[,n])
}
for (n in 1 : 11) {
  hist(datazq[,n], main = paste0("HistVerd (n = ", n, ")"))
}


#Welke variabele kandidaat is voor transformatie
qqnorm(datazq[,8])
qqline(datazq[,8])
shapiro.test(datazq[,3])

bc2X1 <- boxCox(datazq[,10] ~ 1)
lambda2 <- bc2X1$x[which.max(bc2X1$y)]
lambda2

lambdaBC <- powerTransform(datazq[,10] ~ 1)
lambdaBC
bc <- bcPower(datazq[,11], -1)
qqnorm(bc)
qqline(bc)
shapiro.test(bc)
shapiro.test(log(datazq[,2]))


plot(density(datazq[,7]))
boxplot(datazq[,7])
hist(datazq[,7])
shapiro.test(log(datazq$volatile.acidity))

#4,8
sugarDensity <- cbind(datazq[,4],datazq[,8])
densityAlcohol <- cbind(datazq[,8],datazq[,11])
MDda <- mahalanobis(densityAlcohol, colMeans(densityAlcohol), var(densityAlcohol))
MDsd <- mahalanobis(sugarDensity, colMeans(sugarDensity), var(sugarDensity))
qqplot(qchisq(ppoints(600), df = 2), MDda) 
qqline(MDda, distribution = function(pp) {qchisq(pp, df = 2)}, col = "red", lwd = 2)
abline(h = qchisq(0.975, df = 2), col = "gray", lty = 2, lwd = 2)

#Verder gaan met alles behalve pH 
summary(datazq)
str(datazq)
apply(datazq,2, sd)
Sigm <- cov(datazq)
Rr <- cor(datazq)


dataCov <- PcaClassic(datazq, scale = FALSE)
dataCor <- PcaClassic(datazq, scale = TRUE)  
str(dataCor)
dataCor@loadings[, 1]
dataCor@eigenvalues

#??dataCov@loadings[, 1] * sqrt(dataCov@eigenvalues[1]) / sqrt(diag(Sigm)) geen conclusie? 
dataCor@loadings[, 1] * sqrt(dataCor@eigenvalues[1])

#wat kan kan hier gepresenteerd worden? wat is verloren? 
summary(dataCor)
str(dataCor)
screeplot(dataCor, type = "lines")
wineQuality <- factor(mydata$quality)
plot(dataCor, col = wineQuality)
legend("topright", bty = "n", legend = c("5","6","7"),pch = 16,col=c("black", "red", "green"), cex=0.8)
biplot(dataCor)

#hoeveel componenten k is ideaal voor projectie? 
load_datazq_3 <- PcaClassic(datazq, k = 4, scale = TRUE) 
load_datazq_3Rob <- PcaHubert(datazq, k = 4, scale = TRUE) 

#interpreteer eerste twee loadings, maw eerste twe pricipaal componenten
dataCor@loadings[,1]
dataCor@loadings[,2]

scr <- load_datazq_3@scores

#wat is de correlatie tussen load en oorspronkelijke variabele? 
#indien sterke correlatie --> meest informatie is gevat in die variabele 
#load_datazq_3[,1]*sqrt(load_data_3@eigenvalues[1])
#load_datazq_3[,2]*sqrt(load_data_3@eigenvalues[2])

pairs(load_datazq_3)
plot(load_datazq_3 )
plot(load_datazq_3Rob )

#dit is je nieuwe dataset, wat verwachte je? ==> Geen correlatie!!

#robust PCA??

#Clustering. Opgelet, het is al geschaald

summary(scr[,1:2])
plot3d(scr[,1:2])
scores2scale <- scale(scr[,1:2], center=TRUE, scale=TRUE)
apply(scr[,1:2], 2, sd)
kmeans_load2 <- kmeans(scr[,1:2], centers=2)
kmeans_load2
plot(silhouette(kmeans_load2$cluster, dist(scr[,1:2])))
pam_scr <- pam(scr[,1:2], k = 2)
str(pam_scr)
plot(silhouette(pam(scr[,1:2], k = 2)))
clusplot(pam_scr)
plot(scr[,1:2], col = wineQuality)


scr2 <- scr[,1:2]
Cl1 <- subset(scr2, pam_scr$clustering==1)
Cl2 <- subset(scr2, pam_scr$clustering==2)
dim(Cl1)

plot(scr2, col=wineQuality)
rd <- sqrt(qchisq(0.95, df = 2))
mcd1 <- covMcd(scr2, alpha = 0.5)
mcd2 <- covMcd(scr2, alpha = 0.75)
car::ellipse(center = mcd1$center,            shape = mcd1$cov,          radius = rd, col = "blue")
car::ellipse(center = mcd2$center,            shape = mcd2$cov,          radius = rd, col = "magenta")
plot(mcd2)

#Cl2 bevat meest observatie
plot(Cl2, col=wineQuality)
mcd1Cl2 <- covMcd(Cl2, alpha = 0.5)
mcd2Cl2 <- covMcd(Cl2, alpha = 0.75)
car::ellipse(center = mcd1Cl2$center,            shape = mcd1Cl2$cov,          radius = rd, col = "blue")
car::ellipse(center = mcd2Cl2$center,            shape = mcd2Cl2$cov,          radius = rd, col = "magenta")
plot(mcd2Cl2)


p <- ellipse(center = mcd1Cl2$center,            shape = mcd1Cl2$cov,          radius = rd, col = "blue")
pInEllipse <- data.frame(Cl2, in.p = as.logical(point.in.polygon(Cl2[,1], Cl2[,2], p[,1], p[,2])))
Cl2_outliersFiltered <- subset(Cl2, pInEllipse[,3]==TRUE)
