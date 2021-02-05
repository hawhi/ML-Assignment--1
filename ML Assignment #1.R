library(faraway)
head(mammalsleep)
dim(mammalsleep)
sum(is.na(mammalsleep))
sum(is.na(mammalsleep$dream))
sum(is.na(mammalsleep$nondream))
sum(is.na(mammalsleep$lifespan))
?mammalsleep

library(skimr)
skim(mammalsleep)
summary(mammalsleep)
pairs(mammalsleep)

library(faraway)
head(mammalsleep)
dim(mammalsleep)

?mammalsleep

library(skimr)
skim(mammalsleep)
summary(mammalsleep)



#<the one difference>
cleandata<-subset(mammalsleep, select = -c(dream, nondream))
cleandata<-na.omit(cleandata)
#</the one difference>

skim(cleandata)
summary(cleandata)
pairs(cleandata)

#full model 2
fullmodel<-lm(sleep ~., data=cleandata)
summary(fullmodel)

#model search method - BSS

library(leaps)
best = regsubsets(sleep~., cleandata)
summary(best)
results<-summary(best)
names(results)
RSS<-results$rss
r2<-results$rsq
Cp<-results$cp
BIC<-results$bic
Adjr2<-results$adjr2
cbind(RSS, r2, Cp, BIC, Adjr2)

par(mfrow = c(1, 2))
plot(RSS, xlab = "No. Predictors", ylab = "RSS", type = "l", lwd = 2)
plot(r2, xlab = "No. Predictors", ylab = "R-square", type = "l", lwd = 2)

which.min(Cp)
which.min(BIC)
which.max(Adjr2)

#all say 3

par(mfrow = c(1, 3))
plot(Cp, xlab = "Number of Predictors", ylab = "Cp", type = 'l', lwd = 2)
points(3, Cp[3], col = "red", cex = 2, pch = 8, lwd = 2)
plot(BIC, xlab = "Number of Predictors", ylab = "BIC", type = 'l', lwd = 2)
points(3, BIC[3], col = "red", cex = 2, pch = 8, lwd = 2)
plot(Adjr2, xlab = "Number of Predictors", ylab = "Adjusted RSq", type = "l", lwd = 2)
points(3, Adjr2[3], col = "red", cex = 2, pch = 8, lwd = 2)


par(mfrow = c(1,1))
plot(best, scale = "bic")
coef(best,3) # Cp

#gestation, predation and danger

#pcr
#install.packages("pls")
library(pls)

set.seed(420) 
pcr_model<-pcr(sleep~., data=cleandata, scale=TRUE, validation="CV")
summary(pcr_model)

?validationplot
validationplot(pcr_model, val.type = "MSEP")

min.pcr = which.min(MSEP(pcr_model)$val[1,1, ] ) - 1
min.pcr

#it suggests 3 or 4 comps

?mammalsleep

coef(pcr_model, ncomp = min.pcr)

head(predict(pcr_model, ncomp = min.pcr))

#regularisation paths

coef.mat = matrix(NA, 7, 7)
for(i in 1:7){
  coef.mat[,i] = pcr_model$coefficients[,,8-i]
}
plot(coef.mat[1,], type = 'l', ylab = 'Coefficients', xlab = 'Number of components', ylim = c(min(coef.mat), max(coef.mat)))
for(i in 2:7){
  lines(coef.mat[i,], col = i)
}
abline(v = min.pcr, lty = 3)


#ridge regression w/o CV
#ridge regression is good for multicollinearity
#/maybe also if p is close to n (which is the case here?)
  
install.packages("glmnet")
#library(glmnet)
?glmnet
y=cleandata$sleep
x<-model.matrix(sleep~., data=cleandata) [,-1] #excluding intercept
head(x)
ridge<-glmnet(x,y, alpha = 0)
names(ridge)
ridge$lambda
dim(ridge$beta)
ridge$beta[,1:3]
coef(ridge)[,1:3]
par(mfrow=c(1,2))
plot(ridge, xvar = 'lambda')
plot(ridge, xvar = 'norm')

#cv for lambda

set.seed(420)
ridge.cv<-cv.glmnet(x, y, alpha=0)
ridge.cv$lambda.min
#0.4502597
ridge.cv$lambda.1se
#0.5423394

round(cbind(
  coef(ridge.cv, s = "lambda.min"),
  coef(ridge.cv, s = "lambda.1se")), 3)

par(mfrow=c(1,2))
plot(ridge.cv)
plot(ridge, xvar = "lambda")
abline(v = log(ridge.cv$lambda.min), lty = 3)
abline(v = log(ridge.cv$lambda.1se), lty = 3)

#Workshop 3 notes: 'Important note: Any method/technique that relies on validation/cross-validation is subject to variability.
#If we re-run this code under a different random seed results will change.'
#talk about bias-variance

#post-hoc analysis on RR - 'sparsification'

beta.1se = coef(ridge.cv, s = 'lambda.1se')[-1]
rank.coef = sort(abs(beta.1se), decreasing=TRUE)
rank.coef

repetitions = 50
cor.1 = c()
cor.2 = c()
set.seed(420)

?sample
dim(cleandata)



for(i in 1:repetitions){
  # Step (i) data splitting
  training.obs = sample(1:42, 28)
  y.train = cleandata$sleep[training.obs]
  x.train = model.matrix(sleep~., cleandata[training.obs, ])[,-1]
  y.test = cleandata$sleep[-training.obs]
  x.test = model.matrix(sleep~., cleandata[-training.obs, ])[,-1]
  # Step (ii) training phase
  ridge.train1 = cv.glmnet(x.train, y.train, alpha = 0)
  # Step (iii) here we find which coefficients to set equal to zero
  ind = which(abs(coef(ridge.train1)) < 0.02)
  beta.sparse = coef(ridge.train1)
  beta.sparse[ind] = 0
  # Step (iv) generating predictions
  predict.1 = predict(ridge.train1, x.test, s = 'lambda.1se')
  predict.2 = cbind(1, x.test) %*% as.numeric(beta.sparse)
  # Step (v) evaluating predictive performance
  cor.1[i] = cor(y.test, predict.1)
  cor.2[i] = cor(y.test, predict.2)
}
boxplot(cor.1, cor.2, names = c('Standard','Sparse'), ylab = 'Test correlation', col = 3)

#test correlation went down a little after sparsification

repetitions = 100
cor.cp = c()
cor.bic = c()
cor.adjr2 = c()


predict.regsubsets = function(object, newdata, id, ...){
  form = as.formula(object$call[[2]])
  mat = model.matrix(form, newdata)
  coefi = coef(object, id = id)
  xvars = names(coefi)
  mat[, xvars]%*%coefi
}


set.seed(1) # for the results to be reproducable
for(i in 1:repetitions){
  # Step (i) data splitting
  training.obs = sample(1:51, 34)
  cd.train = cleandata[training.obs, ]
  cd.test = cleandata[-training.obs, ]
  dim(cd.train)
  dim(cd.test)
  # Step (ii) training phase
  regfit.train = regsubsets(sleep~., data = cd.train, nvmax = 9)
  min.cp = which.min(results$cp)
  min.bic = which.min(results$bic)
  max.adjr2 = which.max(results$adjr2)
  # Step (iii) generating predictions
  predict.cp = predict.regsubsets(regfit.train, cd.test, min.cp)
  predict.bic = predict.regsubsets(regfit.train, cd.test, min.bic)
  predict.adjr2 = predict.regsubsets(regfit.train, cd.test, max.adjr2)
  # Step (iv) evaluating predictive performance
  cor.cp[i] = cor(cd.test$sleep, predict.cp)
  cor.bic[i] = cor(cd.test$sleep, predict.bic)
  cor.adjr2[i] = cor(cd.test$sleep, predict.adjr2)
}
boxplot(cor.cp, cor.bic, cor.adjr2, names = c('Cp','BIC','adjRsq'), ylab = 'Test correlation', col = 2)


#install.packages("pls")
library(pls)

repetitions = 100
cor.pcr<-c()


set.seed(1) 
for(i in 1:repetitions){
  #split
  training.obs = sample(1:51, 34)
  cd.train = cleandata[training.obs, ]
  cd.test = cleandata[-training.obs, ]
  #train
  pcr_model<-pcr(sleep~., data=cd.train, scale=TRUE, validation="CV")
  min.pcr = which.min(MSEP(pcr_model)$val[1,1, ] ) - 1
  #predict
  predict(pcr_model, cd.test, min.pcr)
  #evaluate
  cor.pcr[i]<-cor(cd.test$sleep, predict.cp)
}
boxplot(cor.cp, xlab="PCR", ylab='Test correlation', col = 3)

