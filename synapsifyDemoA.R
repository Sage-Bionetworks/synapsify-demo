# synapsifyDemoA.R

require(sss)
require(Biobase)
require(ggplot2)
require(breastCancerTRANSBIG)
require(ggplot2)
require(ROCR)

data(transbig)
expressData <- exprs(transbig)
pheno <- phenoData(transbig)

## CREATE TRAINING AND VALIDATION COHORTS
set.seed(031512)
randVec <- rbinom(dim(transbig)[2], size = 1, prob = 0.5)
trainExpress <- expressData[ , randVec == 0]
validExpress <- expressData[ , randVec == 1]
trainScore <- as.numeric(pheno@data$er[randVec == 0])
validScore <- as.numeric(pheno@data$er[randVec == 1])

## BINARY MODEL OF 'SCORE' USING SSS
sssERFit <- sss(trainScore ~ t(trainExpress))

# evaluate y_hat
trainScoreHat <- predict(sssERFit, newdata = t(trainExpress))

trainScoreDF <- as.data.frame(cbind(trainScore, trainScoreHat))
colnames(trainScoreDF) <- c("yTrain", "yTrainHat")
ggplot(trainScoreDF, aes(factor(yTrain), yTrainHat)) + geom_boxplot() +
  geom_jitter(aes(colour = as.factor(yTrain)))

# VALIDATE WITH HELD OUT COHORT
validScoreHat <- predict(sssERFit, newdata = t(validExpress))
validScoreDF <- as.data.frame(cbind(validScore, validScoreHat))
colnames(validScoreDF) <- c("yValid", "yValidHat")
ggplot(validScoreDF, aes(factor(yValid), yValidHat)) + geom_boxplot() +
  geom_jitter(aes(colour = as.factor(yValid)))

# EVALUATE MODEL PERFORMANCE
erPred <- prediction(validScoreHat, validScore)
erPerf <- performance(erPred, "tpr", "fpr")
erAUC <- performance(erPred, "auc")

# FIND YOUDEN'S J POINT AND OPTIMAL SENSITIVITY AND SPECIFICITY
erSSPerf <- performance(erPred, "sens", "spec")
youdensJ <- erSSPerf@x.values[[1]] + erSSPerf@y.values[[1]] - 1
jMax <- which.max(youdensJ)
optCut <- erPerf@alpha.values[[1]][jMax]

optSens <- unlist(erSSPerf@x.values)[jMax]
optSpec <- unlist(erSSPerf@y.values)[jMax]

rankSum <- wilcox.test(validScoreHat[validScore == 0], 
                       validScoreHat[validScore == 1])

# Wilcoxon rank sum test with continuity correction
# 
# data:  validScoreHat[validScore == 0] and validScoreHat[validScore == 1] 
# W = 315, p-value = 1.257e-10
# alternative hypothesis: true location shift is not equal to 0

## CREATE A ROC CURVE USING GGPLOT
dfPerf <- as.data.frame(cbind(unlist(erPerf@x.values),
                              unlist(erPerf@y.values)))
colnames(dfPerf) <- c("FalsePositiveRate", "TruePositiveRate")

ggplot(dfPerf, aes(FalsePositiveRate, TruePositiveRate)) +
  geom_line() + geom_abline(slope = 1, colour = "red")

# AUROC = 0.88

## LET'S INSPECT HIGHLY WEIGHTED FEATURES
erPMP <- sssERFit@postMargProb
pmpHist <- qplot(erPMP, geom = "histogram")

