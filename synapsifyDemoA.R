# synapsifyDemoA.R

require(sss)
require(Biobase)
require(ggplot2)
require(breastCancerTRANSBIG)
require(ggplot2)
require(ROCR)
require(hgu133a.db)
require(synapseClient)


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

## BINARY MODEL OF 'ER Status' USING SSS
sssERFit <- sss(trainScore ~ t(trainExpress))

# EVALUATE AND VISUALIZE TRAINING Y-HAT
trainScoreHat <- predict(sssERFit, newdata = t(trainExpress))

trainScoreDF <- as.data.frame(cbind(trainScore, trainScoreHat))
colnames(trainScoreDF) <- c("yTrain", "yTrainHat")
ggplot(trainScoreDF, aes(factor(yTrain), yTrainHat)) + geom_boxplot() +
  geom_jitter(aes(colour = as.factor(yTrain)))

# VALIDATE & VISUALIZE WITH HELD OUT VALIDATION COHORT
validScoreHat <- predict(sssERFit, newdata = t(validExpress))
validScoreDF <- as.data.frame(cbind(validScore, validScoreHat))
colnames(validScoreDF) <- c("yValid", "yValidHat")
validBoxPlot <- ggplot(validScoreDF, aes(factor(yValid), yValidHat)) + 
  geom_boxplot() +
  geom_jitter(aes(colour = as.factor(yValid)))

# Alternative visualization (density plots)
validDensPlot <- ggplot(validScoreDF, aes(x = validScoreHat, 
                                          fill = factor(validScore))) + 
  geom_density(alpha = 0.3)

# EVALUATE VALIDATION MODEL PERFORMANCE
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

# Sensitivity at Youden's J-point = 0.91
# Specificity at Youden's J-point = 0.79

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
qplot(erPMP[1:20], geom = "histogram")
topERGenes <- as.character(mget(names(erPMP)[1:25], hgu133aSYMBOL, 
                                 ifnotfound = NA))

## PUSH INTERMEDIATE DATA AND MODEL OBJECTS TO SYNAPSE
synapseLogin("synapsify@sagebase.org", "XXXXXXX")
myProject <- Project(list(name = "Synapse Demonstration"))
myProject <- createEntity(myProject)
# An object of class "Project"
# Synapse Entity Name : Synapse Demonstration
# Synapse Entity Id   : 162999
# Parent Id           : 4489

erSSSModelDS <- Dataset(list(name = "ER SSS Model Construction",
                             parentId = properties(myProject)$id))
erSSSModelDS <- createEntity(erSSSModelDS)
An object of class "Dataset"
# Synapse Entity Name : ER SSS Model Construction
# Synapse Entity Id   : 163000
# Parent Id           : 162999
# Version Number      : 1
# Version Label       : 0.0.0

trainDL <- Layer(list(name = "TRANSBIG Training Cohort", 
                      parentId = properties(erSSSModelDS)$id,
                      type = "E"))
trainDL <- createEntity(trainDL)
trainDL <- addObject(trainDL, trainExpress)
trainDL <- storeEntity(trainDL)
# An object of class "ExpressionLayer"
# Synapse Entity Name : TRANSBIG Training Cohort
# Synapse Entity Id   : 163002
# Parent Id           : 163000
# Type                : E
# Version Number      : 1
# Version Label       : 0.0.0
# 
# loaded object(s):
#   [1] "trainExpress" (matrix)

validDL <- Layer(list(name = "TRANSBIG Validation Cohort", 
                      parentId = properties(erSSSModelDS)$id,
                      type = "E"))
validDL <- createEntity(validDL)
validDL <- addObject(validDL, validExpress)
validDL <- storeEntity(validDL)
# An object of class "ExpressionLayer"
# Synapse Entity Name : TRANSBIG Validation Cohort
# Synapse Entity Id   : 163004
# Parent Id           : 163000
# Type                : E
# Version Number      : 1
# Version Label       : 0.0.0
# 
# loaded object(s):
#   [1] "validExpress" (matrix)

modelDL <- Layer(list(name = "TRANSBIG ER SSS Model Object", 
                      parentId = properties(erSSSModelDS)$id,
                      type = "E"))
modelDL <- createEntity(modelDL)
modelDL <- addObject(modelDL, sssERFit)
modelDL <- storeEntity(modelDL)
# An object of class "ExpressionLayer"
# Synapse Entity Name : TRANSBIG ER SSS Model Object
# Synapse Entity Id   : 163006
# Parent Id           : 163000
# Type                : E
# Version Number      : 1
# Version Label       : 0.0.0
# 
# loaded object(s):
#   [1] "sssERFit" (sssBinaryResult)
