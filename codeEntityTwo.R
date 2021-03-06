# codeEntityTwo.R

# Code Entity 2 of the Synapsify Demo Project
# Taking the output from Code Entity 1, generating a model for ER status


## CODE ENTITY TWO: FUNCTION A
buildModel <- function(returnOne){
# SOURCE NECESSARY LIBRARIES
require(ggplot2)
require(randomForest)

## BINARY MODEL OF 'ER Status' USING SSS
sssERFit <- randomForest(t(returnOne$trainExpress), 
                         factor(returnOne$trainScore), 
                         ntree = 50,
                         do.trace = 5)

# EVALUATE AND VISUALIZE TRAINING Y-HAT
trainScoreHat <- predict(sssERFit, t(returnOne$trainExpress),
                         type = "prob")

trainScoreDF <- as.data.frame(cbind(returnOne$trainScore, trainScoreHat))
colnames(trainScoreDF) <- c("yTrain", "yTrainHat")
trainBoxPlot <- ggplot(trainScoreDF, aes(factor(yTrain), yTrainHat)) + 
  geom_boxplot() +
  geom_jitter(aes(colour = as.factor(yTrain)), size = 4) +
  opts(title = "ER SSS Model Training Set Hat") +
  ylab("Training Set ER Prediction") +
  xlab("True ER Status") +
  opts(plot.title = theme_text(size = 14))


return(list("sssERFit " = sssERFit,
            "trainBoxPlot" = trainBoxPlot,
            "returnOne" = returnOne)) # Workaround here
}





## CODE ENTITY TWO: FUNCTION B
validateModel <- function(returnTwo){
  
  # SOURCE LIBRARIES
  require(randomForest)
  require(ggplot2)
  require(ROCR)
  
  # DEFINE VARIABLES
  sssERFit <- returnTwo$sssERFit
  validExpress <- returnTwo$returnOne$validExpress
  validScore <- returnTwo$returnOne$validScore
  
  # VALIDATE & VISUALIZE WITH HELD OUT VALIDATION COHORT
  validScoreHat <- predict(sssERFit, t(validExpress), type = "prob")
  validScoreHat <- validScoreHat[ , 2]
  validScoreDF <- as.data.frame(cbind(validScore, validScoreHat))
  colnames(validScoreDF) <- c("yValid", "yValidHat")
  validBoxPlot <- ggplot(validScoreDF, aes(factor(yValid), yValidHat)) + 
    geom_boxplot() +
    geom_jitter(aes(colour = as.factor(yValid)), size = 4) +
    opts(title = "ER SSS Model Indepenedent Validation Set") +
    ylab("Validation Set ER Prediction") +
    xlab("True ER Status") +
    opts(plot.title = theme_text(size = 14))
  
  # Alternative visualization (density plots)
  validDensPlot <- ggplot(validScoreDF, 
                          aes(yValidHat, 
                              fill = factor(yValid))) + 
                                geom_density(alpha = 0.3) +
                                ylab("Density") +
                                xlab("True ER Status") +
                                opts(plot.title = theme_text(size = 14))
  
  # EVALUATE VALIDATION MODEL PERFORMANCE
  erPred <- prediction(as.numeric(validScoreHat), as.numeric(validScore))
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
  
  ## CREATE A ROC CURVE USING GGPLOT
  dfPerf <- as.data.frame(cbind(unlist(erPerf@x.values),
                                unlist(erPerf@y.values)))
  colnames(dfPerf) <- c("FalsePositiveRate", "TruePositiveRate")
  
  rocCurve <- ggplot(dfPerf, aes(FalsePositiveRate, TruePositiveRate)) +
    geom_line() + 
    geom_abline(slope = 1, colour = "red") +
    opts(title = "Validation Cohort ROC Curve") +
    ylab("False Positive Rate") +
    xlab("True Positive Rate") +
    opts(plot.title = theme_text(size = 14))
  
  ## RETURN
  return(list("validScoreDF" = validScoreDF,
         "validBoxPlot" = validBoxPlot,
         "validDensPlot" = validDensPlot,
         "rankSum" = rankSum,
         "rocCurve" = rocCurve,
         "sensitivity" = optSens,
         "specificity" = optSpec,
         "auc" = erAUC@y.values))
  
}