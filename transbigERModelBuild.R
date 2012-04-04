# transbigERModelBuild.R

# Code Entity 2 of the Synapsify Demo Project
# Taking the output from Code Entity 1, generating a model for ER status

transbigERModelBuild <- function(returnValue){
# SOURCE NECESSARY LIBRARIES
require(ggplot2)
require(sss)

## BINARY MODEL OF 'ER Status' USING SSS
sssERFit <- sss(returnValue$trainScore ~ t(returnValue$trainExpress))

# EVALUATE AND VISUALIZE TRAINING Y-HAT
trainScoreHat <- predict(sssERFit, newdata = t(returnValue$trainExpress))

trainScoreDF <- as.data.frame(cbind(returnValue$trainScore, trainScoreHat))
colnames(trainScoreDF) <- c("yTrain", "yTrainHat")
trainBoxPlot <- ggplot(trainScoreDF, aes(factor(yTrain), yTrainHat)) + 
  geom_boxplot() +
  geom_jitter(aes(colour = as.factor(yTrain)), size = 4) +
  opts(title = "ER SSS Model Training Set Hat") +
  ylab("Training Set ER Prediction") +
  xlab("True ER Status") +
  opts(plot.title = theme_text(size = 14))


return(list("sssERFit " = sssERFit,
            "trainBoxPlot" = trainBoxPlot))
}