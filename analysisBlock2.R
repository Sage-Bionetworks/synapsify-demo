# ....###....##....##....###....##.......##....##..######..####..######.
# ...##.##...###...##...##.##...##........##..##..##....##..##..##....##
# ..##...##..####..##..##...##..##.........####...##........##..##......
# .##.....##.##.##.##.##.....##.##..........##.....######...##...######.
# .#########.##..####.#########.##..........##..........##..##........##
# .##.....##.##...###.##.....##.##..........##....##....##..##..##....##
# .##.....##.##....##.##.....##.########....##.....######..####..######.
# 
# .########..##........#######...######..##....##.....#######.
# .##.....##.##.......##.....##.##....##.##...##.....##.....##
# .##.....##.##.......##.....##.##.......##..##.............##
# .########..##.......##.....##.##.......#####........#######.
# .##.....##.##.......##.....##.##.......##..##......##.......
# .##.....##.##.......##.....##.##....##.##...##.....##.......
# .########..########..#######...######..##....##....#########

require(sss)
require(Biobase)
require(ggplot2)
require(breastCancerTRANSBIG)
require(ggplot2)
require(ROCR)
require(hgu133a.db)
require(synapseClient)

## BINARY MODEL OF 'ER Status' USING SSS
sssERFit <- sss(trainScore ~ t(trainExpress))

# EVALUATE AND VISUALIZE TRAINING Y-HAT
trainScoreHat <- predict(sssERFit, newdata = t(trainExpress))

trainScoreDF <- as.data.frame(cbind(trainScore, trainScoreHat))
colnames(trainScoreDF) <- c("yTrain", "yTrainHat")
trainBoxPlot <- ggplot(trainScoreDF, aes(factor(yTrain), yTrainHat)) + 
  geom_boxplot() +
  geom_jitter(aes(colour = as.factor(yTrain)), size = 4) +
  opts(title = "ER SSS Model Training Set Hat") +
  ylab("Training Set ER Prediction") +
  xlab("True ER Status") +
  opts(plot.title = theme_text(size = 14))

png(file = "trainBoxPlot.png", bg = "transparent", width = 1024, 
    height = 768)
trainBoxPlot
dev.off()