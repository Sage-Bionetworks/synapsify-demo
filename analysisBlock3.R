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
# .##.....##.##.......##.....##.##.......##..##.............##
# .##.....##.##.......##.....##.##....##.##...##.....##.....##
# .########..########..#######...######..##....##.....#######.

require(sss)
require(Biobase)
require(ggplot2)
require(breastCancerTRANSBIG)
require(ggplot2)
require(ROCR)
require(hgu133a.db)
require(synapseClient)


# VALIDATE & VISUALIZE WITH HELD OUT VALIDATION COHORT
validScoreHat <- predict(sssERFit, newdata = t(validExpress))
validScoreDF <- as.data.frame(cbind(validScore, validScoreHat))
colnames(validScoreDF) <- c("yValid", "yValidHat")
validBoxPlot <- ggplot(validScoreDF, aes(factor(yValid), yValidHat)) + 
  geom_boxplot() +
  geom_jitter(aes(colour = as.factor(yValid)), size = 4) +
  opts(title = "ER SSS Model Indepenedent Validation Set") +
  ylab("Validation Set ER Prediction") +
  xlab("True ER Status") +
  opts(plot.title = theme_text(size = 14))

png(file = "validBoxPlot.png", bg = "transparent", width = 1024, 
    height = 768)
validBoxPlot
dev.off()

# Alternative visualization (density plots)
validDensPlot <- ggplot(validScoreDF, aes(x = validScoreHat, 
                                          fill = factor(validScore))) + 
                                            geom_density(alpha = 0.3) +
                                            opts(title = "ER SSS Model Independent Validation Set") +
                                            ylab("Density") +
                                            xlab("True ER Status") +
                                            opts(plot.title = theme_text(size = 14))

png(file = "validDensPlot.png", bg = "transparent", width = 1024, 
    height = 768)
validDensPlot
dev.off()