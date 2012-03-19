# ....###....##....##....###....##.......##....##..######..####..######.
# ...##.##...###...##...##.##...##........##..##..##....##..##..##....##
# ..##...##..####..##..##...##..##.........####...##........##..##......
# .##.....##.##.##.##.##.....##.##..........##.....######...##...######.
# .#########.##..####.#########.##..........##..........##..##........##
# .##.....##.##...###.##.....##.##..........##....##....##..##..##....##
# .##.....##.##....##.##.....##.########....##.....######..####..######.
# 
# .########..##........#######...######..##....##....##.......
# .##.....##.##.......##.....##.##....##.##...##.....##....##.
# .##.....##.##.......##.....##.##.......##..##......##....##.
# .########..##.......##.....##.##.......#####.......##....##.
# .##.....##.##.......##.....##.##.......##..##......#########
# .##.....##.##.......##.....##.##....##.##...##...........##.
# .########..########..#######...######..##....##..........##.

require(sss)
require(Biobase)
require(ggplot2)
require(breastCancerTRANSBIG)
require(ggplot2)
require(ROCR)
require(hgu133a.db)
require(synapseClient)

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

rocCurve <- ggplot(dfPerf, aes(FalsePositiveRate, TruePositiveRate)) +
  geom_line() + 
  geom_abline(slope = 1, colour = "red") +
  opts(title = "Validation Cohort ROC Curve") +
  ylab("False Positive Rate") +
  xlab("True Positive Rate") +
  opts(plot.title = theme_text(size = 14))

png(file = "rocCurve.png", bg = "transparent", width = 768, 
    height = 768)
rocCurve
dev.off()

# AUROC = 0.88