# transbigRandomCohorts.R

# Code Entity 1 of the Synapsify Demo Project
# Taking the full TRANSBIG breast cancer dataset and randomly dividing it into
# "training" and "validation" cohorts and concomitantly dividing the clinical
# "ER Status" phenotype along with the cohorts

# SOURCE IN NECESSARY LIBRARIES
require(synapseClient)

# LOGIN TO SYNAPSE
synapseLogin("public@sagebase.org", "public")

# LOAD THE SYNAPSE ENTITY INTO MEMORY
transbig <- loadEntity(169192)
expressData <- exprs(transbig)
pheno <- phenoData(transbig)

## CREATE TRAINING AND VALIDATION COHORTS
set.seed(031512)
randVec <- rbinom(dim(transbig)[2], size = 1, prob = 0.5)
trainExpress <- expressData[ , randVec == 0]
validExpress <- expressData[ , randVec == 1]
trainScore <- as.numeric(pheno@data$er[randVec == 0])
validScore <- as.numeric(pheno@data$er[randVec == 1])

# END CODE ENTITY 1

