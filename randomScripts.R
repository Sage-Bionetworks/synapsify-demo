# randomScripts.R

# Save TRANSBIG as a Synapse Entity
require(breastCancerTRANSBIG)
require(synapseClient)
data(transbig)

# synProjEnt <- getEntity(162999)
# transbDataset <- Dataset(list(name = "TRANSBIG Data",
#                               parentId = properties(synProjEnt)$id))
# transbDataset <- createEntity(transbDataset)
# 
# transbLayer <- Layer(list(name = "Full TRANSBIG Data",
#                           parentId = properties(transbDataset)$id,
#                           type = "E"))
# transbLayer <- createEntity(transbLayer)
# transbLayer <- addObject(transbLayer, transbig)
# transbLayer <- storeEntity(transbLayer)

transbDataset <- getEntity(163016)
transbLayer <- Layer(list(name = "TRANSBIG Dataset",
                          parentId = properties(transbDataset)$id,
                          type = "E"))
transbLayer <- createEntity(transbLayer)
transbLayer <- addObject(transbLayer, transbig)
transbLayer <- storeEntity(transbLayer)


# SAVE CODE ENTITY ONE "TRANSBIGRANDOMCOHORTS"
codeEntityOneLayer <- Code(list(name = "CodeEntityOne",
                                parentId = properties(synProjEnt)$id))
codeEntityOneLayer <- createEntity(codeEntityOneLayer)
codeEntityOneLayer <- addFile(codeEntityOneLayer, "transbigRandomCohorts.R")
codeEntityOneLayer <- storeEntity(codeEntityOneLayer)