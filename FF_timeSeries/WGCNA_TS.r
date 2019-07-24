# The code is exactly the functions in WGCNA_TS_Documentation.ipynb
source("TimeSeriesSimilarity.r") # used for similarity matrix
library("WGCNA")

# My function using time series distance measure
WGCNA_TS = function(X,len_time,softPower=8,minModuleSize=30,MEDissThres=0.25,
                   type = "cor"){
    # similairity matrix (in terms of time series)
    sim_mat = similarity_TS(X,len_time = len_time,type = type) # time series measure
    # adjacency matrix
    adjacency_mat = sim_mat^softPower
    # TOM
    TOM = TOMsimilarity(adjacency_mat)
    # distance matrix
    dissTOM = 1-TOM
    #cluster using dissTom
    geneTree = hclust(as.dist(dissTOM), method = "average")
    # Module identification using dynamic tree cut (cut the genetree)
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                                pamRespectsDendro = FALSE,minClusterSize = minModuleSize)

    # merge similar modules
    if (length(table(dynamicMods))==1){
        names(dynamicMods) = names(X)
        return (dynamicMods)
    }

    # Calculate eigengenes
    # colors is the current module membership
    MEList = moduleEigengenes(X, colors = dynamicMods)  
    MEs = MEList$eigengenes
    # Calculate dissimilarity of module eigengenes
    MEDiss = 1-cor(MEs)
    # Cluster module eigengenes
    METree = hclust(as.dist(MEDiss), method = "average")
    # We choose a height cut of MEDissThres (= 0.25 for example) corresponding to 
    # correlation of 1-MEDissThres (0.75) to merge
    # Call an automatic merging function; dynamicMods is the membership before merge
    merge = mergeCloseModules(X, dynamicMods, cutHeight = MEDissThres, 
                              verbose = 3)
    moduleColors = merge$colors
    names(moduleColors) = names(X)
    return (moduleColors)
}

# Only for test, regular WGCNA using correlation
WGCNA_regular = function(X,len_time,softPower=6,minModuleSize=30,MEDissThres=0.25){
    # similairity matrix 
    sim_mat = abs(cor(X)) # wgcna's naive way
    # adjacency matrix
    adjacency_mat = sim_mat^softPower
    # TOM
    TOM = TOMsimilarity(adjacency_mat)
    # distance matrix
    dissTOM = 1-TOM
    #cluster using dissTom
    geneTree = hclust(as.dist(dissTOM), method = "average")
    # Module identification using dynamic tree cut (cut the genetree)
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                                pamRespectsDendro = FALSE,minClusterSize = minModuleSize)

    # merge similar modules
    if (length(table(dynamicMods))==1){
        names(dynamicMods) = names(X)
        return (dynamicMods)
    }

    # Calculate eigengenes
    # colors is the current module membership
    MEList = moduleEigengenes(X, colors = dynamicMods)  
    MEs = MEList$eigengenes
    # Calculate dissimilarity of module eigengenes
    MEDiss = 1-cor(MEs)
    # Cluster module eigengenes
    METree = hclust(as.dist(MEDiss), method = "average")
    # We choose a height cut of MEDissThres (= 0.25 for example) corresponding to 
    # correlation of 1-MEDissThres (0.75) to merge
    # Call an automatic merging function; dynamicMods is the membership before merge
    merge = mergeCloseModules(X, dynamicMods, cutHeight = MEDissThres, 
                              verbose = 3)
    moduleColors = merge$colors
    names(moduleColors) = names(X)
    return (moduleColors)
}
