
library("TSdist")
library("proxy")

# Normalize each time series (not a column)
normalize_TS = function(X,len_time){
    normalized_X = matrix(0,dim(X)[1],dim(X)[2])
    flag = 1 # record the first index of the current time series
    for (t in len_time){
        normalized_X[flag:(flag+t-1),] = scale(X[flag:(flag+t-1),])
        flag = flag + t
    }
    return (normalized_X)
}

# The main function                                                 
# compute similarity using different measure
similarity_TS = function(X,len_time,type){
    p = dim(X)[2] # num of features
    n = length(len_time) # number of samples
    
    ## type = cor, no need to normalize ##
    if (type == "cor"){
        # first compute distance matrix
        if (n==1){
            cor_mat = cor(X)
        }else{
            # cor_mat can be computed like a vector! but I have to create it
            cor_mat = cor(X[1:len_time[1],])*len_time[1]
            flag = 1+len_time[1] # the first time index of the current sample
            for (sample in 2:n){
                cor_mat = (cor_mat + 
            cor(X[flag:(flag+len_time[sample]-1),])*len_time[sample])
                
                flag = flag+len_time[sample]
            }
            cor_mat = cor_mat/(sum(len_time))
        }     
           
        return (cor_mat)
        } # end "cor"

    X = normalize_TS(X,len_time) # normalize time series
    #check for missing value (may due to sigma=0)
    tmp = sum(is.na(X))
    if(tmp>0){
        stop(sprintf("There are %d missing values in normalized X",tmp))
    }
    
    # compute similarity pairwisely with different measure
    
    # "L2" using dist(M) which computes dist between rows of M
    if (type == "L2"){
        # first compute distance matrix
        if (n==1){
            dist_obj = dist(t(X),method = "euclidean")
        }else{
            # dist_obj can be computed like a vector! but I have to create it
            dist_obj = dist(t(X[1:len_time[1],]))*len_time[1]
            flag = 1+len_time[1] # the first time index of the current sample
            for (sample in 2:n){
                dist_obj = (dist_obj + 
            dist(t(X[flag:(flag+len_time[sample]-1),]))*len_time[sample])
                
                flag = flag+len_time[sample]
            }
            dist_obj = dist_obj/(sum(len_time))
        }     
    
        # transform dist to similarity
        sim_obj = 1/(1+dist_obj)
        # similarity matrix
        sim_matrix = as.matrix(sim_obj)
        diag(sim_matrix) = 1
        
        return (sim_matrix)
        } # end "L2"
    
    ## type = cor ##
    if (type == "cor"){
        # first compute distance matrix
        if (n==1){
            cor_mat = cor(X)
        }else{
            # cor_mat can be computed like a vector! but I have to create it
            cor_mat = cor(X[1:len_time[1],])*len_time[1]
            flag = 1+len_time[1] # the first time index of the current sample
            for (sample in 2:n){
                cor_mat = (cor_mat + 
            cor(X[flag:(flag+len_time[sample]-1),])*len_time[sample])
                
                flag = flag+len_time[sample]
            }
            cor_mat = cor_mat/(sum(len_time))
        }     
           
        return (cor_mat)
        } # end "cor"
    
    # Other method: use proxy(dist function) and TSdist package
    # Note: They are really slow!
    if (n==1){
            dist_obj = dist(t(X),method="tsDistances",distance=type)
        }else{
            # dist_obj can be computed like a vector! but I have to create it            
            dist_obj = dist(t(X[1:len_time[1],]),
                            method="tsDistances",distance=type)*len_time[1]
            flag = 1+len_time[1] # the first time index of the current sample
            for (sample in 2:n){
                dist_obj = (dist_obj + 
            dist(t(X[flag:(flag+len_time[sample]-1),]),
                   method="tsDistances",distance=type)*len_time[sample])
                
                flag = flag+len_time[sample]
            }
            dist_obj = dist_obj/(sum(len_time))
        }
    
    # transform dist to similarity
    sim_obj = 1/(1+dist_obj)
    # similarity matrix
    sim_matrix = as.matrix(sim_obj)
    diag(sim_matrix) = 1
        
    return (sim_matrix)
    
    
    }