library("proxy")
library("dtw")
library("rucrdtw")

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
    
     ## type = fastDTW, no need to normalize ##
    if (type == "fastDTW"){
        # pre-allocate memory
        dist_mat = matrix(0,p,p)
        if (n==1){
            for (i in 1:(p-1)){
                for (j in (i+1):p){
                    dist_mat[i,j] = ucrdtw_vv(X[,i], X[,j], 0.05)$distance
                }
            }
        }else{
            flag = 1 # the first time index of the current sample
            for (sample in 1:n){
                for (i in 1:(p-1)){
                    for (j in (i+1):p){
                        dist_mat[i,j] = (dist_mat[i,j]
                        +(ucrdtw_vv(X[flag:(flag+len_time[sample]-1),i], 
                                                  X[flag:(flag+len_time[sample]-1),j],
                                                  0.05)$distance)*len_time[sample])
                    }
                }                
                flag = flag+len_time[sample]
            }
            dist_mat = dist_mat/(sum(len_time))
        }
        # complete the matrix (diag=0 automatically)
        dist_mat[lower.tri(dist_mat)] = t(dist_mat)[lower.tri(dist_mat)]
        
        # transform dist to similarity
        sim_matrix = 1/(1+dist_mat)
        
        return (sim_matrix)
        } # end "fastDTW"
    
    #### data Z score normaliztion #####
    X = normalize_TS(X,len_time) # normalize time series
    #check for missing value (may due to sigma=0)
    tmp = sum(is.na(X))
    if(tmp>0){
        stop(sprintf("There are %d missing values in normalized X",tmp))
    }
    ### end normalization ###
    
    
    # "L2" using dist(M) which computes dist between rows of M
    if (type == "L2"|type == "Euclidean"){
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
    
    # "dtw" using dtwDist(M) which computes dist between rows of M
    if (type == "dtw" | type == "DTW"){
        # first compute distance matrix
        if (n==1){
            dist_obj = dtwDist(t(X))
        }else{
            # dist_obj can be computed like a vector! but I have to create it
            dist_obj = dtwDist(t(X[1:len_time[1],]))*len_time[1]
            flag = 1+len_time[1] # the first time index of the current sample
            for (sample in 2:n){
                dist_obj = (dist_obj + 
            dtwDist(t(X[flag:(flag+len_time[sample]-1),]))*len_time[sample])
                
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
        } # end "dtw"
    
    
    stop("The type is not availble")
    
    }