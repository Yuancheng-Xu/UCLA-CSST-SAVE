
library(MASS)

# Input: a matrix x_data; Output, a vector with the f_sim applied to each row of x_data
f_sim = function (X_data){
    y = (5*X_data[,1]+2*X_data[,2]+2*X_data[,3]+5*X_data[,2]*X_data[,3]
         +5*X_data[,301]+2*X_data[,302]+2*X_data[,303]+5*X_data[,302]*X_data[,303])
    return (y)
}

## for all of the following return, the last column is the label
# No time structure; The features are grouped
sim_1 = function (n,T,cor_feature=0.8,var_noise=1){
    p = 400
    p0 = 100
    data = matrix(0,nrow = n*T, ncol = p+1)
    
    #### covariance matrix between features: it is either 0 (independent) or cor_feature ####
    cov_feature = matrix(0,nrow = p, ncol = p)
    # cov within the first three modules
    cov_star = matrix(cor_feature,nrow = p0,ncol = p0)
    diag(cov_star)=1
    # put cov_star into cov_feature
    cov_feature[1:p0,1:p0] = cov_star
    cov_feature[(p0+1):(2*p0),(p0+1):(2*p0)] = cov_star
    cov_feature[(2*p0+1):(3*p0),(2*p0+1):(3*p0)] = cov_star
    cov_feature[(3*p0+1):(4*p0),(3*p0+1):(4*p0)] = diag(p0)
    ####
    
    # Create X matrix
    data[1:(n*T),1:p] = mvrnorm(n=n*T,rep(0,p),cov_feature)
    
    # create label y
    data[1:(n*T),p+1] = (f_sim(data[1:(n*T),1:p]) 
                         + mvrnorm(n=1,rep(0,n*T),diag(x=var_noise,n*T)))
    return (data)
}

# CS time structured on y and grouped features
sim_2 = function(n,T,cor_feature=0.8,var_noise=1,cor_noise=0.8){
    p = 400
    p0 = 100
    data = matrix(0,nrow = n*T, ncol = p+1)
    
    #### covariance matrix between features: it is either 0 (independent) or cor_feature ####
    cov_feature = matrix(0,nrow = p, ncol = p)
    # cov within the first three modules
    cov_star = matrix(cor_feature,nrow = p0,ncol = p0)
    diag(cov_star)=1
    # put cov_star into cov_feature
    cov_feature[1:p0,1:p0] = cov_star
    cov_feature[(p0+1):(2*p0),(p0+1):(2*p0)] = cov_star
    cov_feature[(2*p0+1):(3*p0),(2*p0+1):(3*p0)] = cov_star
    cov_feature[(3*p0+1):(4*p0),(3*p0+1):(4*p0)] = diag(p0)
    ####
    
    #### covariance matrix of noise #####
    # the covariance matrix of noise within a unit
    cov_noise_star = matrix(cor_noise,nrow = T, ncol = T)
    diag(cov_noise_star) = var_noise
    # the overal cov matrix has diagonal block matrix as cov_noise_star
    cov_noise = matrix(0,nrow = n*T, ncol = n*T)
    for (i in 1:n){
        cov_noise[(1+(i-1)*T):(i*T),(1+(i-1)*T):(i*T)] = cov_noise_star
    }
    ###

    # Create X matrix
    data[1:(n*T),1:p] = mvrnorm(n=n*T,rep(0,p),cov_feature)
    
    # create label y
    data[1:(n*T),p+1] = (f_sim(data[1:(n*T),1:p]) 
                         + mvrnorm(n=1,rep(0,n*T),cov_noise))
    
    return(data)
    
}

# Following pdf file: simulation_AR
# x(i+1) = alpha*x(i) + (1-alpha^2)^0.5*std_normal; in the pdf, alpha = 0.8
sim_3 = function(n,T,cor_feature=0.8,var_noise=1,alpha=0.8){
    p = 400
    p0 = 100
    data = matrix(0,nrow = n*T, ncol = p+1)
    
    #### covariance matrix between features: it is either 0 (independent) or cor_feature ####
    cov_feature = matrix(0,nrow = p, ncol = p)
    # cov within the first three modules
    cov_star = matrix(cor_feature,nrow = p0,ncol = p0)
    diag(cov_star)=1
    # put cov_star into cov_feature
    cov_feature[1:p0,1:p0] = cov_star
    cov_feature[(p0+1):(2*p0),(p0+1):(2*p0)] = cov_star
    cov_feature[(2*p0+1):(3*p0),(2*p0+1):(3*p0)] = cov_star
    cov_feature[(3*p0+1):(4*p0),(3*p0+1):(4*p0)] = diag(p0)
    ####
    
    # create x matrix
    tmp = (1-alpha**2)**0.5
    for (i in 1:n){
        data[1+(i-1)*T,1:p] = mvrnorm(n = 1, rep(0, p), cov_feature)
        for (j in 2:T){
            data[j+(i-1)*T,1:p] = (alpha*data[j-1+(i-1)*T,1:p]+
                                    tmp*mvrnorm(n = 1, rep(0, p), cov_feature))
        }
    }
    
    # create y 
    data[1:(n*T),p+1] = ( f_sim(data[1:(n*T),1:p])+ 
                          mvrnorm(n = 1, rep(0,n*T), diag(x=var_noise,n*T)) )
    
    return (data)
}

# CS structure on X
sim_4= function(n,T,cor_feature=0.8,var_noise=1,beta=0.8){
    p = 400
    p0 = 100
    data = matrix(0,nrow = n*T, ncol = p+1)
    
    #### covariance matrix between features: it is either 0 (independent) or cor_feature ####
    cov_feature = matrix(0,nrow = p, ncol = p)
    # cov within the first three modules
    cov_star = matrix(cor_feature,nrow = p0,ncol = p0)
    diag(cov_star)=1
    # put cov_star into cov_feature
    cov_feature[1:p0,1:p0] = cov_star
    cov_feature[(p0+1):(2*p0),(p0+1):(2*p0)] = cov_star
    cov_feature[(2*p0+1):(3*p0),(2*p0+1):(3*p0)] = cov_star
    cov_feature[(3*p0+1):(4*p0),(3*p0+1):(4*p0)] = diag(p0)
    ####
    
    # create X matrix
    tmp = (1-beta**2)**0.5
    for (i in 1:n){
        # note that I put the coefficient beta and tmp in gamma and e directly
        gamma = mvrnorm(n = 1, rep(0, p), cov_feature)*beta
        e = mvrnorm(n = T, rep(0, p), cov_feature)*tmp
        
        for (j in 1:T){
            data[j+(i-1)*T,1:p] = gamma + e[j,] 
        }
    }
    
    # create y 
    data[1:(n*T),p+1] = ( f_sim(data[1:(n*T),1:p])+ 
                          mvrnorm(n = 1, rep(0,n*T), diag(x=var_noise,n*T)) )
    
    return (data)
}
