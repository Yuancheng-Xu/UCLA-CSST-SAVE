library(MASS)

f_sim = function (n,T,X_data,imp_mod){
    y = rep(0,n*T)
    for (mod in imp_mod){
        a1 = (mod-1)*100+1
        a2 = a1 + 1
        a3 = a2 + 1
        y = y+5*X_data[,a1]+2*X_data[,a2]+2*X_data[,a3]+5*X_data[,a2]*X_data[,a3]
    }     
    return (y)
}

sim_time=function(n,T=5,p=400,imp_mod,cor_feature=0.8,var_re=1,
                  var_noise=1,a1=5,a2=-5){
    p0 = 100
    p_mult = p/100 # number of modules
    if (p_mult%%1 !=0){
        stop("p should be a multiple of 100")
    }
    if (p_mult==1){
        cov_feature = diag(p0) # just one independent group
    }else{
        # Now p_mult>1
        # covariance matrix beween features
        cov_feature = matrix(0,nrow = p, ncol = p)
        # cov of correlated modules
        cov_star = matrix(cor_feature,nrow = p0,ncol = p0)
        diag(cov_star)=1
        # all but the last modules are correlated
        for (k in 1:(p_mult-1)){
            cov_feature[((k-1)*p0+1):(k*p0),((k-1)*p0+1):(k*p0)] = cov_star
        }
        # last modules are independent
        k = p_mult
        cov_feature[((k-1)*p0+1):(k*p0),((k-1)*p0+1):(k*p0)] = diag(p0)
        
    }
     # Create X matrix
    data = mvrnorm(n=n*T,rep(0,p),cov_feature) # observations of X are iid
    data <- data.frame(data)
    names(data) = paste("V",1:p,sep="")

    #### random intercept for each patient ####
    # random intercept draw from N(0,1)
    b = mvrnorm(n = 1, rep(0,n), diag(x=var_re,n))
    data$rand_int = rep(b,each = T)
    ### end random intercept

    data$time <- rep(1:T, n) # time
    # treatment 1 or 2 ,categorical type
    data$treatment[1:(n*T/2)] <- 1 
    data$treatment[((n*T/2)+1):(n*T)] <- 2
    data$treatment = factor(data$treatment)

    # patient information
    data$patient = rep(1:n,each = T)
    
    # noise
    noise = mvrnorm(n = 1, rep(0,n*T), diag(x=var_noise,n*T))

    # response y
    med = median(1:T)
    data$y = (f_sim(n=n,T=T,X_data=data[1:p],imp_mod=imp_mod)+ 
        (data$treatment==1)*a1*(data$time-med)^2 + 
        (data$treatment==2)*a2*(data$time-med)^2 + data$rand_int+noise)
    
    return(data)
}

