#Transfer learning self-starting T2 control chart

library(foreach)
library(parallel)
library(iterators)
library(doParallel)
library(miscTools)
library(maxLik)

options(scipen=999) # Disable scientific notation globally
memory.limit(102400)

Transfer_T2_chart <- function(hs, UCL, gamma, tau, K, d0, d1){
  
  #####Parameter settings#####
  # hs represents the number of historical profiles
  # UCL represents the control limit
  # gamma represents the similarity threshold
  # tau represents the starting position of the shift
  # K is the number of simulation runs
  # d0 represents the multiple of standard deviation for beta0 shift
  # d1 represents the multiple of standard deviation for beta1 shift
  
  theta0 <- sqrt(0.0117) # Standard deviation of beta0
  theta1 <- sqrt(0.0207) # Standard deviation of beta1
  
  beta0 <- 3 # Profile parameter beta0 for GLM
  beta1 <- 2 # Profile parameter beta1 for GLM
  M <- 10 # Number of observations in a profile
  X <- 0.1*c(1:M) # Explanatory variable
  
  lt <- 1000 # Total number of profiles in the target domain
  no_core <- 4 # Number of CPU logical processors
  
  #####Simulation process#####
  t1 <- proc.time() # Record start time
  
  # Start parallel computing
  cl <- makeCluster(no_core)
  registerDoParallel(cl)
  
  #####Functions required for transfer learning#####
  
  # Generate Poisson distributed profile response variables, parameters are (b0, b1)
  # For each profile, generate M Poisson random numbers as response variables
  poi_data <- function(b0, b1, M){
    y <- sapply(c(1:M), function(i) rpois(1, exp(b0+b1*0.1*i)))
    return(y)
  }
  
  # Generate Poisson distributed candidate source domain profile response variable dataset, parameters are (b10,b11) and (b20,b21)
  # Generate K sets of candidate source domain profile response variables, each set with M observations
  poidb_source <- function(b10, b11, b20, b21, K, M){
    D <- matrix(nrow = K, ncol = M)   
    # First type of candidate source domain profile response variables, parameters are (b10, b11), half of K
    for(i in 1:round(0.5*K)){
      D[i,] <- poi_data(b0 = b10, b1 = b11, M = M)
    }
    # Second type of candidate source domain profile response variables, parameters are (b10, b11), half of K
    for(i in round(0.5*K+1):K){
      D[i,] <- poi_data(b0 = b20, b1 = b21, M = M)
    }
    return(D)
  }
  
  # Generate Poisson distributed controlled target domain profile response variable dataset
  # Generate K sets of controlled target domain profile response variables, each set with M observations
  poidb_ic <- function(b0, b1, K, M){
    D <- matrix(nrow = K, ncol = M)   
    for(i in 1:K){
      D[i,] <- poi_data(b0 = b0, b1 = b1, M = M)
    }
    return(D)
  }
  
  # Iterative weighted least squares estimation (IWLS) to estimate beta for Poisson GLM
  # Maximum likelihood estimation (MLE) about beta
  # DY as response variable, DX as explanatory variable
  MLE_beta_GLM <- function(DY, DX){
    MLE <- glm(DY ~ ., data = as.data.frame(DX), family = poisson(link = "log"))
    b <- as.numeric(coef(MLE))
    return(b)
  }
  
  # Fit Poisson GLM regression coefficients using IWLS
  # Y as response variable, X as explanatory variable
  GLM_fit <- function(Y, X){
    b <- MLE_beta_GLM(Y,X)  # Use maximum likelihood estimation (MLE) for beta_0
    return(b)
  }
  
  # Loss function calculation
  loss_function <- function(b, Y){
    l1 <- sum((b[1]+b[2]*X)*Y)
    l2 <- sum(exp(b[1]+b[2]*X))
    l3 <- sum(log(sqrt(2*pi*Y))+Y*log(Y)-Y) 
    l <- l1-l2-l3
    return(l)
  }
  
  # Calculate similarity between source domain and target domain profiles
  # D.t represents target domain response variable, D.s represents source domain response variable
  similarity <- function(D.t, D.s){
    loss.t <- loss_function(GLM_fit(D.t, X), D.t)
    loss.s <- apply(D.s, 1, function(Y) loss_function(trans_fit(D.t = D.t, D.s = Y), D.t))
    zeta <- 1/(1+abs(loss.s-loss.t))
    return(zeta)
  }
  
  # Transfer learning parameter estimation method, using both source domain and target domain
  # Maximum likelihood estimation (MLE) about beta
  # DY as response variable
  MLE_beta_trans <- function(DY){
    N <-  nrow(DY) 
    matrix_X <- t(replicate(N, X))
    vector_DY <- as.vector(DY)
    vector_DX <- as.vector(matrix_X)
    MLE <- glm(vector_DY ~ vector_DX, family = poisson(link = "log"))
    b <- as.numeric(coef(MLE))
    return(b)
  }
  
  # Auxiliary function for summing matrices of different profiles during iteration
  # D as response variable, b as coefficients
  AF <- function(D, b){
    eta <- b[1]+b[2]*X
    mu <- exp(eta)
    W <- diag(mu)
    x <- t(rbind(rep(1,M),X))
    a <- t(x)%*%W%*%x
    return(a)
  }
  
  #Parameter estimation method for transfer learning
  #D.s represents the source domain，D.t represents the target domain
  trans_fit <- function(D.t, D.s){
    D <- rbind(D.t, D.s)
    b <- MLE_beta_trans(D)  # Obtain the initial value for iteration through Maximum Likelihood Estimation (MLE)
    return(b)
  }
  
  # Calculate covariance matrix
  # X as explanatory variable, b as profile parameters
  cov_fit <- function(X, b){
    W <- diag(exp(b[1]+b[2]*X))
    x <- t(rbind(rep(1,M),X))
    C <- solve(t(x)%*%W%*%x)
    return(C)
  }
  
  # Use parallel computing to obtain results of K simulation runs
  Tranf_T2 <- foreach(j = 1:K, .combine = rbind)%dopar%{
    library(maxLik)
    
    # Generate 10 groups of profile response variables with parameters (3.10, 1.80) and (3.15, 1.75) as candidate source domain
    D.s <- poidb_source(b10 = 3.10, b11 = 1.80, b20 = 3.15, b21 = 1.75, K = 20, M = M)
    # Generate hs groups of historical controlled sample profile response variables
    D.T <- poidb_ic(b0 = beta0, b1 = beta1, K = hs, M = M)
    # Use the average of historical controlled sample profile parameters as the initial values for self-starting monitoring process parameters
    B0 <- apply(apply(D.T, 1, function(Y) GLM_fit(Y, X)), 1, mean)
    
    RL <- 0  # Run length
    T2_1 <- vector(length = lt)  # Record T2 statistics
    TA1 <- 0
    # Transfer learning self-starting monitoring process
    for(i in 1:lt){
      
      # Generate target domain profiles
      if(i <= tau){
        # Before the shift occurs, generate controlled profile response variables
        D.t <- poi_data(b0 = beta0, b1 = beta1, M = M)
        
      }else{
        # After the shift occurs, generate out-of-control profile response variables
        U0 <- beta0+d0*theta0
        U1 <- beta1+d1*theta1
        D.t <- poi_data(b0 = U0, b1 = U1, M = M)
      }
      # Source domain selection method based on similarity
      # Calculate the similarity between the target domain and each candidate source domain
      sim <- similarity(D.t, D.s)
      # Filter out candidate source domain profiles with similarity higher than the threshold
      TA <- which(sim >= gamma)
      
      # Profile parameter estimation
      if(length(TA)==0){
        # When no suitable source domain for transfer learning is selected
        # Estimate the parameters of the target domain profile directly through IWLS
        B1 <- GLM_fit(D.t, X)
      }else{
        # When a suitable source domain for transfer learning is selected
        # Select profiles that meet the conditions as the final source domain from the candidate source domain
        D <- D.s[TA,] 
        # Estimate the parameters of the target domain profile through transfer learning
        B1 <- trans_fit(D.t = D.t, D.s = D)
      }
      cov1 <- cov_fit(X, B0)  # Calculate covariance matrix
      
      # Calculate the cumulative number of transfer samples
      TA1 <- TA1 + length(TA)
      
      # Self-starting T2 control chart
      T2_1[i] <- t(B1-B0)%*%solve(cov1)%*%(B1-B0)
      
      # Update controlled process parameters
      B0 <- B0+(B1-B0)/(i+1)
      
      # If the statistic exceeds the control limit, assign the run length and end monitoring
      if(UCL < T2_1[i]){
        RL <- i
        TA2 = TA1%/%i # Calculate the average number of transfer samples before the breakpoint
        break
      }
    }
    
    # If the statistic does not exceed the control limit until the end, record the run length as the total number of profiles in the target domain
    if(RL == 0){
      RL <- lt
      # Calculate the average number of transfer samples
      TA2 = TA1%/%lt
    }
    
    # Return run length as a result of parallel computing
    my_list1<- list(RL, TA2)
    names(my_list1) <- c("RL", "TA2")
    return(my_list1)
  }
  stopCluster(cl)  # Terminate parallel computing
  
  t2 <- proc.time()  # Record end time
  t2-t1  # Total running time
  
  RL <- as.numeric(Tranf_T2[,1])
  TA2 <- as.numeric(Tranf_T2[,2])
  
  filtered_T2<- RL[RL > tau] # Keep ARL greater than the shift time tau
  filtered_T2<- filtered_T2 - tau # ARL from the breakpoint to the alarm
  SSARL1 <- mean(filtered_T2, na.rm = TRUE) # Steady-state SSARL1
  
  SSSD <- sd(filtered_T2, na.rm = TRUE) # Steady-state SSSD
  
  MTA2 <- mean(TA2, na.rm = TRUE) # Mean transfer sample size
  
  my_list<- list(SSARL1, SSSD, MTA2)
  names(my_list) <- c("SSARL1", "SSSD", "Average Transfer Sample Size")
  return(my_list)
}
