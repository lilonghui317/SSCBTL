# Non-transfer learning traditional self-starting MEWMA control chart

library(maxLik)
library(foreach)
library(parallel)
library(iterators)
library(doParallel)
library(miscTools)

options(scipen=999) # Disable scientific notation globally
memory.limit(102400)

Non_transfer_MEWMA_chart <- function(hs, UCL, lamda, tau, K, d0, d1){
  ##### Parameter Settings #####
  # hs: number of historical profiles
  # UCL: Upper Control Limit
  # tau: shfit start position
  # K: number of simulation runs
  # d0: shift multiplier for beta0 by its standard deviation
  # d1: shift multiplier for beta1 by its standard deviation
  # lamda: weight of MEWMA statistic
  
  theta0 <- sqrt(0.0117) # Standard deviation of beta0
  theta1 <- sqrt(0.0207) # Standard deviation of beta1
  
  beta0 <- 3  # Generalized linear model profile parameter beta0
  beta1 <- 2  # Generalized linear model profile parameter beta1
  M <- 10  # Number of observations in the profile
  X <- 0.1*c(1:M)  # Explanatory variable
  lt <- 1000 # Total number of profiles in the target domain
  no_core <- 4  # Number of CPU logical processors
  
  ##### Functions required for non-transfer learning method #####
  
  # Generate Poisson distribution candidate source domain profile response variable dataset, with parameters (b10, b11) and (b20, b21)
  # Generate K sets of candidate source domain profile response variables, each with M observations
  poidb_source <- function(b10, b11, b20, b21, K, M){
    D <- matrix(nrow = K, ncol = M)   
    # First type of candidate source domain profile response variable, with parameters (b10, b11), half of the quantity of K
    for(i in 1:round(0.5*K)){
      D[i,] <- poi_data(b0 = b10, b1 = b11, M = M)
    }
    # Second type of candidate source domain profile response variable, with parameters (b20, b21), half of the quantity of K
    for(i in round(0.5*K+1):K){
      D[i,] <- poi_data(b0 = b20, b1 = b21, M = M)
    }
    return(D)
  }
  
  # Generate Poisson distribution profile response variable, with parameters (b0, b1)
  # For each profile, generate M Poisson random numbers as response variables
  poi_data <- function(b0, b1, M){
    y <- sapply(c(1:M), function(i) rpois(1, exp(b0+b1*0.1*i)))
    return(y)
  }
  
  # Generate Poisson distribution controlled target domain profile response variable dataset
  # Generate K sets of controlled target domain profile response variables, each with M observations
  poidb_ic <- function(b0, b1, K, M){
    D <- matrix(nrow = K, ncol = M)   
    for(i in 1:K){
      D[i,] <- poi_data(b0 = b0, b1 = b1, M = M)
    }
    return(D)
  }
  
  # Iterative weighted least squares estimation (IWLS), using a single generalized linear model profile to estimate parameters
  # Maximum likelihood estimation of beta
  # DY: response variable, DX: explanatory variable
  MLE_beta_GLM <- function(DY, DX){
    MLE <- glm(DY ~ ., data = as.data.frame(DX), family = poisson(link = "log"))
    b <- as.numeric(coef(MLE))
    return(b)
  }
  
  # Fitting Poisson distribution GLM regression coefficients using IWLS
  # Y: response variable, X: explanatory variable
  GLM_fit <- function(Y, X){
    b <- MLE_beta_GLM(Y,X)  
    return(b)
  }
  
  # Loss function calculation
  # b: parameters, Y: explanatory variable
  loss_function <- function(b, Y){
    l1 <- sum((b[1]+b[2]*X)*Y)
    l2 <- sum(exp(b[1]+b[2]*X))
    l3 <- sum(log(sqrt(2*pi*Y))+Y*log(Y)-Y)   # Sum of factorial logarithms
    l <- l1-l2-l3
    return(l)
  }
  
  # Calculate covariance matrix
  # X: explanatory variable, b: profile parameters
  cov_fit <- function(X, b){
    W <- diag(exp(b[1]+b[2]*X))
    x <- t(rbind(rep(1,M),X))
    C <- solve(t(x)%*%W%*%x)
    return(C)
  }
  
  ##### Simulation Process #####
  t1 <- proc.time()  # Record start time
  
  # Start parallel computing
  cl <- makeCluster(no_core)
  registerDoParallel(cl)
  
  # Use parallel computing to get K simulation results
  MEWMA <- foreach(j = 1:K, .combine = rbind)%dopar%{
    
    library(maxLik)
    # Generate hs initial controlled sample profiles
    D.t <- poidb_ic(b0 = beta0, b1 = beta1, K = hs, M = M)
    # Use the average of the initial values ​​of the controlled sample profile parameters as the initial value of the self-starting monitoring process parameters
    b0 <- apply(apply(D.t, 1, function(Y) GLM_fit(Y, X)), 1, mean)
    
    RL <- 0  # Run Length
    MEWMA1 <- vector(length = lt)  # Record MEWMA statistic
    
    # Self-starting monitoring process
    for(i in 1:lt){
      
      # Generate target domain profile
      if(i <= tau){
        # Before the shift, generate controlled profiles
        D.t <- poi_data(b0 = beta0, b1 = beta1, M = M)
      }else{
        # After the shift, generate out-of-control profiles
        D.t <- poi_data(b0 = beta0+d0*theta0, b1 = beta1+d1*theta1, M = M)
      }
      
      # Estimate generalized linear model profile parameters using IWLS method
      b1 <- GLM_fit(D.t, X)
      cov1 <- cov_fit(X, b0)  # Calculate covariance matrix
      
      # Self-starting MEWMA control chart
      w <- t(b1-b0)%*%solve(cov1)%*%(b1-b0)
      if(i == 1){
        MEWMA1[i] <- lamda*w
      }else{
        MEWMA1[i] <- (1-lamda)*MEWMA1[i-1] + lamda*w
      }
      
      # Update controlled parameters
      b0 <- b0+(b1-b0)/(i+1)
      
      # If the statistic exceeds the control limit, assign the run length and end monitoring
      if(UCL < MEWMA1[i]){
        RL <- i
        break
      }
    }
    
    # If the statistic has not exceeded the control limit until the end, set the run length to the total number of profiles in the target domain
    if(RL == 0){
      RL <- lt
    }
    
    # Return the run length as a value of parallel computing
    return(RL)
  }
  stopCluster(cl)  # Terminate parallel computing
  
  t2 <- proc.time()  # Record end time
  t2-t1  # Total running time
  
  
  filtered_T2<- MEWMA[MEWMA > tau] # Keep ARL greater than or equal to shift time tau
  filtered_T2<- filtered_T2 - tau # Average run length from change point to alarm
  SSARL1 <- mean(filtered_T2, na.rm = TRUE) # Steady-state SSARL1
  SSSD <- sd(filtered_T2, na.rm = TRUE) # Steady-state SSSD
  my_list<- list(SSARL1, SSSD)
  names(my_list) <- c("SSARL1", "SSSD")
  return(my_list)
}
