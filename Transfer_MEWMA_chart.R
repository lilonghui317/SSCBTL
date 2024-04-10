#Transfer learning self-starting MEWMA control chart
library(maxLik)
library(foreach)
library(parallel)
library(iterators)
library(doParallel)
library(miscTools)

options(scipen=999) # Disable scientific notation globally
memory.limit(102400)

Transfer_MEWMA_chart <- function(hs, UCL, lamda, gamma, tau, K, d0, d1){
  
  ##### Parameter Settings #####
  # hs: Number of historical profiles
  # UCL: Upper Control Limit
  # gamma: Similarity threshold
  # tau: shift start position
  # K: Number of simulation runs
  # d0: shift multiplier for beta0 by its standard deviation
  # d1: shift multiplier for beta1 by its standard deviation
  # lamda: Weight of the MEWMA statistic
  
  theta0 <- sqrt(0.0117) # Standard deviation of beta0 
  theta1 <- sqrt(0.0207) # Standard deviation of beta1 
  
  beta0 <- 3       # Generalized linear model profile parameter beta0
  beta1 <- 2       # Generalized linear model profile parameter beta1
  M <- 10          # Number of observations in the profile
  X <- 0.1*c(1:M)  # Explanatory variable
  
  lt <- 1000    # Total number of profiles in the target domain
  no_core <- 4  # Number of CPU logical processors
  
  ##### Simulation Process #####
  t1 <- proc.time()  # Record start time
  
  ##### Functions required for transfer learning method #####
  
  # Generate Poisson distribution profile response variable, with parameters (b0, b1)
  # For each profile, generate M Poisson random numbers as response variables
  poi_data <- function(b0, b1, M){
    y <- sapply(c(1:M), function(i) rpois(1, exp(b0+b1*0.1*i)))
    return(y)
  }
  
  # Generate Poisson distribution candidate source domain profile response variable dataset, with parameters (b10,b11) and (b20,b21)
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
    b <- coef(MLE)
    return(b)
  }
  
  # Fitting Poisson distribution GLM regression coefficients using IWLS
  # Y: response variable, X: explanatory variable
  GLM_fit <- function(Y, X){
    b <- MLE_beta_GLM(Y,X)  # Estimate beta_0 using maximum likelihood estimation (MLE)
    return(b)
  }
  
  # Loss function calculation
  loss_function <- function(b, Y){
    l1 <- sum((b[1]+b[2]*X)*Y)
    l2 <- sum(exp(b[1]+b[2]*X))
    # Factorial calculation with factorial() function in R cannot handle factorial calculation for values larger than 170
    # Since Y can have values larger than 170, approximate factorial calculation is used instead
    l3 <- sum(log(sqrt(2*pi*Y))+Y*log(Y)-Y) 
    l <- l1-l2-l3
    return(l)
  }
  
  # Calculate similarity between source and target domain profiles
  # D.t: Target domain profile, D.s: Source domain profile
  similarity <- function(D.t, D.s){
    # Loss function for target domain profile
    loss.t <- loss_function(GLM_fit(D.t, X), D.t)
    # Loss function for source domain profiles
    loss.s <- apply(D.s, 1, function(Y) loss_function(trans_fit(D.t = D.t, D.s = Y), D.t))
    # Quantify similarity based on distance between source and target domain profiles
    zeta <- 1/(1+abs(loss.s-loss.t))
    return(zeta)
  }
  
  # Transfer learning parameter estimation method, jointly using source and target domain
  # Maximum likelihood estimation of beta
  # DY: response variable
  MLE_beta_trans <- function(DY){
    N <-  nrow(DY)  # Set the number of columns in the matrix
    matrix_X <- t(replicate(N, X))
    vector_DY <- as.vector(DY)
    vector_DX <- as.vector(matrix_X)
    MLE <- glm(vector_DY ~ vector_DX, family = poisson(link = "log"))
    b <- coef(MLE)
    return(b)
  }
  
  # Auxiliary function for summing matrices of different profiles during iteration
  # D: response variable, b: coefficients
  AF <- function(D, b){
    eta <- b[1]+b[2]*X
    mu <- exp(eta)
    W <- diag(mu)
    x <- t(rbind(rep(1,M),X))
    a <- t(x)%*%W%*%x
    return(a)
  }
  
  CF <- function(D, b){
    eta <- b[1]+b[2]*X
    mu <- exp(eta)
    W <- diag(mu)
    q <- eta+solve(W)%*%(D-mu)
    x <- t(rbind(rep(1,M),X))
    c <- t(x)%*%W%*%q
    return(c)
  }
  
  # Transfer learning parameter estimation method
  # D.s: Source domain, D.t: Target domain
  trans_fit <- function(D.t, D.s){
    D <- rbind(D.t, D.s)
    b <- MLE_beta_trans(D)  # Get initial values through maximum likelihood estimation (MLE)
    e <- 1
    # Combine target and source domain data for iterative regression coefficient estimation
    while(e > 0.00001){
      A <- matrix(rowSums(apply(D, 1, function(Y) AF(Y, b))), 2, 2)
      C <- rowSums(apply(D, 1, function(Y) CF(Y, b)))
      B <- solve(A)%*%C
      e <- norm(B-b, type="2")/norm(b, type="2") # Euclidean norm
      b <- B
    }
    return(b)
  }
  
  # Calculate covariance matrix
  # X: explanatory variable, b: profile parameters
  cov_fit <- function(X, b){
    W <- diag(exp(b[1]+b[2]*X))
    x <- t(rbind(rep(1,M),X))
    C <- solve(t(x)%*%W%*%x)
    return(C)
  }
  
  # Start parallel computing
  cl <- makeCluster(no_core)
  registerDoParallel(cl)
  
  # Use parallel computing to get the results of K simulation runs
  MEWMA <- foreach(j = 1:K, .combine = rbind) %dopar% {
    
    # Generate 25 sets of profile response variables with parameters (3.10, 1.80) and (3.15, 1.75) as candidate source domains
    D.s <- poidb_source(b10 = 3.10, b11 = 1.80, b20 = 3.15, b21 = 1.75, K = 20, M = M)
    # Generate hs sets of historical controlled sample profile response variables
    D.t <- poidb_ic(b0 = beta0, b1 = beta1, K = hs, M = M)
    # Use the average of historical controlled sample profile parameters as the initial value of self-starting monitoring process parameters
    b0 <- apply(apply(D.t, 1, function(Y) GLM_fit(Y, X)), 1, mean)
    
    RL <- 0  # Run length
    MEWMA1 <- vector(length = lt)  # Record MEWMA statistic
    TA1 <- 0
    # Transfer learning self-starting monitoring process
    for(i in 1:lt){
      
      # Generate target domain profiles
      if(i <= tau){
        # Before the shift occurs, generate controlled profile response variables
        D.t <- poi_data(b0 = beta0, b1 = beta1, M = M)
      }else{
        # After the shift occurs, generate out-of-control profile response variables
        D.t <- poi_data(b0 = beta0+d0*theta0, b1 = beta1+d1*theta1, M = M)
      }
      
      # Source domain selection method based on similarity
      # Calculate the similarity between the target domain and each candidate source domain
      sim <- similarity(D.t, D.s)
      # Filter out candidate source domain profiles with similarity higher than the threshold
      TA <- which(sim >= gamma)
      
      # Profile parameter estimation
      if(length(TA)==0){
        # When no source domain suitable for transfer learning is filtered
        # Estimate the parameters of the target domain profile directly using IWLS method
        b1 <- GLM_fit(D.t, X)
      }else{
        # When source domains suitable for transfer learning are filtered
        # Select the profile that meets the conditions from the candidate source domains as the final source domain
        D <- D.s[TA,] 
        # Estimate the parameters of the target domain profile using transfer learning
        b1 <- trans_fit(D.t = D.t, D.s = D)
      }
      cov1 <- cov_fit(X, b0)  # Calculate covariance matrix
      TA1 <- TA1 + length(TA) # Calculate cumulative transferred sample size
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
        TA2 <- TA1%/%i
        break
      }
    }
    
    # If the statistic has not exceeded the control limit until the end, set the run length to the total number of profiles in the target domain
    if(RL == 0){
      RL <- lt
      TA2 <- TA1%/%lt
    }
    
    # Return the run length as a value of parallel computing
    my_list1 <- list(RL, TA2)
    names(my_list1) <- c("RL", "TA2")
    return(my_list1)
  }
  stopCluster(cl)  # Terminate parallel computing
  
  t2 <- proc.time()  # Record end time
  t2-t1  # Total running time
  
  RL <- as.numeric(MEWMA[,1])
  TA2 <- as.numeric(MEWMA[,2])
  
  filtered_RL <- RL[RL > tau] # Keep ARL greater than shift time tau
  filtered_RL <- filtered_RL - tau # Average run length from change point to alarm
  SSARL1 <- mean(filtered_RL, na.rm = TRUE) # Steady-state SSARL1
  SSSD <- sd(filtered_RL, na.rm = TRUE) # Steady-state SSSD
  
  TA2_1 <- mean(TA2, na.rm = TRUE) # Average transferred sample size
  
  my_list <- list(SSARL1, SSSD, TA2_1)
  names(my_list) <- c("SSARL1", "SSSD", "Average Transferred Sample Size")
  return(my_list)
}
