#####################################################
#This R-code can be used to perform a power analysis#
#with the ML-MIMIC model                            #
#####################################################

#Load the required packages
require(psych)
require(MplusAutomation)
require(texreg)


########################
#Start simulation study#
########################

#set.seed(123)
B <- 3

  #Set parameters for power analysis.  

  #Sample size (adjust these in each simulation run to achieve sufficient power
  n2 <- 100
  n1 <- 20
  N  <- n2 * n1  
  
  #Number of indicators (this requires some adaptions when creating the items)
  p <- 5

  #Main within and between effects
  mainw <- 0.4
  mainb <- 0.5
  
  #Interaction effect (effect size: cross-level interaction effect)
  int <- 0.2	
  
  #Factor loadings within and between-level
  lambda <- c(1.0, 0.95, 0.8, 0.9, 0.85)

  #Between-level factor variance
  tau <- 0.25
  
  #Within-level factor variance
  Pi <- 0.75

  #Random slope variance within-level covariate
  slope.var <- 0.50
  
  #Covariance random slope and random intercept
  cov.sl.int <- 0.1		
  
  #Residual variance between-level
  theta_b <- rep(0.05, p)
  
  #Residual variance within-level
  theta_w <- c(0.60, 0.19, 0.36, 0.42, 0.51)

  
  ###################################
  #Set up for loop to simulate power#	  
  ################################### 
  #Set output
  #Admissable solutions
  asr_mlmimic <- 0
  
  #P-value, standard error and estimate interaction effects of ML-MIMIC model
  p.mlmimic       <- rep(NA, B)
  est.cr.mlmimic  <- rep(NA, B)
  se.cr.mlmimic   <- rep(NA, B)
  
  #Start for loop		
  i   <- 1           #Count number of admissable solutions
  tot <- 0           #Count total number of replications
  while (i <= B) {
    #Continue the loop until B admissable solutions
 
    tot <- tot + 1
    
    #Create dataframe
    #Create covariate variables
    
    #Group covariate level 1
    L1g <- rep(sort(rep(c(0, 1), n1 / 2)), n2)
    
    #Group coviariate level 2
    L2g <- sort(rep(c(0, 1), N / 2))
    
    #Create patient ID variable
    id  <- sort(rep(1:n2, n1))
    
    #Create data frame
    df  <- data.frame(cbind(id, L2g, L1g))
    
    #Create model matrix (without intercept) to simulate covariate effects
    X   <- model.matrix(~ -1 + L2g * L1g, data = df)
    
    #Set (residual) variance of the latent variables: dzeta's
    #dzeta at level 2: 2x2 matrix for random intercept (tau) and random slope var and covariance
    dzeta      <- matrix(c(tau, cov.sl.int, cov.sl.int, slope.var), 2, 2)
    dzeta_i    <- MASS::mvrnorm(n2, c(0, 0), dzeta)
    
    #dzeta at level 1 (variance of Pi)
    dzeta_ij   <- matrix(rnorm(N, 0, sqrt(Pi)))
    
    #Repeat the level 2 dzeta's for each level 2 subject, n1 times
    dzeta_i.r <- dzeta_i[rep(seq_len(nrow(dzeta_i)), each = n1), ]
    
    #Create latent variable at within-level
    eta_ij    <-
      matrix(X[, 2] * mainw + X[, 3] * int + X[, 2] * dzeta_i.r[, 2] + dzeta_ij)
    
    #Create latent variable at between-level
    eta_i     <- matrix(X[, 1] * mainb + dzeta_i.r[, 1])
    
    tmp <- cbind(eta_ij, eta_i)
    
    #Create between-level residual variance matrix (for each person)
    tscoresb <- MASS::mvrnorm(n2, rep(0, p), diag(theta_b))
    tscoresb <- tscoresb[rep(seq_len(nrow(tscoresb)), each = n1), ]
    
    #Create within-level residual matrix (for each observation)
    tscoresw <- MASS::mvrnorm(N,  rep(0, p), diag(theta_w))
    
    tmp.data <- data.frame(cbind(tmp, tscoresw, tscoresb))
    
    #Simulate the observed scores y1 ... yp
    y_tmp <- as.data.frame(matrix(NA, ncol=p, nrow=N))
    for (j in 1:p){		
    	y_tmp[,j] <- lambda[j] * tmp[, 1] + tscoresw[, j] + lambda[j] * tmp[, 2] + tscoresb[, j]
    }    
    #Create data
    data  <- cbind(y_tmp, L1g, L2g, id)
    cnames <- paste0("y", 1:p)
    colnames(data) <-
      c(cnames, "L1grp", "L2grp", "id") 

    mplusdata      <- data 
    
    ####################
    #Fit ML-MIMIC model#
    ####################
    
    #Fit ML-MIMIC model in Mplus using MplusAutomation and load parameters into R
    #Be sure that the Mplus code specifies the same number of observed indicators 	    
	
    prepareMplusData(mplusdata, "./Mplus/mlmimic.dat")
    runModels("./Mplus/mlmimic_randomslope.inp")
    mlmimic <- readModels("./Mplus/mlmimic_randomslope.out")
    parammlmimic <- mlmimic$parameters$unstandardized
        
    ###################################
    #Check ASR solutions ML-MMIC model#
    ###################################
    
    #Admissable solution check ML-MIMIC model
    if ((length(mlmimic$errors) == 0) &
        (length(mlmimic$warnings) == 0) &
        (sum(colnames(parammlmimic) == "est") == 1) &
        (sum(colnames(parammlmimic) == "se") == 1)) {
      if ((sum(parammlmimic[, "se"] < 0) == 0) &
          (sum(parammlmimic[parammlmimic$paramHeader == "Residual.Variances", "est"] <
               0) == 0) &
          (parammlmimic[parammlmimic$paramHeader == "New.Additional.Parameters", "est"] >=
           -1 &
           parammlmimic[parammlmimic$paramHeader == "New.Additional.Parameters", "est"] <=
           1))
      {
        c              <- 1
        
        #Count the number of admissable solutions within the B planned replications
        if (tot <= B) {
          asr_mlmimic    <- asr_mlmimic + 1
        }
      }
    } else{
      c              <- 0
      if (tot <= B) {
        asr_mlmimic <- asr_mlmimic
      }
    }
        
    #Store parameters if admissable solution is obtained
    if (c == 1) {
      p.mlmimic[i]         <-
        parammlmimic[parammlmimic$paramHeader == "S1.ON", "pval"]
      est.cr.mlmimic[i]    <-
        parammlmimic[parammlmimic$paramHeader == "S1.ON", "est"]
      se.cr.mlmimic[i]     <-
        parammlmimic[parammlmimic$paramHeader == "S1.ON", "se"]
            
      #Increment i
      i                  <- i + 1
    }
  }
 
  ################# 	
  #Calculate POWER#
  #################
  power <- sum(p.mlmimic < 0.05) / B
   