##################################################
#This R-code was used to run the simulation study#
##################################################

#Load the required packages
require(psych)
require(MplusAutomation)
require(texreg)
require(lme4)
require(lmerTest)

#########################
#Simulation study set-up#
#########################

#Specify simulation conditions
#Factor loadings
loadings  <- c(1, 2, 3)

#Sample size
n2        <- c(50, 100, 200)
n1        <- c(10, 20)

#Between-factor variance
tau       <- c(0.10, 0.25)

#Within-level residual variances
varw      <- c(1, 2)

#Cross-level interaction effect
gamma     <- c(0.0, 0.2, 0.4)

#Combine all conditions
cond  <- expand.grid(gamma, tau, n2, n1, loadings, varw)
cond  <- cond[order(rev(cond$Var1), cond$Var2, cond$Var3, cond$Var4),]

colnames(cond) <- c("gamma", "tau", "n2", "n1", "loadings", "varw")

#Number of replications
B <- 1000

#Number of indicators
p <- 5

#Create data to store output for each model separately
outputmlmimic <-
  cbind(cond, as.data.frame(matrix(
    NA, nrow = nrow(cond), ncol = 8
  )))
outputml      <-
  cbind(cond, as.data.frame(matrix(
    NA, nrow = nrow(cond), ncol = 8
  )))
outputaov     <-
  cbind(cond, as.data.frame(matrix(
    NA, nrow = nrow(cond), ncol = 8
  )))

#Colnames of output files:
#ASR = Admissable solution rate
#Pwr = Power
#cr = cross-level interaction effect
#secr = Standard error cross-level interaction effect
#biascr = Bias cross-level interaction effect
#RBcr = Relative bias cross-level interaction effect
#RSEBcr = Relative standard error bias cross-level interaction effect
#MSEcr = Mean squared error cross-level interaction effect

colnames(outputmlmimic)[7:14] <-
  colnames(outputml)[7:14] <- colnames(outputaov)[7:14] <-
  c("ASR", "Pwr", "cr", "secr", "biascr", "RBcr", "RSEBcr",
    "MSEcr")

########################
#Start simulation study#
########################

set.seed(123)

#Loop through all k combinations
for (k in 1:nrow(cond)) {
  #Set parameters for each run
  
  #Main within and between effects
  mainw <- 0.4
  mainb <- 0.5
  
  #Interaction effect
  if (cond[k, "gamma"] == 0) {
    int          <- 0.0
  } else if (cond[k, "gamma"] == 0.2) {
    int          <- 0.2
  } else if (cond[k, "gamma"] == 0.4) {
    int          <- 0.4
  }
  
  #Factor loadings
  if (cond[k, "loadings"] == 1) {
    lambda <- rep(1.0, p)
  } else if (cond[k, "loadings"] == 2) {
    lambda <- c(1.0, 0.95, 0.8, 0.9, 0.85)
  } else if (cond[k, "loadings"] == 3) {
    lambda <- c(1.0, 0.9, 0.6, 0.8, 0.7)
  }
  
  #Sample size
  n2 <- cond[k, "n2"]
  n1 <- cond[k, "n1"]
  N  <- n2 * n1
  
  #Between-level factor variance
  tau <- cond[k, "tau"]
  
  #Within-level factor variance
  if (tau == 0.10) {
    Pi <- 0.90
  } else if (tau == 0.25) {
    Pi <- 0.75
  }
  
  #Residual variance between-level
  theta_b <- rep(0.05, p)
  
  #Residual variance within-level
  if (cond[k, "varw"] == 1) {
    theta_w <- rep(0.36, p)
  } else if (cond[k, "varw"] == 2) {
    theta_w <- c(0.60, 0.19, 0.36, 0.42, 0.51)
  }
  
  #Data Generation
  #Set output
  #Admissable solutions
  asr_mlmimic <- asr_ml <- asr_aov <- 0
  
  #P-value, standard error and estimate interaction effects of each model
  p.mlmimic       <- rep(NA, B)
  est.cr.mlmimic  <- rep(NA, B)
  se.cr.mlmimic   <- rep(NA, B)
  
  p.ml            <- rep(NA, B)
  est.cr.ml       <- rep(NA, B)
  se.cr.ml        <- rep(NA, B)
  
  p.aov           <- rep(NA, B)
  est.cr.aov      <- rep(NA, B)
  se.cr.aov       <- rep(NA, B)
  
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
    #dzeta at level 2: 2x2 matrix for random intercept (tau) and random slope (0.5) and zero covariance
    dzeta      <- matrix(c(tau, 0.1, 0.1, 0.5), 2, 2)
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
    
    #Simulate the observed scores y1, y2, y3, y4, y5
    y1 <- lambda[1] * tmp[, 1] + tscoresw[, 1] + lambda[1] * tmp[, 2] + tscoresb[, 1]
    y2 <- lambda[2] * tmp[, 1] + tscoresw[, 2] + lambda[2] * tmp[, 2] + tscoresb[, 2]
    y3 <- lambda[3] * tmp[, 1] + tscoresw[, 3] + lambda[3] * tmp[, 2] + tscoresb[, 3]
    y4 <- lambda[4] * tmp[, 1] + tscoresw[, 4] + lambda[4] * tmp[, 2] + tscoresb[, 4]
    y5 <- lambda[5] * tmp[, 1] + tscoresw[, 5] + lambda[5] * tmp[, 2] + tscoresb[, 5]
    
    #Create data
    y_tmp <- as.data.frame(cbind(y1, y2, y3, y4, y5))
    data  <- cbind(y_tmp, L1g, L2g, id)
    colnames(data) <-
      c("y1", "y2", "y3", "y4", "y5", "L1grp", "L2grp", "id")
    mplusdata      <- data
    
    #Sum score for regular multilevel model
    data$ysum      <- rowSums(data[, 1:p])
    
    #Aggregated data for ANOVA
    agg.y <-
      aggregate(cbind(ysum, L2grp) ~ id + L1grp,
                FUN = mean,
                data = data)
    agg.y <- agg.y[order(agg.y$id, agg.y$L1grp), ]
    
    ####################
    #Fit ML-MIMIC model#
    ####################
    
    #Fit ML-MIMIC model in Mplus using MplusAutomation and load parameters into R
    prepareMplusData(mplusdata, "./Mplus/mlmimic.dat")
    runModels("./Mplus/mlmimic_randomslope.inp")
    mlmimic <- readModels("./Mplus/mlmimic_randomslope.out")
    parammlmimic <- mlmimic$parameters$unstandardized
    
    ######################
    #Fit multilevel model#
    ######################
    
    #Fit a multilevel model on sum score (ysum)
    fit1 <-
      suppressMessages(lmer(
        ysum ~ L1grp * L2grp  + (1 + L1grp | id),
        REML = FALSE,
        data = data
      ))
    
    #Store the covariance matrix of the random effects and extract correlation random slope random intercept
    var.mat.ml  <- as.data.frame(VarCorr(fit1))[, 4][-3]
    cor <- as.data.frame(VarCorr(fit1))[3, 5]
    
    #################
    #Fit ANOVA model#
    #################
    
    #Fit ANOVA on aggregated dataset (agg.y)--> Analogous BWS ANOVA when REML = TRUE.
    fit2  <-
      suppressMessages(lmer(
        ysum ~ L1grp * L2grp + (1 | id),
        REML = TRUE,
        data = agg.y
      ))
    
    ####################################
    #Check ASR solutions for each model#
    ####################################
    
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
    
    #Admissable solutions check multilevel model
    if ((any(
      grepl("failed to converge", fit1@optinfo$conv$lme4$messages)
    ) == FALSE) &
    (isSingular(fit1) == FALSE) &
    (sum(summary(fit1)$coef[, 2] < 0) == 0) &
    (sum(var.mat.ml < 0) == 0) &
    (!is.na(cor))) {
      if (cor <= 1 & cor >= -1)
      {
        c              <- c + 1
        if (tot <= B) {
          asr_ml         <- asr_ml + 1
        }
      }
    } else{
      c              <- 0
      if (tot <= B) {
        asr_ml         <- asr_ml
      }
    }
    
    #Admissable solutions check ANOVA
    if ((any(
      grepl("failed to converge", fit2@optinfo$conv$lme4$messages)
    ) == FALSE) &
    (isSingular(fit2) == FALSE) &
    (sum(summary(fit2)$coef[, 2] < 0) == 0) &
    (sum(as.data.frame(VarCorr(fit2))[, 4] < 0) == 0)) {
      c              <- c + 1
      if (tot <= B) {
        asr_aov        <- asr_aov + 1
      }
    } else{
      c               <- 0
      if (tot <= B) {
        asr_aov         <- asr_aov
      }
    }
    
    #Store parameters if 3 admissable solutions are obtained (each for one model): p-value, estimate and standard error
    if (c == 3) {
      p.mlmimic[i]         <-
        parammlmimic[parammlmimic$paramHeader == "S1.ON", "pval"]
      est.cr.mlmimic[i]    <-
        parammlmimic[parammlmimic$paramHeader == "S1.ON", "est"]
      se.cr.mlmimic[i]     <-
        parammlmimic[parammlmimic$paramHeader == "S1.ON", "se"]
      
      p.ml[i]              <- summary(fit1)$coef[4, 5]
      est.cr.ml[i]         <- summary(fit1)$coef[4, 1] / p
      se.cr.ml[i]          <- summary(fit1)$coef[4, 2] / p
      
      p.aov[i]             <- summary(fit2)$coef[4, 5]
      est.cr.aov[i]        <- summary(fit2)$coef[4, 1] / p
      se.cr.aov[i]         <- summary(fit2)$coef[4, 2] / p
      
      #Increment i
      i                  <- i + 1
    }
  }
  #Calculate output per model
  #ASR
  outputmlmimic[k, "ASR"] <- asr_mlmimic / B
  outputml[k, "ASR"]      <- asr_ml / B
  outputaov[k, "ASR"]     <- asr_aov / B
  
  #Power
  outputmlmimic[k, "Pwr"] <- sum(p.mlmimic < 0.05) / B
  outputml[k, "Pwr"]      <- sum(p.ml < 0.05) / B
  outputaov[k, "Pwr"]     <- sum(p.aov < 0.05) / B
  
  #Point estimates
  outputmlmimic[k, "cr"] <- mean(est.cr.mlmimic)
  outputml[k, "cr"]      <- mean(est.cr.ml)
  outputaov[k, "cr"]     <- mean(est.cr.aov)
  
  #SE's
  outputmlmimic[k, "secr"] <- mean(se.cr.mlmimic)
  outputml[k, "secr"]      <- mean(se.cr.ml)
  outputaov[k, "secr"]     <- mean(se.cr.aov)
  
  #Bias
  outputmlmimic[k, "biascr"] <- mean(est.cr.mlmimic) - int
  outputml[k, "biascr"]      <- mean(est.cr.ml) - int
  outputaov[k, "biascr"]     <- mean(est.cr.aov) - int
  
  if (int != 0) {
    #Relative Bias
    outputmlmimic[k, "RBcr"] <- (mean(est.cr.mlmimic) - int) / int
    outputml[k, "RBcr"]      <- (mean(est.cr.ml) - int) / int
    outputaov[k, "RBcr"]     <- (mean(est.cr.aov) - int) / int
  }
  
  #Relative standard error Bias
  outputmlmimic[k, "RSEBcr"] <-
    (mean(se.cr.mlmimic) - sd(est.cr.mlmimic)) / sd(est.cr.mlmimic)
  outputml[k, "RSEBcr"]      <-
    (mean(se.cr.ml) - sd(est.cr.ml)) / sd(est.cr.ml)
  outputaov[k, "RSEBcr"]     <-
    (mean(se.cr.aov) - sd(est.cr.aov)) / sd(est.cr.aov)
  
  #MSE
  outputmlmimic[k, "MSEcr"] <-
    (mean(est.cr.mlmimic) - int) ^ 2 + sd(est.cr.mlmimic) ^ 2
  outputml[k, "MSEcr"]     <-
    (mean(est.cr.ml) - int) ^ 2 + sd(est.cr.ml) ^ 2
  outputaov[k, "MSEcr"]    <-
    (mean(est.cr.aov) - int) ^ 2 + sd(est.cr.aov) ^ 2
  
  write.table(outputmlmimic,
              "./output/outputmlmimic.txt",
              sep = "\t",
              dec = ".")
  write.table(outputml,
              "./output/outputml.txt",
              sep = "\t",
              dec = ".")
  write.table(outputaov,
              "./output/outputaov.txt",
              sep = "\t",
              dec = ".")
}
