############################################################
#This R-code is used to create the tables in the manuscript#
############################################################

require(xtable)

#Read tables
outputmlmimic <- read.table("./Output data/outputmlmimic.txt", sep = "\t", dec = ".")
outputml      <- read.table("./Output data/outputml.txt", sep = "\t", dec = ".")
outputaov     <- read.table("./Output data/outputaov.txt", sep = "\t", dec = ".")

outputmlmimic <- outputmlmimic[order(outputmlmimic$gamma,
                                     outputmlmimic$tau,
                                     outputmlmimic$n2,
                                     outputmlmimic$n1), ]
outputml      <- outputml[order(outputml$gamma, 
                                outputml$tau, 
                                outputml$n2, 
                                outputml$n1), ]
outputaov     <- outputaov[order(outputaov$gamma, 
                                 outputaov$tau, 
                                 outputaov$n2, 
                                 outputaov$n1), ]

############################################################
#Function to make the latex table of the simulation results#
# This function produces latex code that need to be copied #
# in the manuscript.                                       #
# Commands:                                                #
# - inter = whether the interaction effect must be included#
# - select.vars = which variables to put in the columns    #
# - theta   = which variance within settings?              #
############################################################

createtable <- function(inter, select.vars) {
  merge.vars <- select.vars[-length(select.vars)]
  comb       <- 1:3
  tablecomb <- list(0)
  
  for (i in comb) {
    if (inter == FALSE) {
      table1.1 <- subset(outputmlmimic,
                         (gamma == 0.0 & loadings == i),
                         select = c(select.vars))
      
      table1.2 <- subset(outputml,
                         (gamma == 0.0 & loadings == i),
                         select = c(select.vars))
      
      table1.3 <- subset(outputaov,
                         (gamma == 0.0 & loadings == i),
                         select = c(select.vars))
    } else if (inter == TRUE) {
      table1.1 <- subset(outputmlmimic,
                         (gamma != 0.0 & loadings == i),
                         select = c(select.vars))
      
      table1.2 <- subset(outputml,
                         (gamma != 0.0 & loadings == i),
                         select = c(select.vars))
      
      table1.3 <- subset(outputaov,
                         (gamma != 0.0 & loadings == i),
                         select = c(select.vars))
    }
    tablecomb[[i]] <-
      merge(
        merge(
          table1.1,
          table1.2,
          by = merge.vars,
          all = TRUE,
          sort = F,
          suffixes = c("mimic", "ml")
        ),
        table1.3,
        by = merge.vars,
        sort = F,
        all = TRUE
      )
  }
  tablecomb2 <-
    merge(
      merge(
        tablecomb[[1]],
        tablecomb[[2]],
        by = merge.vars,
        sort = F,
        all = T,
        suffixes = c("1-1", "2-1")
      ),
      tablecomb[[3]],
      by = merge.vars,
      sort = F,
      all = T
    )
  
  tablecomb2 <- tablecomb2[order(tablecomb2$varw), ]
  
  tablecomb2$empty  <- NA
  tablecomb2$N      <- paste0(tablecomb2$n2, "/", tablecomb2$n1)
  tablecomb2$tau    <- formatC(tablecomb2$tau, digits = 2, format = "f")
  
  if ("gamma" %in% select.vars) {
    tablecomb2$gamma <- formatC(tablecomb2$gamma, digits = 2, format = "f")
    tablecomb2$tau[duplicated(tablecomb2[, c("gamma", "varw", "tau")])] <-
      NA
    tablecomb2$gamma[duplicated(tablecomb2[, c("gamma", "varw")])] <-
      NA
    tablecomb2$varw[duplicated(tablecomb2$varw)] <- NA
    tablecomb2$varw   <- ifelse(tablecomb2$varw == 1,
                                "EQ",
                                ifelse(tablecomb2$varw == 2, "NEQ", tablecomb2$varw))
    
    tablefinal <- tablecomb2[, c(
      1,
      2,
      3,
      ncol(tablecomb2),
      6,
      7,
      8,
      ncol(tablecomb2) - 1,
      9,
      10,
      11,
      ncol(tablecomb2) - 1,
      12,
      13,
      14
    )]
  } else{
    tablecomb2$tau[duplicated(tablecomb2[, c("varw", "tau")])] <- NA
    tablecomb2$varw[duplicated(tablecomb2$varw)] <- NA
    tablecomb2$varw   <- ifelse(tablecomb2$varw == 1,
                                "EQ",
                                ifelse(tablecomb2$varw == 2, "NEQ", tablecomb2$varw))
    tablefinal <- tablecomb2[, c(
      1,
      2,
      ncol(tablecomb2),
      5,
      6,
      7,
      ncol(tablecomb2) - 1,
      8,
      9,
      10,
      ncol(tablecomb2) - 1,
      11,
      12,
      13
    )]
  }
  print(
    xtable(tablefinal, digits = 3),
    include.rownames = F,
    include.colnames = F,
    only.content = T
  )
}

#Table 2
asr       <-
  createtable(inter = T,
              select.vars = c("varw", "gamma", "tau", "n2", "n1", "ASR"))

#Table 3
typ1error <-
  createtable(inter = F,
              select.vars = c("varw", "tau", "n2", "n1", "Pwr"))

#Table 4
power     <-
  createtable(inter = T,
              select.vars = c("varw", "gamma", "tau", "n2", "n1", "Pwr"))

#Table 5
se        <-
  createtable(inter = T,
              select.vars = c("varw", "gamma", "tau", "n2", "n1", "secr"))

#Table 6
rb        <-
  createtable(inter = T,
              select.vars = c("varw", "gamma", "tau", "n2", "n1", "RBcr"))

#Table 7
rseb      <-
  createtable(
    inter = T,
    select.vars = c("varw", "gamma", "tau", "n2", "n1", "RSEBcr")
  )

#Table 8
mse       <-
  createtable(inter = T,
              select.vars = c("varw", "gamma", "tau", "n2", "n1", "MSEcr")
              