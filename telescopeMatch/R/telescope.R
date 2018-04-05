
############################
### Baseline telescope matching routine - wraps the Matching package
###########################
telescopeMatch <- function(outcome, outcome.formula, treatment, mediator, pre.treatment, post.treatment, data, K.fs=3, K.ss = 3, BiasCorrect=T, nBoot=100){

  ### Sanity check - no pre-treatment covariates in post-treatment and vice-versa.
  if (any(pre.treatment %in% post.treatment)|any(post.treatment %in% pre.treatment)){
    stop("Warning: pre.treatment and post.treatment covariates overlap")
  }

  all.covariates <- c(pre.treatment, post.treatment)

  ### W/in Treated group, match M = 0 to M = 1
  reg.treated = lm(outcome.formula, data=data[data[[treatment]]==1&data[[mediator]]==0,])
  reg.control = lm(outcome.formula, data=data[data[[treatment]]==0&data[[mediator]]==0,])

  ### Predicted Y
  data$pred.Y.M0 <- NA
  data[data[[treatment]]==1,]$pred.Y.M0 <- predict(reg.treated, newdata=data[data[[treatment]]==1,][,all.covariates])
  data[data[[treatment]]==0,]$pred.Y.M0 <- predict(reg.control, newdata=data[data[[treatment]]==0,][,all.covariates])

  fs.1 = Match(Y = data[[outcome]][data[[treatment]] == 1], Tr = data[[mediator]][data[[treatment]] == 1],
               X = data[data[[treatment]] == 1,][,all.covariates], M=K.fs, BiasAdjust = BiasCorrect)

  ### W/in Control group, Match M = 0 to M = 1
  fs.0 = Match(Y = data[[outcome]][data[[treatment]] == 0], Tr = data[[mediator]][data[[treatment]] == 0],
               X = data[data[[treatment]] == 0,][,all.covariates], M=K.fs, BiasAdjust = BiasCorrect)

  fs.1.save <- fs.1
  fs.0.save <- fs.0

  ### Impute the outcome using each match set using the bias-correction
  data$YM0 = NA
  ### for units with observed mediator 0, this is just the observed outcome
  data$YM0[data[[treatment]] == 1&data[[mediator]] == 0] <- data[[outcome]][data[[treatment]] == 1&data[[mediator]] == 0]
  data$YM0[data[[treatment]] == 0&data[[mediator]] == 0] <- data[[outcome]][data[[treatment]] == 0&data[[mediator]] == 0]

  #### For each treated unit w/ M=1, match K.fs units with M=0 to M=1 - subtract off the regression bias correction
  impute.treatment <- tapply(fs.1$index.control, fs.1$index.treated, function(x) mean(data[data[[treatment]] == 1,][[outcome]][x] - data[data[[treatment]] == 1,]$pred.Y.M0[x]))
  data$YM0[data[[treatment]] == 1&data[[mediator]] == 1] <- data$pred.Y.M0[data[[treatment]] == 1&data[[mediator]] == 1] + impute.treatment

  ### Do same for control
  impute.control <- tapply(fs.0$index.control, fs.0$index.treated, function(x) mean(data[data[[treatment]] == 0,][[outcome]][x] - data[data[[treatment]] == 0,]$pred.Y.M0[x]))
  data$YM0[data[[treatment]] == 0&data[[mediator]] == 1] <- data$pred.Y.M0[data[[treatment]] == 0&data[[mediator]] == 1] + impute.control

  ### Second-stage, match on de-mediated outcome, only pre-treatment covariates


  ss <- Match(Y = data$YM0, Tr = data[[treatment]], X= data[,pre.treatment], M = K.ss, BiasAdjust = BiasCorrect)
  estimate <- ss$est

  ss.save <- ss
  
  if(nBoot > 1){
    boot.ests <- rep(NA, nBoot)
    for (boot in 1:nBoot){

      boot.data <- data[sample(1:nrow(data), nrow(data), replace=T),]
      boot.ests[boot] <- tryCatch({

      ### W/in Treated group, match M = 0 to M = 1
      boot.reg.treated = lm(outcome.formula, data=boot.data[boot.data[[treatment]]==1&boot.data[[mediator]]==0,])
      boot.reg.control = lm(outcome.formula, data=boot.data[boot.data[[treatment]]==0&boot.data[[mediator]]==0,])

      ### Predicted Y
      boot.data$pred.Y.M0 <- NA
      boot.data[boot.data[[treatment]]==1,]$pred.Y.M0 <- predict(boot.reg.treated, newdata=boot.data[boot.data[[treatment]]==1,][,all.covariates])
      boot.data[boot.data[[treatment]]==0,]$pred.Y.M0 <- predict(boot.reg.control, newdata=boot.data[boot.data[[treatment]]==0,][,all.covariates])

      fs.1 = Match(Y = boot.data[[outcome]][boot.data[[treatment]] == 1], Tr = boot.data[[mediator]][boot.data[[treatment]] == 1],
                   X = boot.data[boot.data[[treatment]] == 1,][,all.covariates], M=K.fs, BiasAdjust = BiasCorrect)

      ### W/in Control group, Match M = 0 to M = 1
      fs.0 = Match(Y = boot.data[[outcome]][boot.data[[treatment]] == 0], Tr = boot.data[[mediator]][boot.data[[treatment]] == 0],
                   X = boot.data[boot.data[[treatment]] == 0,][,all.covariates], M=K.fs, BiasAdjust = BiasCorrect)


      ### Impute the outcome using each match set using the bias-correction
      boot.data$YM0 = NA
      ### for units with observed mediator 0, this is just the observed outcome
      boot.data$YM0[boot.data[[treatment]] == 1&boot.data[[mediator]] == 0] <- boot.data[[outcome]][boot.data[[treatment]] == 1&boot.data[[mediator]] == 0]
      boot.data$YM0[boot.data[[treatment]] == 0&boot.data[[mediator]] == 0] <- boot.data[[outcome]][boot.data[[treatment]] == 0&boot.data[[mediator]] == 0]

      #### For each treated unit w/ M=1, match K.fs units with M=0 to M=1 - subtract off the regression bias correction
      impute.treatment <- tapply(fs.1$index.control, fs.1$index.treated, function(x) mean(boot.data[boot.data[[treatment]] == 1,][[outcome]][x] - boot.data[boot.data[[treatment]] == 1,]$pred.Y.M0[x]))
      boot.data$YM0[boot.data[[treatment]] == 1&boot.data[[mediator]] == 1] <- boot.data$pred.Y.M0[boot.data[[treatment]] == 1&boot.data[[mediator]] == 1] + impute.treatment

      ### Do same for control
      impute.control <- tapply(fs.0$index.control, fs.0$index.treated, function(x) mean(boot.data[boot.data[[treatment]] == 0,][[outcome]][x] - boot.data[boot.data[[treatment]] == 0,]$pred.Y.M0[x]))
      boot.data$YM0[boot.data[[treatment]] == 0&boot.data[[mediator]] == 1] <- boot.data$pred.Y.M0[boot.data[[treatment]] == 0&boot.data[[mediator]] == 1] + impute.control


      ss.boot <- Match(Y = boot.data$YM0, Tr = boot.data[[treatment]], X= boot.data[,pre.treatment], M = K.ss, BiasAdjust = BiasCorrect)
      ss.boot$est
      },
      error=function(cond) {
        message(paste("Warning: Bootstrap iteration", boot, "encountered an error"))
        message("Here's the original error message:")
        message(cond)
        # Choose a return value in case of error
        return(NA)
      },
      warning=function(cond) {
        message(paste("Warning: Bootstrap iteration", boot, "encountered a warning"))
        message("Here's the original error message:")
        message(cond)
        # Choose a return value in case of warning
        return(NA)
      })
    }

    ### Calculate bootstrap SE
    standardError <- sd(boot.ests, na.rm=T)

  }else{
    boot.ests <- NA
    standardError <- NA
  }

  ### Return Match object from the second stage
  return(results = list(est = estimate, se = standardError, boot = boot.ests, stage2match = ss.save, stage1match0 = fs.1.save, stage1match1 = fs.0.save))

}
