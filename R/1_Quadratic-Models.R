#' Takes germination trial data,
#' runs mixed effect logistic regression models for each species*location (Grid.ID),
#' then creates temeprature estimates and their standard errors for:
#' Tmax - Temperature where model predicts 5% overall germination from model (upper),
#' Tmin - Temperature where model predicts 5% overall germination from model (lower),
#' Topt - Temperature where model predicts maximum germination,
#' Topt.upp - Temperature where model predicts 95% of maximum germination from model (upper),
#' Topt.low - Temperature where model predicts 95% of maximum germination from model (lower),
#' Tbreadth - Tmax - Tmin
#'
#' Standard Errors for each of the above are calculated using the delta method
#'
#' For Seed Germination over Latitude Project
#'
#+ quadratic-models, eval = FALSE


getQuadModels <-  function() {
  #Load packages
  require(tidyr)
  require(ggplot2)
  require(dplyr)
  require(lme4)
  require(broom)
  require(purrr)
  
  # Steps for pre-model filtering:
  # - Open cleaned data file and make it a table.
  # - Data file must include at least:
  # "Grid.ID", "Test.Temp", "NumSown", "NumGerm", "Grid.Lat"
  # - Add in Proportion germinated
  # - Find maximum germination percentage in group
  # - Filter out any groups with no germination lower than
  #  0.25 of maximum germiantion % (too flat)
  # - Filter out any groups with no germination above 10% (too low)
  # - Filter out any groups with only one value above 0 germination (not enough data)
  
  MSBP <- read.csv("./Outputs/Germ_Filtered.csv") %>%
    tbl_df() %>%
    rename(Observation = X) %>% #add a term for observation for models below
    mutate(propgerm = NumGerm / NumSownOrig)
  
  #Find maximum germination percentage in group
  MSBP <- MSBP %>%
    group_by(Grid.ID)  %>%
    mutate(maxgerm = max(propgerm))
  
  #Filter out any groups with no germination lower than
  # 0.25 of maximum germiantion % (too flat)
  MSBP <- MSBP %>%
    group_by(Grid.ID)  %>%
    filter(any(propgerm < 0.25 * maxgerm))
  
  #Filter out any groups with no germination higher than 10% (too low)
  MSBP <- MSBP %>%
    group_by(Grid.ID)  %>%
    filter(any(propgerm < 0.1))
  
  #Filter out any groups with only one value above 0 germination (not enough data)
  MSBP <- MSBP %>%
    group_by(Grid.ID)  %>%
    add_tally(propgerm > 0) %>%
    filter(n > 1)
  
  # Estimating Tmin, Tmax and Topt (and Topt.upp and Topt.low) using Quadratic Model
  #
  # Steps:
  # - Group by Grid.ID and perform quadratic binomial models on each group
  # - This is done 'safely' to avoid errors (no model returned if error occurs)
  # - Filter out any concave up fits (can't predict cardinal values)
  # - Calulate Topt, Tmin, Tmax
  # - Calculate standard error using coefficents and variance covariance matrix
  
  #Method for creating a list of models, or errors if they not not converge
  
  ctrl <-
    glmerControl(optimizer = "bobyqa", 
                 optCtrl = list(maxfun = 100000)) #control parameters for glmer
  
  #continues if there is an error (returns the error if there is one) package:purrr
  safe_glm <- safely(glmer) 
  
  MSBP_Models <- MSBP %>%
    group_by(Grid.ID)  %>%
    do(germ.model = safe_glm(
      cbind(NumGerm, NumSownOrig - NumGerm) ~
        I(Test.Temp)  + I((Test.Temp) ^ 2) + (1 | Observation),
      family = binomial,
      control = ctrl,
      nAGQ = 0,
      data = .
    ))
  
  
  #Extract coefficients a, b, and c from models
  MSBP_MaxMinOpt <- MSBP_Models %>%
    rowwise() %>%
    mutate(
      a = ifelse(is.null(germ.model[[1]]), NA, fixef(germ.model[[1]])[3]),
      b = ifelse(is.null(germ.model[[1]]), NA, fixef(germ.model[[1]])[2]),
      c = ifelse(is.null(germ.model[[1]]), NA, fixef(germ.model[[1]])[1])
    )
  
  #Remove concave up fits (a<0) and models without fits (a is NA)
  MSBP_MaxMinOpt <- filter(MSBP_MaxMinOpt, a < 0)
  
  #Add in Topt, Tmax, Tmin estimates from coefficients
  MSBP_MaxMinOpt <- MSBP_MaxMinOpt %>%
    mutate(
      Topt = (-b / (2 * a)),
      Logit.Max = a * Topt * Topt + b * Topt + c,
      #Maximum germination from model at Topt (logit)
      Model.Max = (exp(Logit.Max) / (1 + exp(Logit.Max))),
      #convert to actual maximum germiantion
      Topt.upp = (-b - (sqrt((
        b * b - 4 * a * (c - log(0.95 * Model.Max / (1 - 0.95 * Model.Max)))
      )))) / (2 * a),
      Topt.low = (-b + (sqrt((
        b * b - 4 * a * (c - log(0.95 * Model.Max / (1 - 0.95 * Model.Max)))
      )))) / (2 * a),
      Tmax = (-b - (sqrt((
        b * b - 4 * a * (c - log(0.05 / (1 - 0.05)))
      )))) / (2 * a),
      Tmin = (-b + (sqrt((
        b * b - 4 * a * (c - log(0.05 / (1 - 0.05)))
      )))) / (2 * a),
      Tbreadth = Tmax - Tmin,
      Toptbreadth = Topt.upp - Topt.low
    )
  
  # Delta method  - standard error estimate
  se.matrix.fun <-
    function(a,
             b,
             c,
             Model.Max,
             prop,
             upp.low,
             germ.model,
             is.opt) {
      if (is.null(germ.model[[1]])) {NA}
      else {
        x <-
          ifelse(is.opt == T, 
                 log(prop * Model.Max / (1 - prop * Model.Max)), #If calulating Topt, use model maximum
                 log(prop / (1 - prop))) #If calculating Min/Max, use total proportion
        
        delta <- b * b - 4 * c * a + 4 * x * a
        
        #Original error calculations
    pa <- ifelse(upp.low == "upp",
      1 / (2 * a ^ 2) * (b + (delta ^ (1 / 2)) + 2 * a * (c - x) * (delta ^ (-1 / 2))), #Topt.upp/Tmax
      1 / (2 * a ^ 2) * (b - (delta ^ (1 / 2)) - 2 * a * (c - x) * (delta ^ (-1 / 2)))) #Topt.low/Tmin
    
    pb <- ifelse(upp.low == "upp", 
                 -1 / (2 * a) * (1 + b * (delta ^ (-1 / 2))), #Topt.upp/Tmax
                 -1 / (2 * a) * (1 - b * (delta ^ (-1 / 2))))#Topt.low/Tmin
  
    pc <- ifelse(upp.low == "upp", (delta ^ (-1 / 2)), #Topt.upp/Tmax
                 -(delta ^ (-1 / 2)))  #Topt.low/Tmin

    vector.v <- as.vector(c(pc, pb, pa))
    vector.h <- t(vector.v) #transposed vector
    vcov.matrix <- as.matrix(vcov(germ.model[[1]]))
    variance <- vector.h %*% vcov.matrix %*% vector.v
    variance <- ifelse(variance < 0, NA, sqrt(as.numeric(variance)))
      }
    }
  
  se.opt.fun <- function(a, b, germ.model) {
    if (is.null(germ.model[[1]])) {
      NA
    }
    else {
      pa <- b / (2 * a * a)
      pb <- -1 / (2 * a)
      pc <- 0
      
      vector.v <- as.vector(c(pc, pb, pa))
      vector.h <- t(vector.v) #transposed vector
      vcov.matrix <- as.matrix(vcov(germ.model[[1]]))
      vcov.matrix[1, ] <- 0
      vcov.matrix[, 1] <- 0
      variance <- vector.h %*% vcov.matrix %*% vector.v
      variance <- ifelse(variance < 0, NA, sqrt(as.numeric(variance)))
    }
  }
  
  se.breadth.fun <- function(a, b, c, Model.Max, prop, germ.model, is.opt) {
    if (is.null(germ.model[[1]])) {
      NA
    }
    else {
      x <-
        ifelse(is.opt == T, 
               log(prop * Model.Max / (1 - prop * Model.Max)), #If calulating Topt breadth, use model maximum
               log(prop / (1 - prop))) #If calculating breadth, use total proportion

      delta <- b * b - 4 * c * a + 4 * x * a
      
      ## Corrected formula
      pa <- (-1 / (a * a)) * (sqrt(delta) + (2 * a * (c - x) / (sqrt(delta))))
      pb <- b / (a * sqrt(delta))
      pc <- -2 / sqrt(delta)
      
      vector.v <- as.vector(c(pc, pb, pa))
      vector.h <- t(vector.v) #transposed vector
      vcov.matrix <- as.matrix(vcov(germ.model[[1]]))
      variance <- vector.h %*% vcov.matrix %*% vector.v
      variance <- ifelse(variance < 0, NA, sqrt(as.numeric(variance)))
    }
  }
  
  MSBP_MaxMinOpt <- MSBP_MaxMinOpt %>%
    rowwise() %>%
    mutate(
      SEmax = se.matrix.fun(a, b, c, Model.Max, 0.05, "upp", germ.model, F),
      SEmin = se.matrix.fun(a, b, c, Model.Max, 0.05, "low", germ.model, F),
      SEopt.upp = se.matrix.fun(a, b, c, Model.Max, 0.95, "upp", germ.model, T),
      SEopt.low = se.matrix.fun(a, b, c, Model.Max, 0.95, "low", germ.model, T),
      SEopt = se.opt.fun(a, b, germ.model),
      SEbreadth = se.breadth.fun(a, b, c, Model.Max, 0.05, germ.model, F),
      SEoptbreadth = se.breadth.fun(a, b, c, Model.Max, 0.95, germ.model, T)
    )
  
  #Remove models
  MSBP_MaxMinOpt <- dplyr::select(MSBP_MaxMinOpt,-germ.model)
  
  MSBP_MaxMinOpt <- left_join(MSBP, MSBP_MaxMinOpt, "Grid.ID")
  
  #If there is no germination proportion above Topt that is less than 25% maxgerm, 
  #ignore Tmax, Topt.upp, Topt.low and Topt
  #If there is no germination proportion below Topt that is less than 25% maxgerm,
  #ignore Tmin, Topt.upp, Topt.low and Topt
  MSBP_MaxMinOpt <- MSBP_MaxMinOpt %>%
    group_by(Grid.ID)  %>%
    mutate(
      Tmax = ifelse(any((Test.Temp > Topt & propgerm < 0.25 * maxgerm)), Tmax, NA),
      Tmin = ifelse(any((Test.Temp < Topt & propgerm < 0.25 * maxgerm)), Tmin, NA),
      Topt = ifelse((any((Test.Temp > Topt & propgerm < 0.25 * maxgerm)) &
        any((Test.Temp < Topt & propgerm < 0.25 * maxgerm))), Topt, NA),
      Topt.upp = ifelse((any((Test.Temp > Topt & propgerm < 0.25 * maxgerm)) &
        any((Test.Temp < Topt & propgerm < 0.25 * maxgerm))), Topt.upp, NA),
      Topt.low = ifelse((any((Test.Temp > Topt &  propgerm < 0.25 * maxgerm)) &
        any((Test.Temp < Topt & propgerm < 0.25 * maxgerm))), Topt.low, NA),
      Tbreadth = ifelse((any((Test.Temp > Topt & propgerm < 0.25 * maxgerm)) &
        any((Test.Temp < Topt & propgerm < 0.25 * maxgerm))), Tbreadth, NA),
      Toptbreadth = ifelse((any((Test.Temp > Topt & propgerm < 0.25 * maxgerm)) &
        any((Test.Temp < Topt & propgerm < 0.25 * maxgerm))), Tbreadth, NA)
    )
  
  
  write.csv(MSBP_MaxMinOpt, file = "./Outputs/MSBP_MaxMinOpt.csv")
}

getQuadModels()
