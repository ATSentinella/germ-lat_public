#' Calculating approximate R2 from metafor objects
#+ metafor-R2, eval = FALSE


#We are using marginal  R2 here (fix/fix+random)
R2 <- function(model){
  
  # fixed effect variance
  fix <- var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))
  
  # marginal
  R2m <- fix / (fix + sum(model$sigma2))
  
  # conditional
  R2c <- (fix + sum(model$sigma2) - model$sigma2[length(model$sigma2)]) / 
    (fix + sum(model$sigma2))
  
  R2s <- c(R2_marginal = R2m, R2_conditional = R2c)
  
  return(R2s)
}
