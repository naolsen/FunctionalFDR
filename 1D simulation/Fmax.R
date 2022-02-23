###################################################################################################################################################
# Fmax test
# modified to perform F test with permutation of responses only!
###################################################################################################################################################


Fmax <- function(formula,
                  B = 1000,
                  method = 'residuals',
                  dx=NULL,
                  recycle=TRUE){
  
  extract_residuals <- function(regr){
    return(regr$residuals)
  }
  extract_fitted <- function(regr){
    return(regr$fitted)
  }
  
  env <- environment(formula)
  variables = all.vars(formula)
  y.name = variables[1]
  covariates.names <- colnames(attr(terms(formula),"factors"))
  #data.all = model.frame(formula)
  cl <- match.call()
  data <- get(y.name,envir = env)
  if(is.fd(data)){ # data is a functional data object
    rangeval <- data$basis$rangeval
    if(is.null(dx)){
      dx <- (rangeval[2]-rangeval[1])*0.01
    }
    abscissa <- seq(rangeval[1],rangeval[2],by=dx)
    coeff <- t(eval.fd(fdobj=data,evalarg=abscissa))
  }else if(is.matrix(data)){
    coeff <- data
  }else{
    stop("First argument of the formula must be either a functional data object or a matrix.")
  }
  
  dummynames.all <- colnames(attr(terms(formula),"factors"))
  formula.const <- deparse(formula[[3]],width.cutoff = 500L) #extracting the part after ~ on formula. this will not work if the formula is longer than 500 char
  
  formula.discrete <- as.formula(paste('coeff ~',formula.const),env=environment())
  design_matrix = model.matrix(formula.discrete)
  mf = model.frame(formula.discrete)
  
  nvar <- dim(design_matrix)[2] - 1
  var_names <- colnames(design_matrix)
  p <- dim(coeff)[2]
  n <- dim(coeff)[1]
  # Univariate permutations
  regr0 <- lm.fit(design_matrix, coeff)
  # Test statistics
  Sigma <- chol2inv(regr0$qr$qr)
  resvar <- colSums(regr0$residuals ^ 2) / regr0$df.residual
  se <- sqrt(matrix(diag(Sigma), nrow = nvar + 1, ncol = p, byrow = FALSE) 
             * matrix(resvar, nrow = nvar + 1, ncol = p, byrow = TRUE))
  T0_part <- abs(regr0$coeff / se)^2
  if (nvar > 0) {
    T0_glob <- colSums((regr0$fitted - matrix(colMeans(regr0$fitted),
                                              nrow = n, ncol = p, 
                                              byrow = TRUE)) ^ 2) / (nvar * resvar)
  } else {
    method <- 'responses' 
    T0_glob <- numeric(p)
    T0_part <- t(as.matrix(T0_part))
  }

  # CMC algorithm
  T_glob <- matrix(ncol = p,nrow = B)
  T_part <- array(dim = c(B, nvar + 1, p))
  for (perm in 1:B) {

    permutazioni <- sample(n)
    coeff_perm <- coeff[permutazioni, ]
    regr_perm <- lm.fit(design_matrix, coeff_perm)
    Sigma <- chol2inv(regr_perm$qr$qr)
    resvar <- colSums(regr_perm$residuals ^ 2) / regr_perm$df.residual
    
    T_glob[perm, ] <- colSums((regr_perm$fitted - matrix(colMeans(regr_perm$fitted), 
                                                           nrow = n, ncol = p,
                                                           byrow=TRUE)) ^ 2)/ (nvar * resvar)
  }
  pval_glob <- numeric(p)
  pval_glob_adj <- numeric(p)
  maxF_glob = apply(T_glob,1,max)
  
  
  for (i in 1:p) {
    pval_glob[i] <- sum(T_glob[, i] >= T0_glob[i]) / B
    pval_glob_adj[i] = sum(maxF_glob >= T0_glob[i]) / B
  }
  
  
  
  result <- list(call=cl,
                 unadjusted_pval_F=pval_glob,
                 adjusted_pval_F=pval_glob_adj)
  
  return(result)
}
