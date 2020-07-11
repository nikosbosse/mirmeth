generatetoydata <- function(groups = 1, ip = T, n = 6, N = 50, size = 2,
                            mu = NULL, mode = "simple", longformat = F,
                            ipeffect = 4, librarysize = "increased", 
                            meth_prob = 0.25, zero_prob = NULL) {
  toycounts <- matrix(NA, nrow = N, ncol = n)
  if (mode == "simple") {}
  if (is.null(mu)) {
    mu <- rep(1:groups, each = n/groups) * 10
  }
  if (ip == T){
    inputip <- rep(c(F, T), each = n/2)
  } else if (ip == F){
    inputip <- rep(c(F), times = n)
  }
  
  if (librarysize == "increased") {
    toycounts[] <- rnbinom(n = n*N, mu = mu, size = size)
    select <- sample(c(T, F), size = N, replace = T, prob = c(meth_prob, 1- meth_prob))
    
    #methylate some counts
    toycounts[select, inputip] <- rnbinom(n = sum(select) * sum(inputip), 
                                          mu = mu * ipeffect, size = size)
    
    #make some counts zero
    if (!is.null(zero_prob)){
      select0 <- sample(c(T, F), size = (N), replace = T, prob = c(zero_prob, 1- zero_prob))
      toycounts[select0, inputip] <- 0  
    }
    
    
    # for (i in 1:n) {
    #   toycounts[, i] <- rnbinom(n = N, size = size,
    #                            mu = mu[i] * (ipeffect ^ inputip[i]))
    # }
  } else if (librarysize == "same") {
    for (j in 1:n) {
      for (i in 1:N){
        toycounts[i,j] <- rnbinom(n = 1, size = size,
                                  mu = mu[j] * ipeffect ^
                                    (inputip[j] * sample(c(-1,1),
                                                         size = 1,
                                                         prob = c(0.5, 0.5))))
      }
    }
  }
  
  rownames(toycounts) <- paste("mirna", 1:nrow(toycounts), sep = "")
  # add a mark for methylated mirnas
  rownames(toycounts)[select] <- paste(rownames(toycounts[select, ]), "_m", 
                                       sep = "")
  if (!is.null(zero_prob)){
    rownames(toycounts)[select0] <- paste(rownames(toycounts[select0, ]), "_0", 
                                          sep = "")
  }
  
  columnames <- paste("", rep(LETTERS[1:groups], each = n/groups),
                      sep = "")
  
  if (ip == T) {
    colnames(toycounts) <- paste(columnames,
                                 rep(c("_inp_", "_ip_"), each = n/2),
                                 rep(rep(1:(n/(2 * groups)), each = 2), times = groups),
                                 sep = "")
  } else {
    colnames(toycounts) <- columnames
  }
  
  if(longformat == T && ip == T) {
    toycounts <- data.frame(mirna = rep(as.factor(rownames(toycounts)), each = n),
                            batch = rep(as.factor(colnames(toycounts)), times = N),
                            group = as.factor(rep(rep(LETTERS[1:groups],
                                                      each = n/groups),
                                                  times = N)),
                            methylation = as.factor(rep(c("input", "ip"), times = n/2 * N)),
                            counts = as.vector(t(toycounts)))
  } else if(longformat == T && ip == F) {
    toycounts <- data.frame(mirna = rep(as.factor(rownames(toycounts)), each = n),
                            batch = rep(as.factor(colnames(toycounts)), times = N),
                            group = as.factor(rep(rep(LETTERS[1:groups],
                                                      each = n/groups),
                                                  times = N)),
                            methylation = as.factor(rep(c("input"), times = n * N)),
                            counts = as.vector(t(toycounts)))
  }
  
  
  return(toycounts)
}







generatemixedcounts <- function(n = 6, N = 50, phi = 1, mu1 = 20, 
                                nestedinflation = F, zeroinflation = F, mixed = F, mu2 = 5){
  m <- matrix(NA, nrow = N, ncol = n)
  m[] <- rnbinom(n = n*N, mu = mu1, size = phi)
  select <- sample(c(T,F, F, F, F), size = n * N, replace = T)
  tmp <- rnorm(n = sum(select), mean = 0, sd = 1)
  
  if (nestedinflation == T){
    m[select] <- ifelse(tmp > 0, rnbinom(n = 1, mu = mu2, size = phi), 0)
  } else if (zeroinflation == T) {
    m[select] <- 0
  } else if (mixed == T) {
    m[select] <- rnbinom(n = sum(select), mu = mu2, size = phi)
  } else {
    print("some option must be specified, dumbo!")
  }
  return(m)
}