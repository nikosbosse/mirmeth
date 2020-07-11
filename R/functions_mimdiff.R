#' @title Normalize Genes According to GMPR
#'
#' @description
#' 1st step: calculate pairwise r_jk for two samples j and k
#' only look at those genes in the samples for which
#' count_ij * count_ik != zero
#' for those were it is not zero, look at the ratio
#' count_ij / count_ik. Then take the median of those.
#' 
#' @details
#' 
#' @param counts A data.frame with counts
#' 
#' @return MISSING
#' @author Nikos Bosse \email{nikosbosse@gmail.com}
#' @examples
#'  
#' counts <- mirmeth::example_data_wide
#' s_gmpr(counts)
#'
#' @export


gmpr_normalization <- function(counts) {
  counts <- counts[, !(colnames(counts) == "id")]
  n = ncol(counts)
  s_gmpr <- vector(mode = 'numeric', length = n)
  for (j in 1:n) {
    r <- vector(mode = 'numeric', length = n)
    for (k in 1:n) {
      v <- counts[,j]
      w <- counts[,k]
      nonzero <- v != 0 & w
      r[k] <- median(v[nonzero] / w[nonzero])
    }
    s_gmpr[j] <- exp((1/n) * sum(log(r)))
  }
  return(s_gmpr)
}

#===============================================================================


#' @title Generate Stanfit Objcet from Counts
#'
#' @description
#' 
#' @details
#' 
#' @param counts A data.frame with counts
#' 
#' @return MISSING
#' @author Nikos Bosse \email{nikosbosse@gmail.com}
#' @examples
#' # 
#' # missing
#'
#'
#' @export


generatestanobject <- function(counts, s_to_one = F, g_to_one = T){

  n <- ncol(counts)
  N <- nrow(counts)
  counts <- as.matrix(counts)
  ## Error handling. 
  ## warnings if some miRNAs are zero across all samples. 
  ## will be repeated below. Needs to be changed for efficiency reasons. 
  meanmatrix <- matrix(NA, nrow = N, ncol = 3)
  colnames(meanmatrix) <- c("input", "ip", "delta")
  meanmatrix[,1] <- rowMeans(counts[, 1:(n/2)])
  meanmatrix[,2] <- rowMeans(counts[, (1 + (n/2)):n])
  meanmatrix[,3] <- meanmatrix[,2] / meanmatrix[,1]

  if ( sum(meanmatrix[,1] ==  0) > 0 | sum(meanmatrix[,2] ==  0) < 0 ){
    warning("Some miRNAs had zero counts across all input samples and/or across all IPed samples. These were removed.")
      counts <- counts[(meanmatrix[,1] != 0 & meanmatrix[,2] != 0),]
      N = nrow(counts)

  } 


  #vector indicating whether we have IP or not
  ip <- rep(c(0,1), each = n/2)

  #compute size factor for the samples
  if (s_to_one == T){
    s <- rep(1, times = n)
  } else {
    s  <- vector(mode = "numeric", length = n)
    for (i in 1: n){
      s[i] <- sum(counts[,i])/sum(counts)
    }
  }

  #compute scaling factor for the bins w
  if (g_to_one == T){
    g <- rep(1, times = N)
  } else {
    g <- vector(mode = "numeric", length = N)
    g <- rowMeans(counts[,ip == 0])
  }

  ##############################################
  # alternative: DESeq2 normalization approach
  rowprodsfun <- function(x) {
    prod(x[x !=0]) ^ (1/n) # over sum(ip)?
  }
  Ki_R <- apply(counts[, ip == 0], MARGIN = 1, FUN = rowprodsfun)
  # Ki_R <- matrixStats::rowProds(counts[, ip == 0])
  s_deseq <- counts[] / Ki_R # counts[, ip == 0]?
  s_deseq <- matrixStats::colMedians(as.matrix(s_deseq)) 
  ##############################################

  ##############################################
  # alternative: GMPR normalization approach
  # 1st step: calculate pairwise r_jk for two samples j and k
  # only look at those genes in the samples for which
  # count_ij * count_ik != zero
  # for those were it is not zero, look at the ratio 
  # count_ij / count_ik. Then take the median of those. 

  # s_gmpr <- vector(mode = 'numeric', length = n)
  # for (j in 1:n) {
  #   r <- vector(mode = 'numeric', length = n)
  #   for (k in 1:n) {
  #     v <- counts[,j]
  #     w <- counts[,k]
  #     nonzero <- v != 0 & w
  #     r[k] <- median(v[nonzero] / w[nonzero])
  #   }
  #   s_gmpr[j] <- exp((1/n) * sum(log(r)))
  # }

  s_gmpr <- gmpr_normalization(counts)

  ##############################################


  ##############################################
  # Calulculations for empirical Bayes Priors
  transformedcounts <- counts / g
  for (i in 1:nrow(counts)) {
    transformedcounts[i,] <- transformedcounts[i,] / s_gmpr
  }

  ## means
  ## for colmeans, take out the ones that are zero, 
  ## because they are looked at independently anyway

  # means <- rep(0, n)
  # means[1:(n/2)] <- colMeans(transformedcounts[, 1:(n/2)])
  # nonzero <- !(rowMeans(transformedcounts[, (1 + n/2):n]) == 0)
  # means[(1 + n/2):n] <- colMeans(transformedcounts[nonzero, (1 + n/2):n])
  # deltas <-means[(1 + n/2):n] /  means[1:(n/2)]

  meanmatrix <- matrix(NA, nrow = N, ncol = 3)
  colnames(meanmatrix) <- c("input", "ip", "delta")
  #only use the input values to compute base-mean
  meanmatrix[,1] <- rowMeans(transformedcounts[, 1:(n/2)])
  #only use the ip values to compute base-mean of ip
  meanmatrix[,2] <- rowMeans(transformedcounts[, (1 + (n/2)):n])
  #calculate the deltas
  meanmatrix[,3] <- meanmatrix[,2] / meanmatrix[,1]

  mean <- mean(meanmatrix[,1])
  varmean <- var(meanmatrix[,1])

  ## mean and variance of delta
  meandelta <- mean(meanmatrix[meanmatrix[,3] != 0, 3])
  vardelta <- var((meanmatrix[meanmatrix[,3] != 0, 3]))


  ## calculate appropriate values for gamma priors.
  ## Set mu = empirical mean and var = 10 * empirical variance
  ## mean = alpha / beta, var = alpha / beta^2 
  ## beta = mu / var
  ## alpha = mu / beta
  beta_mu <- mean / (10 * varmean)
  alpha_mu <- mean * beta_mu
  
  beta_delta <- meandelta / (10 * vardelta)
  alpha_delta <- meandelta * beta_delta




  # X <- transformedcounts[rowSums(transformedcounts) != 0,]
  # MASS::glm.nb(data=X[1,], formula=~ -1 + group + ip)
  # group <- rep(1:(n/2), times = 2)


  # group <- factor(rep(1:4, 2))
  # ip <- factor(rep(0:1, each = 4))

  # d <- as.data.frame(cbind(group, ip))
  # X <- model.matrix(~ -1 + as.factor(group) + as.factor(ip), data = d)
  ## https://stats.stackexchange.com/questions/11182/when-to-use-an-offset-in-a-poisson-regression

  #######################################################
  ## calculate preliminary dispersion estimates in order to create
  ## an empirical Bayes Prior for the dispersion estimate

  ## https://stats.stackexchange.com/questions/11182/when-to-use-an-offset-in-a-poisson-regression 
  ## here we only use the data in input. We could also use the data
  ## in ip if we replaced exp(basemean) by exp(basemean %*% ip)
  ## unsure whether the betas in xi'beta also include the groups or 
  ## that effect is in theory accounted for by normalization. 

  neg_ll <- function(params, countmatrix, s, gg) {
    l = 0.0
    N = nrow(countmatrix)
    n = ncol(countmatrix)
    basemean <- params[1:(N)] 
    dispersion <- params[(1 + N):(2 * N)] 
    for (i in 1:N) {
      for (j in 1:(n/2)) {
        offset <- s[j] * gg[i]
        #if (countmatrix[i,j] == 0) next
        l = l - dnbinom(countmatrix[i,j], mu = exp(basemean[i]) * offset, size = dispersion[i], log = T)
      }
    }
    return(l)
  }

  # pars <- optim(par = rep(1, (2 * N)), fn = neg_ll, s = s, gg = g, countmatrix = counts)$par
  # dispersions <- pars[(length(pars) / 2 + 1):(length(pars))]

  # mean_disp <- mean(dispersions)
  # var_disp <- var(dispersions)
  # beta_phi <- mean_disp / (10 * var_disp)
  # alpha_phi <- mean_disp * beta_phi

  beta_phi <- 0.1
  alpha_phi <- 0.1
  if (beta_phi == Inf) beta_phi <- 0.1
  if (alpha_phi == Inf) alpha_phi <- 0.1

  #######################################################
  

  ##############################################
  return(list(y = counts,
              transformedcounts = transformedcounts,
                   N = N,
                   n = n,
                   s = s,
                   g = g,
                   ip = ip,
                   mean = mean,
                   varmean = varmean,
                   delta = meandelta,
                   vardelta = vardelta,
                   beta_mu = beta_mu,
                   alpha_mu = alpha_mu, 
                   beta_delta = beta_delta,
                   alpha_delta = alpha_delta,
                   beta_phi = beta_phi,
                   alpha_phi = alpha_phi,                    
                   s_deseq = s_deseq, 
                   s_gmpr = s_gmpr))

}

#===============================================================================

getmlestimates <- function(counts, summary = F){
  n <- ncol(counts)
  N <- nrow(counts)
  mlest <- matrix(nrow = N, ncol = 2)
  colnames(mlest) <- c("phi", "mu")
  for (i in 1:N){
    mlest[i,] <- MASS::fitdistr(counts[i,], "Negative Binomial")$estimate
  }

  if (summary == T){
    return(list(estimates = mlest,
                colMeans = colMeans(mlest),
                medianphi = median(mlest[,1]),
                medianmean = median(mlest[,2]),
                perc_phiestgreaterphi = sum(mlest[,1] >= 2)/nrow(mlest)))
  }
  else {
    return(mlest)
  }
}

#===============================================================================
extractparameters <- function(stanfitobject){
  temp <- rstan::extract(stanfitobject)
  parameters <- list()
  for (i in 1:(length(temp) - 1)){
    parameters[[i]] <- colMeans(temp[[i]])
  }

  parameternames <- names(temp)[-length(temp)] #without the last item, log likelihood
  names(parameters) <- parameternames

  #compute means for all
  parameters$mean_estimates <- vector(length = length(parameternames))
  for (i in 1:(length(parameternames))){
    parameters$mean_estimates[i] <- mean(parameters[[i]])
  }
  #name vector items of last list item
  names(parameters$mean_estimates) <- parameternames

  #compute medians for all

  parameters$median_estimates <- vector(length = length(parameternames))
  for (i in 1:(length(parameternames))){
    parameters$median_estimates[i] <- median(parameters[[i]])
  }
  #name vector items of last list item
  names(parameters$median_estimates) <- parameternames

  return(parameters)
}

#===============================================================================

comparefits <- function(stanfitobject, countsdata, mode = "both", log = F){
  fits <- list()
  bayes <- extractparameters(stanfitobject)
  if (log == T) {
    for (i in 1:(length(bayes)-1))
    bayes[[i]] <- exp(bayes[[i]])
  }

  ml <- getmlestimates(countsdata)

  if (mode == "mean"){
    fits$bayes_mean <- bayes$mean_estimates
    fits$maxlik_mean <- colMeans(ml)
    fits$data <- colMeans(countsdata)

  } else if (mode == "median"){
    fits$bayes_median <- bayes$median_estimates
    fits$maxlik_median <- robustbase::colMedians(ml)
    fits$data <- robustbase::colMedians(countsdata)

  } else if (mode == "both"){
    fits$bayes_mean <- bayes$mean_estimates
    fits$bayes_median <- bayes$median_estimates
    fits$maxlik_mean <- colMeans(ml)
    fits$maxlik_median <- robustbase::colMedians(ml)
    fits$data_mean <- colMeans(countsdata)
    fits$data_median <- robustbase::colMedians(countsdata)
  }
  return(fits)
}

#===============================================================================

get_predicted_group_means <- function(stanfitobject, standata, d0w = T, scalefactors = T){
  n <- standata$n
  N <- standata$N
  params <- extractparameters(stanfitobject = stanfitobject)
  meanmatrix <- matrix(NA, nrow = N, ncol = n)

  if (scalefactors == F){
    #meanmatrix[i,j] <- mu[i] * delta[i] ^ ip[j]
  } else if(d0w == T){
    for (j in 1:n){
      for (i in 1:N){
        meanmatrix[i,j] <- standata$s[j] * standata$g[i] * standata$d0w[i] * params$delta[i] ^ standata$ip[j]
      }
    }
  } else if (d0w == F){
    meanmatrix <- matrix(NA, nrow = N, ncol = n)
    for (j in 1:n){
      for (i in 1:N){
        meanmatrix[i,j] <- standata$s[j] * standata$g[i] * params$delta[i] ^ standata$ip[j]
        #meanmatrix[i,j] <- d0w[i] * delta_w[i] ^ ip[j]
      }
    }
  }
  colnames(meanmatrix) <- paste("means_group_", 1:n)
  return(meanmatrix)
}

#===============================================================================
makecoldata <- function(counts, origin = NULL, batch = c(1,1,1,1,2,2,2,2)){
  if (origin == '2019') {  
  coldata <- data.frame(methylation = ifelse(grepl("input", colnames(counts)),
                                       "input", "iped"),
                        condition = ifelse(grepl("Old", colnames(counts)),
                                       "old", "young"),
                        #condition = gsub("[^ABCDE]","",colnames(counts)),
                        id = as.factor(gsub(".*?([0-9]+).*", "\\1", colnames(counts))),
                        batch = paste("batch", batch), sep = "")
  rownames(coldata) <- colnames(counts)
  return(coldata)
  }

  if (origin == '2018') {  
  coldata <- data.frame(methylation = ifelse(grepl("input", colnames(counts)),
                                       "input", "iped"),
                        condition = ifelse(grepl("Old", colnames(counts)),
                                       "old", "young"),
                        #condition = gsub("[^ABCDE]","",colnames(counts)),
                        id = as.factor(gsub(".*?([0-9]+).*", "\\1", colnames(counts))))
                        #batch = paste("batch", batch), sep = ""
  rownames(coldata) <- colnames(counts)
  return(coldata)
  }

  coldata <- data.frame(methylation = ifelse(grepl("input", colnames(counts)),
                                       "input", "iped"),
                        condition = gsub("[^ABCDE]","",colnames(counts)),
                        id = as.factor(gsub(".*?([0-9]+).*", "\\1", colnames(counts))),
                        batch = paste("batch", c(1,1,1,2,1,2,1,2,1,2,2)), sep = "")
  rownames(coldata) <- colnames(counts)
  return(coldata)
}

#===============================================================================
posteriordensities <- function(stanfit, counts, draws = 50, genes = 10){

  y_post_pred <- rstan::extract(stanfit)$y_pred

  # === density for real data ==== #
  if (genes == "all"){
    genes <- 1:nrow(counts)
  }
  else if (genes == 1) {
    genes <- sample(x = 1:nrow(counts), size = genes, replace = F)
      dens <- density(counts[genes,], from = 0)
      plot(NA, ylim = range(dens$y), xlim = range(dens$x),
           ylab = "density",
           xlab = "counts")
      color_real <- rainbow(length(genes))
      color_pred <- colorspace::lighten(color_real, amount = 0.7)
      rand <- sample(x = 1:dim(y_post_pred)[1], draws, replace = F)
      for (i in rand){
        dens2 <- density(y_post_pred[i,genes,], from = 0)
        lines(dens2, col=color_pred)
        #mapply(lines, dens2, col=color_pred)
      }
      lines(dens, col=color_real, lwd = 2)
      legend("topright", legend= paste("mirna", genes), fill=color_real)
      return()
  }
  else {
    genes <- sample(x = 1:nrow(counts), size = genes, replace = F)
  }
  dens <- apply(counts[genes,], 1, density, from = 0)
  plot(NA, xlim=range(sapply(dens, "[", "x")),
       ylim=range(sapply(dens, "[", "y")),
       ylab = "density",
       xlab = "counts")

  # ==== make suitable colors === #
  color_real <- rainbow(length(genes))
  color_pred <- colorspace::lighten(color_real, amount = 0.7)
  #color_real <- colorspace::darken(color_real, amount = 0.3)

  # === density for posterior predictive data ==== #
  rand <- sample(x = 1:dim(y_post_pred)[1], draws, replace = F)
  for (i in rand){
    dens2 <- apply(y_post_pred[i,genes,], 1, density, from = 0)
    mapply(lines, dens2, col=color_pred)
  }
  invisible(mapply(lines, dens, col=color_real, lwd = 2)) #avoids printing unnecessary output
  legend("topright", legend=names(dens), fill=color_real)
}


#===============================================================================
makepvalueplot <- function(res_not_na){
  use <- res_not_na$baseMean > metadata(res_not_na)$filterThreshold
  h1 <- hist(res_not_na$pvalue[!use], breaks=0:50/50, plot=FALSE)
  h2 <- hist(res_not_na$pvalue[use], breaks=0:50/50, plot=FALSE)
  colori <- c(`do not pass`="khaki", `pass`="powderblue")
  barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
          col = colori, space = 0, main = "", ylab="frequency")
  text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
       adj = c(0.5,1.7), xpd=NA)
  legend("topright", fill=rev(colori), legend=rev(names(colori)))
  return()
}


#===============================================================================
getmapestimates <- function(stanfile, standata, seed = 1,
                            pars = c("phi", "mu1", "mu2")){

  parameters <- vector("list", length = length(pars))
  c <- optimizing(stan_model(file = stanfile),
                  data = standata,
                  iter = 1000, seed = seed)
  for (i in 1:length(pars)){
    parameters[[i]] <- c$par[grepl(pars[i], names(c$par))
                             & !grepl("log", names(c$par))
                             & !grepl("reciprocal", names(c$par))]
  }
  names(parameters) <- pars
  map <- rep(NA, times = length(parameters) * 2)
  for (i in 1:length(parameters)){
    map[2*i - 1] <- median(parameters[[i]])
    map[2*i] <- mean(parameters[[i]])
  }
  names(map) <- c(paste(c("median_", "mean_"), rep(names(parameters), each = 2), sep = ""))
  return(map)
}

#===============================================================================
getoverlapping <- function(allgenes, genelist1, genelist2, genelist3, genelist4){
  allmirnas <- rownames(allgenes)
  overlap <- (allmirnas %in% rownames(genelist1) & allmirnas %in% rownames(genelist2) & allmirnas %in% rownames(genelist3) & allmirnas %in% rownames(genelist4))

  return(allmirnas[overlap])
}


#===============================================================================

parsfoldednormal <- function(mu, sd) {

    mu_true <- sd * sqrt(2/pi) * exp(-mu^2 / (2 * sd^2)) + mu * (1 - 2 * pnorm((- mu / sd), mean = 0, sd = 0)) 

    true_var <- mu^2 + sd^2 - mu_true^2  

    return(list(mu_true = mu_true, sd_true = sqrt(true_var)))   
}


#===============================================================================

prop_zeroinfl <- function(countmatrix){
    n <- ncol(countmatrix)
    zero <- countmatrix[, (n/2 + 1):n] == 0
    return(sum(zero) / length(zero))
}

prob_zeroinfl <- function(countmatrix){
    n <- ncol(countmatrix)
    zero <- countmatrix[, (n/2 + 1):n] == 0
    return(sum(zero) / length(zero))
}

#===============================================================================


#===============================================================================
# basic summary statistic and plotting functions
#===============================================================================

make_library_size_plot <- function(counts) {
  librarysizes <- colSums(counts)
  color <- grepl("IP|Iped|ped", colnames(counts)) + 1
  a <- barplot(colSums(counts), log = "y", col = color, #ylim = c(0, 30000000),
        ylab = "Total Read Counts (Log scale)", main = "Library Sizes")
  text(a, 1.3 * librarysizes, labels=librarysizes, xpd=TRUE)
  return('plot generated')
}

#===============================================================================

make_nonzero_plot <- function(counts) {
  nonzeros <- rep(NA, length(counts))
  for (i in 1:length(counts)){
    names <- colnames(counts)
    nonzeros[i] <- sum(counts[,i] != 0)
  }
  a <- barplot(nonzeros, col = (grepl("ip|IP", colnames(counts))+ 1),
          ylab = "nonzero reads (samples combined)", main = "Number of miRNA detected (non-zero read counts)")
  text(a, -50, labels=nonzeros, xpd=TRUE)
  text(a, -100, labels=colnames(counts), xpd=TRUE)
  return('plot generated')  
}

#===============================================================================

my_make_raw_counts_plot <- function(counts, IP = TRUE, mfrow = c(2,2), zeros = 'included') {
  par(mfrow = mfrow)
  if (IP == TRUE) {
    ipeds <- colnames(counts)[grepl("IP|ip|ped", colnames(counts))]
    for (name in ipeds){
      if (zeros == 'included') {
        a <- hist(log10(counts[,colnames(counts) == name]+1), breaks = 100, plot = F)
        a$counts <- a$counts / sum(a$counts)
        plot(a, main = paste(c("Histogram of log10 counts of"), name),
             xlab = c("log10 counts"))
      } else {
        t <- counts[,colnames(counts) == name]
        a <- hist(log10(t[t != 0]), breaks = 100, plot = F)
        a$counts <- a$counts / sum(a$counts)
        plot(a, main = paste(c("Histogram of log10 counts (without zeros) of"), name),
          xlab = c("log10 counts"))
      }
    }
    par(mfrow = c(1,1))
    return('IP Plot Generated')
  } else {
    inputs <- colnames(counts)[grepl("input", colnames(counts))]
    for (name in inputs){
      if (zeros == 'included') {
        a <- hist(log10(counts[,colnames(counts) == name]+1), breaks = 100, plot = F)
        a$counts <- a$counts / sum(a$counts)
        plot(a, main = paste(c("Histogram of log10 counts of"), name),
             xlab = c("log10 counts"))
      } else {
        t <- counts[,colnames(counts) == name]
        a <- hist(log10(t[t != 0]), breaks = 100, plot = F)
        a$counts <- a$counts / sum(a$counts)
        plot(a, main = paste(c("Histogram of log10 counts (without zeros) of"), name),
          xlab = c("log10 counts"))
      }
    } 
    par(mfrow = c(1,1))
    return('Input Plot Generated')
  }
}

#===============================================================================

RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}


#===============================================================================


my_posterior_probs <- function(stanfitobject, cutoff = 1) {
  mus <- extract(stanfitobject)$mu
  deltas <- extract(stanfitobject)$delta

  len = nrow(mus)
  postprob <- apply(deltas, MARGIN = 2, FUN = function(x) {sum(x > cutoff) / len})
  #postmedian <- apply(mus, MARGIN = 2, FUN = function(x) {sum(x > 1) / len})

  # ls <- list(post_prob <- postbrob, 
  #   something = something)
  return(postprob)
}

#===============================================================================

my_methylation_prob_qnb <- function(stanfitobject){
  mus <- extract(stanfitobject)$mu
  deltas <- extract(stanfitobject)$delta

  ipreads <- mus * deltas
  methprob <- ipreads / (mus + ipreads) 
  methprob <- apply(methprob, MARGIN = 2, FUN = mean)
  return(methprob)
}


#===============================================================================

my_display_results <- function(stanfitobject, counts, cutoff = 1) {

  summary <- summary(stanfitobject)$summary
  ## extract pi
  pi <- summary[1,]#summary[grepl('pi', rownames(summary)), ]

  ## extract mus and deltas
  mus <- summary[grepl('mu', rownames(summary)), ]
  deltas <- summary[grepl('delta', rownames(summary)), ]

  ## make data frame
  results <- data.frame(row.names = rownames(counts))
  results$baseMean <- mus[,1]
  results$se_baseMean <- mus[,3] 
  results$delta <- deltas[,1]
  results$se_delta <- deltas[,3]
  results$post_prob <- my_posterior_probs(stanfitobject, cutoff)
  results$meth_prob <- my_methylation_prob_qnb(stanfitobject)
  results <- results[order(results$post_prob, decreasing = T), ]

  ## calculate significant top results
  p = 1
  i = 1
  ind = 0
  while(TRUE) {
    p <- p * results$post_prob[i]
    i <- i + 1
    if (p <= 0.95) break
    ind = ind + 1
  }
  significant_results <- results[1:ind, ]
  return(list(pi = pi, significant_results = significant_results, results = results))
}

#===============================================================================

mirmeth_setup <- function(){
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
}

