# ================ setup ================ #
setwd("/home/nikos/Desktop/Doktorarbeit/mirmeth/R")
source("../R/functions_mimdiff.R")
mirmeth_setup()

# ================ data ================= #
## read data (and set correct row names if necessary)
counts <- read.csv("../data/sample_counts.csv", header = T)
rownames(counts) <- counts[, 1]
counts <- counts[, -1]     

# counts[nrow(counts) + 1, ] 

# ============== analysis =============== #
stanlist <- generatestanobject(counts)

# look at normalized counts if you like
stanlist$transformedcounts

# stanfit <- rstan::stan(file = "../src/stan_files/mirmeth_improved.stan",
#                        data = stanlist,
#                        iter = 5000, warmup = 1000, thin = 1, 
#                        control = list(adapt_delta = 0.97))

# res <- my_display_results(stanfit, counts)