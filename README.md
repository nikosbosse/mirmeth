# mirmeth

## Goal and Scope
The code allows you to analyze the methylation of miRNA (or other small RNA). It uses a Bayesian Framework that is implemented in Stan. 

## Model



## Application
as input you need a data.frame with counts like this: 

|                 | input1  | input2 | input3  | input4 | IPed1 | IPed2 | IPed3 | IPed4 |
|-----------------|---------|--------|---------|--------|-------|-------|-------|-------|
| mmu-let-7a-1    | 0       | 0      | 0       | 0      | 0     | 0     | 0     | 0     |
| mmu-let-7a-2    | 0       | 0      | 0       | 0      | 0     | 0     | 0     | 0     |
| mmu-let-7a-2-3p | 27      | 2      | 5       | 6      | 0     | 0     | 0     | 0     |
| mmu-let-7a-5p   | 1716116 | 896169 | 1045513 | 741220 | 1033  | 24128 | 25384 | 31075 |

Let's assume the counts are stored in a data.frame called counts. One should probably remove the rows where there are too few counts in the input samples. 

Call 
```R
source("../R/functions_mimdiff.R")
mirmeth_setup()
```
to load the necessary libraries and the implemented helper functions. You may have to change the path to locate the helper functions. 

By calling
```R
stanlist <- generatestanobject(counts)
```
we can then create a list of elements that can readily be plugged into Stan. If you are interested in any of the elements, e.g. the normalized counts, you can just access them via `stanlist$transformedcounts`. See `str(stanlist)` for the contents of that list. 

Let us now do the analysis: 
```R
stanfit <- rstan::stan(file = "../src/stan_files/mirmeth.stan",
                       data = stanlist,
                       iter = 5000, warmup = 1000, thin = 1, 
                       control = list(adapt_delta = 0.97))
```

The results can then most easily be accessed using
```R
res <- my_display_results(stanfit, counts)
```
this returns a list with three items: `res$pi ` gives you an estimate of the overall proportion of unmethylated miRNAs (i.e. counts don't appear under the IPed condition). `res$results` returns a table like this: 

|                 | baseMean      | se_baseMean   | delta       | se_delta     | post_prob  | meth_prob  |
|-----------------|---------------|---------------|-------------|--------------|------------|------------|
| mmu-let-7e-3p   | 4.443232e+04  | 1.190093e+04  | 7.27594898  | 1.995242329  | 0.9999375  | 0.87136811 |
| mmu-miR-142a-5p | 4.892523e+02  | 1.602365e+02  | 8.37142622  | 2.494422938  | 0.9999375  | 0.88456823 |
| mmu-miR-186-5p  | 2.257546e+02  | 7.801288e+01  | 8.36968734  | 2.553275998  | 0.9999375  | 0.88387288 |
| mmu-miR-1264-3p | 6.417675e+01  | 1.890099e+01  | 6.58779248  | 1.966591423  | 0.9998125  | 0.85835484 |


- baseMean is the same as in DESeq2, so it is the baseline mean we estimate for the input sample.

- se_baseMean is the standard error of that baseMean, i.e. an estimate for its standard deviation that tells us how uncertain our estimation of the baseMean is.

- delta is the estimated difference between the mean in input vs. ip. If delta is 5, then we assume that the mean in IP is 5 times the mean in input.

- se_delta is the standard error of delta, again an estimate of how certain we are that delta is correct.

- post_prob is the posterior probability that delta is larger than 1

- meth_prob is an estimate for the percentage of methylated miRNA similar to what was proposed in the QNB paper (Liu et al. DOI: 10.1186/s12859-017-1808-4). Take this with a grain of salt. It is calculated as baseMean * delta / (baseMean + baseMean * delta). 

`res$significant_results` just returns the miRNAs that together (i.e. multiplied) have a posterior probability of delta > 1 that is larger than 0.95. 

You can check out the file Example.R 
