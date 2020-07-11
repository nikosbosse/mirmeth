#' Example Data Wide Format
#'
#' A toy dataset for a probability forecast of a binary outcome variable
#'
#' @format A data.frame with 31 rows and 9 variables:
#' \describe{
#'   \item{id}{unique identifier for true observed values}
#'   \item{input_1}{counts in input sample 1}
#'   \item{input_2}{counts in input sample 2}
#'   \item{input_3}{counts in input sample 3}
#'   \item{input_4}{counts in input sample 4}
#'   \item{iped_1}{counts in immunoprecipitated sample 1}
#'   \item{iped_2}{counts in immunoprecipitated sample 2}
#'   \item{iped_3}{counts in immunoprecipitated sample 3}
#'   \item{iped_4}{counts in immunoprecipitated sample 4}
#' }

"example_data_wide"




#' Example Data Long Format
#'
#' A toy dataset for a probability forecast of a binary outcome variable
#'
#' @format A data.frame with 248 rows and 4 variables:
#' \describe{
#'   \item{id}{unique identifier for true observed values}
#'   \item{count}{counts observed in corresponding sample}
#'   \item{sample_type}{input or immunoprecipitated sample}
#'   \item{sample_nr}{sample number that denotes corresponding pairs of input
#'   and iped sample}
#' }

"example_data_long"

