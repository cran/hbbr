
.onLoad <- function(...) {
  if (!requireNamespace("R2jags", quietly = TRUE)) {
    stop("Package 'R2jags' needed for this package to work. Please install it.",
         call. = FALSE)
  }
}

.onAttach <- function(...) {
  packageStartupMessage(paste("Package hbbr is loaded that implements ...",
                              "\nhierarchical Bayesian benefit-risk assessment using DCE"))
}

utils::globalVariables(c("ilogit", "Del", "eps"))