#' Install package environment
#' 
#' @import magrittr
init_env <- function() {
  path_to_check <- c(
    "./inst/requirements",
    ifelse(requireNamespace("ionshed", quietly = TRUE), system.file(package = "ionshed", "requirements"), "")
  )
  path_exists <- file.exists(path_to_check)
  if (all(!path_exists)) {
    stop("Can't init environment, requirement file not found on system")
  }
  packages_to_get <- readLines(path_to_check[which(path_exists)[1]]) %>% grep(pattern = "(^#|^$)", ., value = TRUE, invert = TRUE)
  ## if renv not present, install
  rlang::check_installed("renv")
  renv::restore()
  renv::install(packages_to_get)
}