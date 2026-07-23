#' Create XCMS object with ionshed peaks
#'
#' @param peak_info peak info obtained from `ionshed()`
#' @inheritParams SpMat_DescendingROI
#' @import data.table magrittr
#' @return
#' A `xcms object` with manual peak integration.
#' @export
get_xcmsobj <- function(peak_info, msnobj, msLevel = 1L) {
  if (!requireNamespace("xcms", quietly = TRUE)) {
    stop("xcms package needs to be installed to use this function")
  }
  obj_out <- xcms::manualChromPeaks(
    msnobj,
    chromPeaks = peak_info[
      ,
      .(
        mzmin = mzmin - 0.00005,
        mzmax = mzmax + 0.00005,
        rtmin = rtmin - 0.001,
        rtmax = rtmax + 0.001
      )
    ] %>%
      as.matrix(),
    msLevel = msLevel
  )
  obj_out@.processHistory <- list(xcms:::ProcessHistory())
  return(obj_out)
}