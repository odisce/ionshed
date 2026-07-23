#' Create XCMS object with ionshed peaks
#'
#' @param peak_info peak info obtained from `ionshed()`
#' @param msLevel MS level to assign peaks in the xcms object
#' @inheritParams SpMat_DescendingROI
#' @inheritParams ionshed
#' @import data.table magrittr
#' @return
#' A `xcms object` with manual peak integration.
#' @export
get_xcmsobj <- function(peak_info, data, msLevel = 1L) {
  rlang::check_installed("xcms")
  obj_out <- xcms::manualChromPeaks(
    data,
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
  return(obj_out)
}