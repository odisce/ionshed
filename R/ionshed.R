#' Get raw data from MSnexp object
#'
#' @inheritParams ionshed
#'
#' @import data.table
#' @importFrom rlang check_installed
#'
#' @return A `data.table` object with:
#'   - `file`: file number
#'   - `scan`: scan index
#'   - `rt`: retention time in seconds
#'   - `mz`: m/Z
#'   - `i`: Intensity
#' @export
#' @examples
#' get_rawdata(data = MSnbase::filterFile(MSnExpObj, 1))
get_rawdata <- function(data) {
  if (class(data) %in% c("XcmsExperiment", "MsExperiment")) {
    rlang::check_installed("MsExperiment")
    rlang::check_installed("Spectra")
    spctr <- MsExperiment::spectra(data)
    mz_values <- Spectra::mz(spctr)
    int_values <- Spectra::intensity(spctr)
    rt_values <- Spectra::rtime(spctr)
    file_value <- Spectra::dataOrigin(spctr)
    mslevel_value <- Spectra::msLevel(spctr)
    isowin_value <- Spectra::isolationWindowTargetMz(spctr)
  } else {
    rlang::check_installed("MSnbase")
    spctr <- MSnbase::spectra(data)
    mz_values <- lapply(spctr, MSnbase::mz)
    int_values <- lapply(spctr, MSnbase::intensity)
    rt_values <- lapply(spctr, MSnbase::rtime)
    file_value <- lapply(spctr, MSnbase::fromFile)
    mslevel_value <- lapply(spctr, MSnbase::msLevel)
    isowin_value <- lapply(spctr, function(x) {
      out <- as.numeric(NA)
      if ("Spectrum2" %in% class(x)) {
        out <- MSnbase::precursorMz(x)
      }
      return(out)
    })
  }
  output <- data.table::data.table(
    "file" = rep(file_value, lengths(mz_values)) %>% unlist(),
    "scan" = rep(seq_len(length(mz_values)), lengths(mz_values)),
    "mslevel" = rep(mslevel_value, lengths(mz_values)) %>% unlist(),
    "isowin" = rep(isowin_value, lengths(mz_values)) %>% unlist() %>% as.numeric(),
    "rt" = rep(unlist(rt_values), lengths(mz_values)),
    "mz" = unlist(mz_values),
    "i" = unlist(int_values)
  )
  data.table::setkey(output, file, scan, mslevel, rt, mz)
  return(output[])
}

#' Check available ms levels
#' @inheritParams ionshed
#' @return An integer with the available MSlevels
#' @export
check_mslevels <- function(data) {
  out <- switch(
    class(data)[1],
    "MsExperiment" = {
      rlang::check_installed("MsExperiment")
      MsExperiment::spectra(data) %>% Spectra::uniqueMsLevels()
    },
    "MSnExp" = {
      rlang::check_installed("MSnbase")
      MSnbase::msLevel(data) %>% unique()
    },
    "OnDiskMSnExp" = {
      rlang::check_installed("MSnbase")
      MSnbase::msLevel(data) %>% unique()
    },
    "data.table" = {
      if ("mslevel" %in% names(data)) {
        data[, unique(mslevel)]
      } else {
        NULL
      }
    },
    stop(sprintf("Class: %s not implemented", class(data)))
  )
  sort(out)
}

#' Extract signals coordinates and map on a sparse grid
#'
#' @param file_ind File index to subset
#' @param mz_res m/Z resolution in Dalton to create the grid
#' @param rt_res rt resolution in second to create the grid
#' @param mslevels choose msLevels to extract (default to mslevel=1)
#' @param isowin_mz `numeric` value to filter a specific isolation window
#' @param rt_range `numeric` vector `c(rtmin, rtmax)` to define the rt range to extract. Default `NULL`
#' @inheritParams ionshed
#'
#' @import data.table magrittr
#'
#' @return A `list` with:
#'   - `resolution`: a `list` with `mz` and `rt` resolutions
#'   - `data`: a `data.table` with:
#'     - `idmz`: m/Z index
#'     - `idrt`: rt index
#'     - `i`: intensity
#'     - `rt`: retention time in seconds
#'     - `mz`: m/Z
#' @export
#' @examples
#' fun_create_map_single(data = MSnbase::filterFile(MSnExpObj, 1), file_ind = 1)
fun_create_map_single <- function(
  data,
  file_ind = 1,
  mz_res = 0.00005,
  rt_res = 0.1,
  mslevels = 1,
  isowin_mz = NULL,
  rt_range = NULL
) {
  if (!is.null(rt_range)) {
    if (length(rt_range) < 2 || rt_range[2] < rt_range[1]) {
      stop("rt_range must be a numeric vector with two element: start, end)")
    }
  }
  if (!is.null(isowin_mz)) {
    if (mslevels == 1) {
      warning("isowin_mz argument can work only on MSn > 1: revert to default (NULL)")
    }
  }
  ms_dt <- switch(
    class(data)[1],
    "MsExperiment" = {
      rlang::check_installed("MsExperiment")
      mslev_sel <- intersect(
        mslevels,
        check_mslevels(data)
      )
      output <- data[file_ind] %>% MsExperiment::filterSpectra(., filter = Spectra::filterMsLevel, msLevel. = mslev_sel)
      if (!is.null(rt_range)) {
        output <- MsExperiment::filterSpectra(output, filter = Spectra::filterRt, rt = rt_range)
      }
      if (!is.null(isowin_mz)) {
        output <- MsExperiment::filterSpectra(output, filter = Spectra::filterIsolationWindow, mz = isowin_mz)
      }
      get_rawdata(output)
    },
    "MSnExp" = {
      rlang::check_installed("MSnbase")
      mslev_sel <- intersect(
        mslevels,
        check_mslevels(data)
      )
      output <- MSnbase::filterFile(data, file_ind) %>% MSnbase::filterMsLevel(., msLevel. = mslev_sel)
      if (!is.null(rt_range)) {
        output <- MSnbase::filterRt(output, rt_range)
      }
      if (!is.null(isowin_mz)) {
        output <- MSnbase::filterIsolationWindow(output, mz = isowin_mz)
      }
      get_rawdata(output)
    },
    "OnDiskMSnExp" = {
      rlang::check_installed("MSnbase")
      mslev_sel <- intersect(
        mslevels,
        check_mslevels(data)
      )
      output <- MSnbase::filterFile(data, file_ind) %>% MSnbase::filterMsLevel(., msLevel. = mslev_sel)
      if (!is.null(rt_range)) {
        output <- MSnbase::filterRt(output, rt_range)
      }
      if (!is.null(isowin_mz)) {
        output <- MSnbase::filterIsolationWindow(output, mz = isowin_mz)
      }
      get_rawdata(output)
    },
    "data.table" = {
      if (all(c("rt", "mz", "i") %in% names(data))) {
        if (!"file" %in% names(data)) {
          warning("file column not found in data, file_ind won't be used: using the whole table as one sample")
          data[, file := file_ind]
        } else {
          data <- data[file %in% file_ind, ]
        }
        if (!"isowin" %in% names(data)) {
          data[, isowin := as.integer(NA)]
        }
        if (!"mslevel" %in% names(data)) {
          data[, mslevel := as.integer(NA)]
        }
        if (!is.null(isowin_mz)) {
          data <- data[isowin %in% isowin, ]
        }

        data <- data[mslevel %in% mslevels, ]
        if (!is.null(rt_range)) {
          data <- data[rt %between% rt_range, ]
        }
        data
      } else {
        stop("missing colum in data, need:\n  mandatory: rt, mz, i\n  optional: mslevel, isowin, file")
      }
    },
    stop(sprintf("data class not recognized: %s", class(data)))
  )
  ## Generate map using mz and rt res
  ms_dt[, idmz := as.integer(mz * 1 / mz_res)]
  ms_dt[, idrt := as.integer(rt * 1 / rt_res)]
  ## Check duplicates
  dup_pos <- ms_dt[, .(valn = .N), by = .(idmz, idrt, mslevel, file, isowin)][valn > 1, ]
  if (dup_pos[, .N] > 0) {
    warning(
      sprintf(
        "%i duplicated %s with %s mz res and %s rt res: Keeping maximum intensity",
        dup_pos[, .N],
        ifelse(dup_pos[, .N] > 1, "positions", "position"),
        mz_res,
        rt_res
      )
    )
    ms_dt <- rbind(
      ms_dt[!dup_pos, on = c("idmz", "idrt", "mslevel", "file", "isowin")],
      ms_dt[dup_pos[, -c("valn")], on = c("idmz", "idrt", "mslevel", "file", "isowin")][, .SD[which.max(i)], by = .(idmz, idrt, mslevel, file, isowin)]
    )    
  }
  return(
    list(
      "resolution" = list(mz = mz_res, rt = rt_res),
      "data" = ms_dt[, .(file_nb = file_ind, mslevel, isowin, idrt, idmz, rt, mz, i)]
    )
  )
}

#' Run ionshed segmentation
#'
#' @param data OnDiskMSnExp, MSnObject, MsExperiment or a data.table with (file, rt, mz, i)
#' @param subset_dt An optional `data.table` to restrict search space with rtmin, rtmax, mzmin, mzmax
#' @param by_isowin `Logical` to run extraction for each isolation windows separately (DIA)
#' @inheritParams get_rawdata
#' @inheritParams fun_create_map_single
#' @inheritParams SpMat_DescendingROI
#' @import data.table magrittr
#'
#' @return A `list` with:
#'   - `roi_data`: A `data.table` with:
#'     - `rt_sc`: retention time index
#'     - `mz_sc`: m/Z index
#'     - `i`: intensity
#'     - `roi`: roi index
#'     - `roiType`: roi type
#'   - `peak_info`: A `data.table` with:
#'     - `roi`: roi index
#'     - `sc_length`: segment lenght
#'     - `rtid_start`: retention time index starting location
#'     - `rtid_end`: retention time index ending location
#'     - `mzid_start`: m/Z index starting location
#'     - `mzid_end`: m/Z index ending location
#'     - `rtmean`: retention time mean in seconds
#'     - `rtmin`: retention time minimum value in seconds
#'     - `rtmax`: retention time maximum value in seconds
#'     - `mzmin`: m/Z minimum value
#'     - `mzmax`: m/Z maximum value
#'     - `rtmed`: retention time median in seconds
#'     - `rtmeanw`: retention intensity weighted mean in seconds
#'     - `mzmed`: m/Z median value
#'     - `mzmeanw`: m/Z intensity weighted mean
#'     - `imax`: intensity maximum value
#'     - `imin`: intensity minimum value
#'     - `idiff`: intensity difference between max and min
#'     - `imedian`: intensity median value
#'     - `imean`: intensity mean
#'     - `into`: sum of intensities inside the segment
#'     - `mzdiff`: m/Z maximum difference
#'     - `rtdiff`: retention time maximum difference in seconds
#' @export
#' @examples
#' res <- ionshed(data = MSnExpObj, file_ind = 1, rttol = 10, ppm = 10)
#' str(res)
ionshed <- function(
  data,
  file_ind = 1,
  mz_res = 0.0001,
  rt_res = 0.01,
  rttol = 5,
  ppm = 5,
  subset_dt = NULL,
  mslevels = 1,
  rt_range = NULL,
  by_isowin = TRUE,
  debugL = FALSE
) {
  mslev_exp <- check_mslevels(data)
  if (!is.null(mslevels)) {
    mslev_sel <- intersect(mslevels, mslev_exp)
    if (any(!mslevels %in% mslev_sel)) {
      warning(
        sprintf(
          "The following asked mslevels was not found in data: %s\n Extraction done one: %s",
          paste0(mslevels[!mslevels %in% mslev_sel], collapse = "-"),
          paste0(mslev_sel, collapse = "-")
        )
      )
    }
  } else {
    mslev_sel <- mslev_exp
  }
  output <- list(
    "roi_data" = data.table(),
    "peak_info" = data.table()
  )
  for (mslev_i in mslev_sel) {
    temp_dt <- fun_create_map_single(
      data = data,
      file_ind = file_ind,
      mz_res = mz_res,
      rt_res = rt_res,
      mslevels = mslev_i,
      rt_range = rt_range
    )
    data_mslev_i <- temp_dt$data[mslevel == mslev_i, ]

    if (!is.null(subset_dt)) {
      ## get id to keep and subset them
      loc_to_keep <- c()
      for (i in subset_dt[, seq_len(.N)]) {
        rt_range <- subset_dt[i, c(rtmin, rtmax)]
        mz_range <- subset_dt[i, c(mzmin, mzmax)]
        loc_to_keep <- c(
          loc_to_keep,
          data_mslev_i[,
            which(
              rt %between% rt_range &
                mz %between% mz_range
            )
          ]
        )
      }
      if (length(unique(loc_to_keep)) <= 0) {
        warning("No remaining locations after using subset_dt")
        return(NULL)
      }
      data_mslev_i <- data_mslev_i[unique(loc_to_keep), ]
    }
    isowin_iter <- as.numeric(NA)
    if (by_isowin) {
      ## If isolation windows, run on each separately
      iso_win_vi <- data_mslev_i[!is.na(isowin), unique(isowin)]
      if (length(iso_win_vi) > 0) {
       isowin_iter <- iso_win_vi
      }
    }
    for (isowin_i in isowin_iter) {
      if (!is.na(isowin_i)) {
        data_mslev_i_isoi <- data_mslev_i[isowin %in% isowin_i, ]
      } else {
        data_mslev_i_isoi <- data_mslev_i
      }
      if (debugL == TRUE) {
        print(
          sprintf(
            "Detection in: mslevel: %i isowin: %0.1f",
            mslev_i,
            isowin_i
          )
        )
      }
      ## Get tolerance grid
      rt_grid <- ceiling(rttol / (temp_dt$resolution$rt))
      ## Set data format (mz in column for faster sparse matrix iteration)
      loc <- data_mslev_i_isoi[, .(idrt, idmz)] %>%
        as.matrix() %>%
        t()
      val <- data_mslev_i_isoi$i
      roi_table <- SpMat_DescendingROI(loc, val, rttol = rt_grid, ppm = ppm) %>%
        as.data.table()
      setnames(roi_table, c("rt_sc", "mz_sc", "i", "roi", "roiType"))
      setnames(data_mslev_i_isoi, c("idrt", "idmz"), c("rt_sc", "mz_sc"))
      roi_table <- merge(
        data_mslev_i_isoi[, -c("i")] %>% unique(),
        roi_table,
        by = c("rt_sc", "mz_sc"),
        all.y = TRUE
      )
      peak_info <- summarize_by_class(roi_table[roi != 1, ]) %>% as.data.table()
      peak_info[, mslevel := mslev_i]
      peak_info[, isowin := isowin_i]
      setkey(peak_info, roi, mslevel, isowin)
      setcolorder(peak_info)
      ## create unique integer roi ids
      roi_start <- tryCatch({output$roi_data[, max(roi)]}, error = function(e) {0})
      roi_table[, roi := roi + roi_start]
      peak_info[, roi := roi + roi_start]
      
      output$roi_data <- rbind(output$roi_data, roi_table)
      output$peak_info <- rbind(output$peak_info, peak_info)
    }

  }
  return(output)
}

#' Extract XICs
#'
#' @param roi_ls roi list obtained with `ionshed()`
#' @inheritParams SpMat_DescendingROI
#' @import data.table magrittr
#'
#' @return
#' A `data.table` object with:
#'   - `$xic_roi` `INTEGER`: roi extracted (XIC id)
#'   - `$rt` `NUMERIC`: retention time
#'   - `$mz` `NUMERIC`: m/Z
#'   - `$i` `NUMERIC`: intensity
#'   - `$roi` `INTEGER`: original roi
#'   - `$assignment` `LOGICAL`: is point inside current roi (`xic_roi`)
#' @export
extract_XICs <- function(
  roi_ls,
  rttol = 0,
  debugL = FALSE
) {
  ## get rt scale
  if (is.null(rttol)) {
    rttol <- 0
  }
  setkey(roi_ls$roi_data, 'roi')
  rttol_val <- roi_ls$roi_data[rt_sc > 1,][1, rttol / (rt / rt_sc)] %>%
    floor() %>%
    as.integer()
  out <- SpMat_getEICs(
    locations = roi_ls$roi_data[, .(rt_sc, mz_sc)] %>%
      as.matrix() %>%
      t(),
    values = roi_ls$roi_data[, .(rt, mz, i, roi)] %>%
      as.matrix(),
    eics_dt = roi_ls$peak_info[
      roi != 1,
      ][, .(roi, rtid_start, rtid_end, mzid_start, mzid_end)] %>%
      as.matrix(),
    rttol = rttol_val,
    debugL = debugL
  ) %>%
    as.data.table()
  setkey(out, 'xic_roi')
  setcolorder(out)
  return(out)
}