skip_if_not_installed("MSnbase")

testthat::test_that("ion-shed on MSnExperiments", {
  ion_res <- ionshed(
    data = msnexp_obj,
    file_ind = 1,
    rttol = 5,
    ppm = 3,
    mz_res = 0.0001,
    rt_res = 0.1
  )
  testthat::expect_true("list" %in% class(ion_res))
  testthat::expect_true(length(ion_res) == 2)
  testthat::expect_true(all(c("roi_data", "peak_info") %in% names(ion_res)))
  testthat::expect_true(all(sapply(ion_res, is.data.table)))
  testthat::expect_true(ion_res$peak_info[duplicated(roi), ][, .N] == 0)
  testthat::expect_true(
    all(
      ion_res$roi_data[, unique(rt)] %between%
        range(MSnbase::rtime(MSnbase::filterFile(msnexp_obj, 1)))
    )
  )
  testthat::expect_true(
    all(
      ion_res$peak_info[, unique(c(mzmed, mzmin, mzmax))] %between%
        range(MSnbase::mz(MSnbase::filterFile(msnexp_obj, 1)))
    )
  )
  testthat::expect_true(
    ion_res$roi_data[, .N] <= {
      MSnbase::filterFile(msnexp_obj, 1) |>
        MSnbase::spectra() |>
        lapply(MSnbase::intensity) |>
        unlist() |>
        length()
    }
  )
})

test_that("xic extraction with sparse matrix", {
  ion_res <- ionshed(
    data = msnexp_obj,
    file_ind = 1,
    rttol = 5,
    ppm = 3,
    mz_res = 0.0001,
    rt_res = 0.1
  )
  for (rttol in c(0, 10, 100)) {
    res <- extract_xics(
      roi_ls = ion_res,
      debug = FALSE,
      rttol = rttol
    )
    temp_dt <- merge(
      ion_res$peak_info[, .(roi, sc_length, rtmin, rtmax)],
      res[,
        .(
          scan_nb = .N,
          rtmin_xic = min(rt),
          rtmax_xic = max(rt)
        ),
        by = .(roi = xic_roi)
      ],
      by = "roi"
    )
    expect_true(temp_dt[, all(scan_nb >= sc_length), ])
    expect_true(all(ion_res$peak_info[, unique(roi)] %in% res[, unique(roi)]))
    expect_true(
      all(
        temp_dt[, {
          out <- all(
            rtmin_xic >= rtmin - rttol,
            rtmax_xic <= rtmax + rttol
          )
          .(test = out)
        }, by = .(roi)][, test]
      )
    )
  }
})


testthat::test_that("From .mzML with MSnbase, MsExperiment and data.table", {
  testthat::skip_if_not_installed("MSnbase")
  mzml_path <- sapply(
    c(
      system.file("test.mzml", package = "ionshed"),
      "./inst/test.mzml"
    ),
    file.exists
  )
  if (all(mzml_path == FALSE)) {
    testthat::skip("test.mzml file not found")
  } else {
    mzml_path <- names(mzml_path[mzml_path == TRUE][1])
    for (msleveli in c(1, 2)) {
      for (format_i in c("MSnbase", "MsExperiment", "table")) {
        msnexp <- switch(
          format_i,
          "MSnbase" = MSnbase::readMSData(mzml_path, msLevel. = msleveli),
          "MsExperiment" = MsExperiment::readMsExperiment(mzml_path),
          "table" = get_rawdata(
            MsExperiment::readMsExperiment(mzml_path)
          )[, .(file = 1, mslevel, rt, mz, i)]
        )
        rt_range <- c(40, 45)
        ion_res <- ionshed(
          data = msnexp,
          mz_res = 0.0001,
          rt_res = 0.01,
          rttol = 5,
          ppm = 5,
          mslevels = msleveli,
          rt_range = rt_range
        )
        testthat::expect_true("list" %in% class(ion_res))
        testthat::expect_true(length(ion_res) == 2)
        testthat::expect_true(all(c("roi_data", "peak_info") %in% names(ion_res)))
        testthat::expect_true(all(sapply(ion_res, is.data.table)))
        testthat::expect_true(ion_res$peak_info[duplicated(roi), ][, .N] == 0)
        testthat::expect_true(ion_res$peak_info[, all(rtmin >= min(rt_range))])
        testthat::expect_true(ion_res$peak_info[, all(rtmax <= max(rt_range))])
        testthat::expect_true(ion_res$roi_data[, all(rt >= min(rt_range))])
        testthat::expect_true(ion_res$roi_data[, all(rt <= max(rt_range))])
        testthat::expect_true(ion_res$roi_data[, all(mslevel == msleveli)])
      }
    }
  }
})
