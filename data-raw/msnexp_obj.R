## code to prepare `msnexp_obj` dataset goes here

mzml_files <- list.files(
  "barb_data/",
  pattern = ".*10ng.*\\.mzML",
  full.names = TRUE
)[1:3]
data_in <- MSnbase::readMSData(
  files = mzml_files,
  mode = "inMemory",
  msLevel. = 1L
)
msnexp_obj <- data_in |>
  MSnbase::filterMz(c(90, 160)) |>
  MSnbase::filterRt(c(10, 60)) |>
  MSnbase::filterZero()
msnexp_obj |> object.size() |> format("Kb")
usethis::use_data(msnexp_obj, overwrite = TRUE)
