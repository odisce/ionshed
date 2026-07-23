## code to prepare `MSnExpObj` dataset goes here

mzml_files <- list.files(
  "/spi/scidospace/data/mzml/",
  pattern = ".mzML",
  full.names = TRUE
)[1:3]

data_in <- MSnbase::readMSData(files = mzml_files, mode = "inMemory", msLevel. = 1L)
require(magrittr)
MSnExpObj <- data_in %>%
  filterMz(., c(100,120)) %>%
  filterRt(., c(10,60))
MSnExpObj %>% object.size() %>% format("Kb")
usethis::use_data(MSnExpObj, overwrite = TRUE)
