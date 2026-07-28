<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/odisce/ionshed/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/odisce/ionshed/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/odisce/ionshed/graph/badge.svg)](https://app.codecov.io/gh/odisce/ionshed)
<!-- badges: end -->

# **ionshed**

**ionshed** is an R package for **LC-MS raw data segmentation**, inspired by the Watershed algorithm and based on sparse matrix segmentation. Originally designed for LC-MS data (MS1 and DIA), it parses 3D coordinates (retention time, m/z, intensity) by decreasing intensities and clusters data points based on their proximity in the 2 dimensions (retention time and m/z).

The algorithm can also be adapted to other 3D datasets by providing coordinates in a table format (see [Usage](#usage)).

### **Key Parameters**
- **`rttol`**: Chromatographic resolution (in seconds). Defines the minimal distance between two signals of the similar intensities to separate chromatographic peaks **at a given intensity level**.
- **`ppm`**: Mass resolution (in Daltons). Defines the distance between two signals to differentiate specific ions.


## **Requirements**
- **R >= 4.1.0**
- **[Rtools](https://cloud.r-project.org/)** (to compile from source)


## **Installation**

| method | command |
|--|--|
| [renv](https://cran.r-project.org/web/packages/renv/) | `renv::install("github::odisce/ionshed")` |
| [pak](https://cran.r-project.org/web/packages/pak/) | `pak::pkg_install("odisce/ionshed")` |

## Usage

### From mz(X)ML data

Load `.mzML` files using one of the following methods:

| format | command |
| :--- | :--- |
| `MSnbase` | `msn_in <- MSnbase::readMSData("path/to/my/.mzml", mode = "onDisk", msLevel. = 1L)` |
| `MsExperiment` | `msn_in <- MsExperiment::readMsExperiment(spectraFiles = "path/to/my/.mzml")` |


Run the `ionshed` algorithm on the first file (`file_ind = 1`):

```r
ion_res <- ionshed::ionshed(
  data = msn_in,
  file_ind = 1,
  rttol = 5,
  ppm = 3
)
```

The results are stored in a list with two levels:
  - `ion_res["roi_data"]`: Original segmented coordinates
    ```r
    str(ion_res$roi_data)
    ```
  - `ion_res["peak_info"]`: Summary of each segment
    ```r
    str(ion_res$peak_info)
    ```


### From a table

To perform the segmentation on a table, provide a table with the following mandatory columns:
  - `rt`: Retention time values (seconds)
  - `mz`: Mass-to-charge ratio values (Daltons)
  - `i`: Intensity values

Optional columns:  
  - `file`: File ID as an integer
  - `mslevel`: MS level
  - `isowin`: Isolation window

### Extract raw XICs

Extract Extracted Ion Chromatograms (XICs) for each segment and nearby signals by extending the retention time range (by `rttol` seconds). This step ensures that all relevant signals are captured for downstream analysis, such as peak integration or quantification.

```r
xics_res <- extract_xics(
  roi_ls = ion_res,
  debugL = FALSE,
  rttol = 10
)
```
