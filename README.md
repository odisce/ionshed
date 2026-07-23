<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/odisce/ionshed/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/odisce/ionshed/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**ionshed** is an R package for sparse matrix segmentation inspired by watershed and implemented in C++. It works on .mz(X)ML data or directly on 3D tables.
It will parse 3D coordinates by decreasing intensities and segment (cluster) the data points using their proximity in two dimensions.

The algorithm was developped to segment LC-MS data (MS1 and DIA), hence the segmentation is dictated by the following parameters:
  - `rttol`: corresponds to the chromatographic resolution in seconds. It's the minimal distance between two points to be able to separate two chromatogrphic peaks.
  - `ppm`: corresponds to the mass resolution in Da. The distance between two signals to differentiate specific ions.

The algorithm can be used on any 3D datasets by giving their coordinates in a table (see #).

## Requirements

- To use in R: [R>=4.0.0](https://cloud.r-project.org/)
- To compile from source: [Rtools](https://cloud.r-project.org/)

## Installation

The latest **ionshed** version can be installed from source using your prefered method (need [Rtools](https://cloud.r-project.org/)):

| method | command |
|--|--|
| [remotes](https://cran.r-project.org/web/packages/remotes/index.html) | `remotes::install_github("odisce/ionshed")` |
| [devtools](https://cran.r-project.org/web/packages/devtools/) | `devtools::install_github("odisce/ionshed")` |
| [renv](https://cran.r-project.org/web/packages/renv/) | `renv::install("github::odisce/ionshed")` |
| [pak](https://cran.r-project.org/web/packages/pak/) | `pak::pkg_install("odisce/ionshed")` |

## Usage

### From mz(X)ML data

Load `.mzML` files:

| format | command |
| :--- | :--- |
| using `MSnbase` | `msn_in <- MSnbase::readMSData("path/to/my/.mzml", mode = "onDisk", msLevel. = 1L)` |
| using `MsExperiment` | `msn_in <- MsExperiment::readMsExperiment(spectraFiles = "path/to/my/.mzml")` |


Run the `ionshed` algorithm on the first file `file_ind = 1`:

```r
ion_res <- ionshed::ionshed(
  data = msn_in,
  file_ind = 1,
  rttol = 5,
  ppm = 3
)
```

The results are stored in a list with two levels:
  - original segmented coordinates `ion_res["roi_data"]`:
    ```r
    str(ion_res$roi_data)
    ```
  - a summary of each segments `ion_res["peak_info"]`:
    ```r
    str(ion_res$peak_info)
    ```


### From 3D coordinates

To perform the segmentation on coordinates, the function needs a table with the following mandatory columns:
  - `rt`: rt values in seconds
  - `mz`: mz values in Daltons
  - `i`: intensity values

And some optional:  
  - `file`: file id as an integer
  - `mslevel`: MS level
  - `isowin`: Isolation window

### Extract raw XICs

A derived function allow to extract XICs for each segments and close signals by extending the retention time range by `rttol` seconds:

```r
xics_res <- extract_xics(
  roi_ls = ion_res,
  debugL = FALSE,
  rttol = 10
)
```
