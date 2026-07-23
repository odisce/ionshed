#include <RcppArmadillo.h>
#include "utils.h"
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::depends(RcppArmadillo)]]

//' Function to attribute ROI by taking the highest intensity points
//'
//' This function returns a matrix with m/Z groups
//' based on a tolerance. This version is the final one.
//'
//' @param locations A matrix (2xN) with row (1xN) and column (1xN)
//' @param values A vector of N values
//' @param rttol Tolerance in rt dimension (grid points)
//' @param ppm Tolerance in ppm (variable grid)
//' @param debugL Boolean to print debug messages
//' @param eachiter Boolean to return iteration steps details
//' @export
//'
// [[Rcpp::export]]
arma::mat SpMat_DescendingROI(
  const arma::umat locations,
  const arma::vec values,
  const arma::uword rttol = 10,
  const double ppm = 10,
  const bool debugL = false,
  const bool eachiter = false
) {
  using namespace std;
  using namespace arma;

  // Generate output table
  int outcol = 5;
  if (eachiter) {
    outcol = 14;
  }
  arma::mat out(locations.n_cols, outcol);
  for(arma::uword i=0;i!=locations.n_cols;++i) {
    out(i,0) = locations(0,i); //rt
    out(i,1) = locations(1,i); //mz
    out(i,2) = values(i);
  }

  // Sort reference table by intensity
  const arma::uvec sort_ref_int = sort_index(out.col(2), "descend");
  // For each, get submat and assign ROI
  const arma::uvec values_roi(values.n_elem, arma::fill::ones);
  arma::sp_umat XROI = arma::sp_umat(locations, values_roi);
  if (debugL) {
    Rcpp::Rcout << "cols: " << XROI.n_cols << " rows: " << XROI.n_rows << endl;
  }
  const int amaxmz = XROI.n_cols; //mz in col
  const int amaxrt = XROI.n_rows; //rt in row
  const double ppm_val = ppm/1E6;

  int roi_cnt = 2;

  for(arma::uword i=0 ; i!=sort_ref_int.n_rows ; ++i) {
    int pointtype = 0;
    arma::uword roi_sel = 0;
    int row_i = sort_ref_int(i);
    const arma::uword int_rowind = out(row_i,0); //rt
    const arma::uword int_colind = out(row_i,1); //mz
    arma::vec artrange(2);
    arma::vec amzrange(2);
    arma::uword i_sub_row_ref = rttol;
    double mztol_dbl = ppm_val*int_colind;
    arma::uword mztol = ceil(mztol_dbl);

    if (int_colind <= mztol) {
      amzrange[0] = 0;
    } else {
      amzrange[0] = int_colind-mztol;
    }
    if (int_colind >= (amaxmz-1-mztol)) {
      amzrange[1] = amaxmz-1;
    } else {
      amzrange[1] = int_colind+mztol;
    }
    if (int_rowind <= rttol) {
      artrange[0] = 0;
      i_sub_row_ref = int_rowind;
    } else {
      artrange[0] = int_rowind-rttol;
    }
    if (int_rowind >= (amaxrt-1-rttol)) {
      artrange[1] = amaxrt-1;
    } else {
      artrange[1] = int_rowind+rttol;
    }

    if (eachiter) {
      out(row_i,5) = amzrange[0];
      out(row_i,6) = amzrange[1];
      out(row_i,7) = artrange[0];
      out(row_i,8) = artrange[1];
      out(row_i,9) = mztol;
      out(row_i,10) = ppm_val;
      out(row_i,11) = mztol_dbl;
      out(row_i,12) = ppm;
      out(row_i,13) = i;
    }
    arma::uword block_nzero = XROI.submat(artrange[0], amzrange[0], artrange[1], amzrange[1]).n_nonzero;
    if (block_nzero == 1) {
      // only one point, go to next
      roi_sel = 1;
      pointtype = 1;
    } else if (block_nzero > 1) {
      // Iterate over block element
      // We should not compute unique value,
      // find a way to get the uniques values here (or at least the closest one only)
      arma::sp_umat roi_block = XROI.submat(artrange[0], amzrange[0], artrange[1], amzrange[1]);
      arma::uvec roi_block_val(roi_block.n_nonzero);
      arma::uvec roi_block_row(roi_block.n_nonzero);
      arma::uword only_one=1;
      arma::uword only_one_cnt=0;
      arma::uword only_one_other=1;
      arma::uword only_one_other_val=0;
      int z=0;
      for(auto it=roi_block.begin() ; it!=roi_block.end() ; ++it, ++z) {
        if (*it == 1) {
          roi_block_val(z) = 0;
        } else if (*it > 1) {
          only_one=0;
          if (only_one_cnt==0) {
            only_one_cnt=1;
            only_one_other_val=*it;
          } else {
            if (*it != only_one_other_val) {
              only_one_other=0;
            }
          }
          roi_block_row(z) = it.row();
          roi_block_val(z) = *it;
        }
      }

      // if only_one == 1, then only roi 1 in block = initiate new roi
      // if only_one_other == 0 means multiple roi != 1
      // if only_one_other == 1 means only one roi != 1 and value in only_one_other_val

      if (only_one == 1) {
        roi_sel = roi_cnt;
        roi_cnt = roi_cnt + 1;
        pointtype = 2;
      } else if (only_one_other == 1) {
        roi_sel = only_one_other_val;
        pointtype = 3;
      } else if (only_one_other == 0) {
        // Multiple roi, search closest roi != 1 in rt (nrow)
        arma::uword min_val_index = rttol+1;
        for (arma::uword iz=0;iz!=roi_block_row.n_elem;++iz) {
          if (roi_block_val(iz) > 0) {
            arma::sword diffval = roi_block_row(iz) - i_sub_row_ref;
            arma::uword absdiffval = abs(diffval);
            if (absdiffval < min_val_index) {
              min_val_index = abs(diffval);
              roi_sel = roi_block_val(iz);
            }
            // if (debugL) {
            //   Rcpp::Rcout << "diffval: " << diffval << " ";
            //   Rcpp::Rcout << "absdiffval: " << absdiffval << " ";
            //   Rcpp::Rcout << "min_val_index: " << min_val_index << " ";
            //   Rcpp::Rcout << "roi_sel: " << roi_sel << " " << "\n";
            // }
          }
        }
        pointtype=4;
      } else {
        roi_sel = 1;
        pointtype = 5;
      }
    }
    // Update SpMat
    XROI(int_rowind, int_colind) = roi_sel;
    // Export Results
    out(row_i,3) = roi_sel;
    out(row_i,4) = pointtype;
  }
  return out;
}

//' Function to extract EIcs using Sparse Matrix
//'
//' This function returns a matrice with EICs corresponding to asked range(s).
//'
//' @param locations A matrix (2xN) with row (1xN) and column (1xN)
//' @param values A matrix (ValuesxN) with values in the same order as locations
//' @param eics_dt coordinate as a matrix with id, mzmin, mzmax, rtmin, rtmax
//' @param rttol (optional) to extend the rt range
//' @param debugL Boolean to print debug message
//' @export
//'
// [[Rcpp::export]]
Rcpp::DataFrame SpMat_getEICs(
  const arma::umat locations,
  const arma::dmat values,
  const arma::umat eics_dt,
  const arma::uword rttol = 10,
  const bool debugL = false
) {
  using namespace std;
  using namespace arma;
  using namespace std::chrono;
  auto timeA = high_resolution_clock::now();
  auto timeB = high_resolution_clock::now();
  auto timeC = high_resolution_clock::now();
  auto timeD = high_resolution_clock::now();
  // initalize output array
  if (debugL) {
    timeA = high_resolution_clock::now();
    Rcpp::Rcout << "-- Starting XICs extraction --";
  }
  arma::field<std::vector<double>> V(6);
  // Create Sparse Matrix
  arma::field<arma::sp_mat> F(4);
  F(0) = arma::sp_mat(locations, values.col(0));
  F(1) = arma::sp_mat(locations, values.col(1));
  F(2) = arma::sp_mat(locations, values.col(2));
  F(3) = arma::sp_mat(locations, values.col(3));
  arma::vec rtrange(2);
  arma::vec mzrange(2);
  // Iterate over eics_dt
  // const arma::uword maxcols = F(0).n_cols;
  const arma::uword maxrows = F(0).n_rows;
  if (debugL) {
    timeB = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(timeB - timeA);
    Rcpp::Rcout << "Initialize input: " << duration.count() << "ms\n";
  }
  for (arma::uword i=0 ; i!=eics_dt.n_rows ; ++i) {
    arma::uword row_start_new = eics_dt(i, 1), row_end_new = eics_dt(i, 2);
    arma::uword col_start_new = eics_dt(i, 3), col_end_new = eics_dt(i, 4);
    const arma::uword roi_i = eics_dt(i, 0);
    if (rttol > 0) {
      if (row_start_new < rttol) {
        row_start_new = 0;
      } else {
        row_start_new = (row_start_new - rttol);
      }
      if ( (row_end_new + rttol) > (maxrows - 1) ) {
        row_end_new = maxrows - 1;
      } else {
        row_end_new = (row_end_new + rttol);
      }
    }
    for (arma::uword y=0; y != F.size(); y++) {
      arma::sp_dmat XROI_submatrix = F(y).submat(row_start_new, col_start_new, row_end_new, col_end_new);
      // Try iterators
      arma::sp_mat::const_iterator it     = XROI_submatrix.begin();
      arma::sp_mat::const_iterator it_end = XROI_submatrix.end();
      for(; it != it_end; ++it) {
        V(y).push_back(*it);
        if (y == 3) {
          V(4).push_back(roi_i);
          if (*it == roi_i) {
            V(5).push_back(1);
          } else {
            V(5).push_back(0);
          }
        }
      }
    }
  }
  if (debugL) {
    timeC = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(timeC-timeB);
    Rcpp::Rcout << "Extracted XICs: " << duration.count() << "ms\n";
  }
  Rcpp::DataFrame output = Rcpp::DataFrame::create(
    Rcpp::Named("rt") = V(0),
    Rcpp::Named("mz") = V(1),
    Rcpp::Named("i") = V(2),
    Rcpp::Named("roi") = V(3),
    Rcpp::Named("xic_roi") = V(4),
    Rcpp::Named("assignment") = V(5)
  );
  if (debugL) {
    timeD = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(timeD-timeC);
    Rcpp::Rcout << "Formating results: " << duration.count() << "ms\n";
  }
  return output;
}