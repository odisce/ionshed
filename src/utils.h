#ifndef __UTILS_MINDIFF__
#define __UTILS_MINDIFF__

double min_diff(std::vector<double> mz_vec);

#endif // __UTILS_MINDIFF__


#ifndef __UTILS_SORTEDINDEX__
#define __UTILS_SORTEDINDEX__

std::vector<size_t> sorted_index(Rcpp::NumericVector& indt);

#endif // __UTILS_SORTEDINDEX__

#ifndef __UTILS_SORTEDINDEX_INT__
#define __UTILS_SORTEDINDEX_INT__

std::vector<size_t> sorted_index_int(Rcpp::IntegerVector& indt);

#endif // __UTILS_SORTEDINDEX_INT__

#ifndef __UTILS_MINDIFFBYGRP__
#define __UTILS_MINDIFFBYGRP__

double min_diff_by_group(Rcpp::NumericVector in_vec, Rcpp::NumericVector group_vec);

#endif // __UTILS_MINDIFFBYGRP__