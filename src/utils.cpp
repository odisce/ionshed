#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double min_diff(std::vector<double> mz_vec) {
  using namespace Rcpp;
  using namespace std;
  // Get unique values
  // Search min diff
  std::sort(mz_vec.begin(), mz_vec.end());
  mz_vec.erase(unique( mz_vec.begin(), mz_vec.end() ), mz_vec.end());
  double mz_min = 1000;
  double mz_diff = 1000;
  for(size_t it = 0; it != mz_vec.size()-1; it++) {
    mz_diff = mz_vec[it+1] - mz_vec[it];
    if (mz_diff < mz_min) {
      mz_min = mz_diff;
    }
  }
  return mz_min;
}

template <typename T>
void sort_indices(std::vector<T> &data, std::vector<size_t> &indices){
    std::sort(indices.begin(), indices.end(), [&data](size_t a, size_t b){ return data[a] < data[b]; });
}

std::vector<size_t> sorted_index(Rcpp::NumericVector& indt) {
  using namespace Rcpp;
  using namespace std;
  std::vector<double> indt_v = as<std::vector<double> >(indt);
  std::vector<size_t> indt_v_index(indt_v.size());
  for (size_t it = 0 ; it != indt_v_index.size() ; it++) {
    indt_v_index[it] = it;
  }
  sort_indices(indt_v, indt_v_index);
  return indt_v_index;
}

std::vector<size_t> sorted_index_int(Rcpp::IntegerVector& indt) {
  using namespace Rcpp;
  using namespace std;
  std::vector<int> indt_v = as<std::vector<int> >(indt);
  std::vector<size_t> indt_v_index(indt_v.size());
  for (size_t it = 0 ; it != indt_v_index.size() ; it++) {
    indt_v_index[it] = it;
  }
  sort_indices(indt_v, indt_v_index);
  return indt_v_index;
}

// [[Rcpp::export]]
double min_diff_by_group(Rcpp::NumericVector in_vec, Rcpp::NumericVector group_vec) {
  using namespace Rcpp;
  using namespace std;
  std::vector<double> mz_vec = as<std::vector<double>>(in_vec);
  // sort groups indexes
  std::vector<size_t> grp_index = sorted_index(group_vec);
  double cur_grp = group_vec[grp_index[0]];
  std::vector<double> grp_val;
  std::vector<double> output_diffs;
  for(size_t grp_ind = 0; grp_ind != grp_index.size(); grp_ind++) {
    size_t grpi = grp_index[grp_ind];
    if (group_vec[grpi] == cur_grp) {
      // same group, add to grp_val
      grp_val.push_back(in_vec[grpi]);
    }
    if (grp_ind == grp_index.size() - 1 || group_vec[grpi] != cur_grp) {
      // different group, calculate previous min diff and add to vector
      double diff_res = min_diff(grp_val);
      output_diffs.push_back(diff_res);
      cur_grp = group_vec[grpi];
      grp_val.clear();
    }
  }
  std::vector<double>::iterator res = std::min_element(output_diffs.begin(), output_diffs.end());
  return *res;
}
