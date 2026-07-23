#include <Rcpp.h>
using namespace Rcpp;
//' Function to compute the ZigZag index
//'
//' This function returns the ZigZag index as a numeric
//' value of a single trace.
//'
//' @param intPts Numeric vector of intensity pre-ordered by time
//' @export
//'
// [[Rcpp::export]]
double zigzag_index(std::vector<double> intPts) {
  const size_t eic_l = intPts.size();
  double zig_zag_index = NAN;
  std::vector<double>::iterator max_it;
  max_it = std::max_element(intPts.begin(), intPts.end());
  if (eic_l > 4) {
    double EPI=*max_it-(intPts[1]+intPts[2]+intPts[eic_l]+intPts[eic_l-1])/4.0;
    double zig_zag_sum = 0.0;
    for (size_t i=2;i<eic_l;i++) {
      double local_zig_zag=std::pow((2*intPts[i]-intPts[i-1]-intPts[i+1]), 2.0);
      zig_zag_sum=zig_zag_sum+local_zig_zag;
    }
    zig_zag_index = zig_zag_sum/(std::pow(EPI, 2.0)*eic_l);
  }
  return zig_zag_index;
}