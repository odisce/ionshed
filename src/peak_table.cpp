#include <Rcpp.h>
#include "utils.h"
#include <chrono>
using namespace Rcpp;

// Define a struct to store summary statistics
// struct SummaryStats {
//     int roi, sc_length, rtid_start, rtid_end, mzid_start, mzid_end;
//     double rtmean, rtmin, rtmax, mzmin, mzmax, rtmed, rtmeanw, mzmed, mzmeanw;
//     double imax, imin, idiff, imedian, imean, into;
//     double mzdiff, rtdiff;
// };

double weightedMean(NumericVector x, NumericVector y) {
    double sum = 0, numWeight = 0;
    for (int i = 0; i!=x.size(); ++i) {
        numWeight = numWeight + x[i] * y[i];
        sum += y[i];
    }
    double out = numWeight / sum;
    return out;
}

// Define a struct to store summary statistics
// struct SummaryStat {
//     std::map<std::string, double> value;
    
//     SummaryStat() {
//         value["roi"] = 0;
//         value["sc_length"] = 0;
//         value["rtid_start"] = 0;
//         value["rtid_end"] = 0;
//         value["mzid_start"] = 0;
//         value["mzid_end"] = 0;
//         value["rtmean"] = 0;
//         value["rtmin"] = 0;
//         value["rtmax"] = 0;
//         value["mzmin"] = 0;
//         value["mzmax"] = 0;
//         value["rtmed"] = 0;
//         value["rtmeanw"] = 0;
//         value["mzmed"] = 0;
//         value["mzmeanw"] = 0;
//         value["imax"] = 0;
//         value["imin"] = 0;
//         value["idiff"] = 0;
//         value["imedian"] = 0;
//         value["imean"] = 0;
//         value["into"] = 0;
//         value["mzdiff"] = 0;
//         value["rtdiff"] = 0;
//     }

//     SummaryStat(
//         int roi,
//         NumericVector rt_nv,
//         NumericVector i_nv,
//         NumericVector mz_nv,
//         IntegerVector mzsc_iv,
//         IntegerVector rtsc_iv
//     ) {
//         value["roi"] = roi;
//         value["sc_length"] = rt_nv.size();
//         value["rtid_start"] = min(rtsc_iv);
//         value["rtid_end"] = max(rtsc_iv);
//         value["mzid_start"] = min(mzsc_iv);
//         value["mzid_end"] = max(mzsc_iv);
//         value["rtmean"] = mean(rt_nv);
//         value["rtmin"] = min(rt_nv);
//         value["rtmax"] = max(rt_nv);
//         value["mzmin"] = min(mz_nv);
//         value["mzmax"] = max(mz_nv);
//         value["rtmed"] = median(rt_nv);
//         value["rtmeanw"] = weightedMean(rt_nv, i_nv);
//         value["mzmed"] = median(mz_nv);
//         value["mzmeanw"] = weightedMean(mz_nv, i_nv);
//         value["imax"] = max(i_nv);
//         value["imin"] = min(i_nv);
//         value["idiff"] = value["imax"] - value["imin"];
//         value["imedian"] = median(i_nv);
//         value["imean"] = mean(i_nv);
//         value["into"] = sum(i_nv);
//         value["mzdiff"] = value["mzmax"] - value["mzmin"];
//         value["rtdiff"] = value["rtmax"] - value["rtmin"];
//     }
// };

// // Function to calculate summary statistics for a NumericVector
// SummaryStats calculate_stats(
//     NumericVector rt_nv,
//     NumericVector i_nv,
//     NumericVector mz_nv,
//     IntegerVector mzsc_iv,
//     IntegerVector rtsc_iv
// ) {
//     SummaryStats stats;
//     stats.sc_length = rt_nv.size();
//     stats.rtmean = mean(rt_nv);
//     stats.rtmin = min(rt_nv);
//     stats.rtmax = max(rt_nv);
//     stats.rtmed = median(rt_nv);
//     stats.rtmeanw = weightedMean(rt_nv, i_nv);
//     stats.rtdiff = stats.rtmax - stats.rtmin;
//     stats.rtid_start = min(rtsc_iv);
//     stats.rtid_end = max(rtsc_iv);
//     stats.mzmin = min(mz_nv);
//     stats.mzmax = max(mz_nv);
//     stats.mzmed = median(mz_nv);
//     stats.mzmeanw = weightedMean(mz_nv, i_nv);
//     stats.mzdiff = stats.mzmax - stats.mzmin;
//     stats.mzid_start = min(mzsc_iv);
//     stats.mzid_end = max(mzsc_iv);
//     stats.imax = max(i_nv);
//     stats.imin = min(i_nv);
//     stats.imedian = median(i_nv);
//     stats.imean = mean(i_nv);
//     stats.idiff = stats.imax - stats.imin;
//     stats.into = sum(i_nv);
//     return stats;
// }

std::map<std::string, std::vector<double>> SummaryObj(int size) {
    std::map<std::string, std::vector<double>> results;
    std::vector<std::string> metrics = {
        "roi",
        "sc_length",
        "rtid_start",
        "rtid_end",
        "mzid_start",
        "mzid_end",
        "rtmean",
        "rtmin",
        "rtmax",
        "mzmin",
        "mzmax",
        "rtmed",
        "rtmeanw",
        "mzmed",
        "mzmeanw",
        "imax",
        "imin",
        "idiff",
        "imedian",
        "imean",
        "into",
        "mzdiff",
        "rtdiff",
    };
    for (auto meti : metrics) {
        results[meti] = std::vector<double> (size, 0);
    }
    return results;
}

void UpdateObj(
    std::map<std::string, std::vector<double>> &value,
    int vec_index,
    double roi,
    NumericVector rt_nv,
    NumericVector i_nv,
    NumericVector mz_nv,
    IntegerVector mzsc_iv,
    IntegerVector rtsc_iv
) {
    value["roi"][vec_index] = roi;
    value["sc_length"][vec_index] = rt_nv.size();
    value["rtid_start"][vec_index] = min(rtsc_iv);
    value["rtid_end"][vec_index] = max(rtsc_iv);
    value["mzid_start"][vec_index] = min(mzsc_iv);
    value["mzid_end"][vec_index] = max(mzsc_iv);
    value["rtmean"][vec_index] = mean(rt_nv);
    value["rtmin"][vec_index] = min(rt_nv);
    value["rtmax"][vec_index] = max(rt_nv);
    value["mzmin"][vec_index] = min(mz_nv);
    value["mzmax"][vec_index] = max(mz_nv);
    value["rtmed"][vec_index] = median(rt_nv);
    value["rtmeanw"][vec_index] = weightedMean(rt_nv, i_nv);
    value["mzmed"][vec_index] = median(mz_nv);
    value["mzmeanw"][vec_index] = weightedMean(mz_nv, i_nv);
    value["imax"][vec_index] = max(i_nv);
    value["imin"][vec_index] = min(i_nv);
    value["idiff"][vec_index] = (max(i_nv) - min(i_nv));
    value["imedian"][vec_index] = median(i_nv);
    value["imean"][vec_index] = mean(i_nv);
    value["into"][vec_index] = sum(i_nv);
    value["mzdiff"][vec_index] = (max(mz_nv) - min(mz_nv));
    value["rtdiff"][vec_index] = (max(rt_nv) - min(rt_nv));
}


// [[Rcpp::export]]
Rcpp::DataFrame summarize_one_class(
    double roi,
    NumericVector rt_nv,
    NumericVector i_nv,
    NumericVector mz_nv,
    IntegerVector mzsc_iv,
    IntegerVector rtsc_iv
) {
    std::map<std::string, std::vector<double>> results = SummaryObj(1);
    UpdateObj(
        results,
        0,
        roi,
        rt_nv,
        i_nv,
        mz_nv,
        mzsc_iv,
        rtsc_iv
    );
    // Convert to Data.Frame
    Rcpp::DataFrame output;
    for (auto const& resi : results) {
        output.push_back( resi.second, resi.first );
    }
    return output;
}

// [[Rcpp::export]]
Rcpp::DataFrame summarize_by_class(
    DataFrame df,
    bool debugL = false
) {
    using namespace std::chrono;
    // Init chronos
    auto timeA = high_resolution_clock::now();
    auto timeB = high_resolution_clock::now();
    auto timeBB = high_resolution_clock::now();
    auto timeBBB = high_resolution_clock::now();
    auto timeC = high_resolution_clock::now();
    auto timeD = high_resolution_clock::now();
    auto timeE = high_resolution_clock::now();
    // Extract columns from the DataFrame
    NumericVector col_mz = df["mz"];
    NumericVector col_rt = df["rt"];
    NumericVector col_i = df["i"];
    IntegerVector col_mzrt = df["mz_sc"];
    IntegerVector col_rtsc = df["rt_sc"];
    IntegerVector col_roi = df["roi"];
    
    if (debugL) {
        timeB = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(timeB - timeA);
        Rcpp::Rcout << "Initialize input: " << duration.count() << "ms\n";
    }
    // Get unique classes
    std::vector<size_t> class_indexes = sorted_index_int(col_roi);
    if (debugL) {
        timeBB = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(timeBB - timeB);
        Rcpp::Rcout << "Initialize class vector: " << duration.count() << "ms\n";
    }
    // First count unique rois
    int roi_nb = 0;
    int prev_i = 0;
    for (int i=0;i!=col_roi.size();i++) {
        const int roi_i = col_roi[class_indexes[i]];
        if (roi_i == prev_i) {
            continue;
        } else {
            prev_i = roi_i;
            roi_nb++;
        }
    }
    if (debugL) {
        timeBBB = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(timeBBB - timeBB);
        Rcpp::Rcout << "Count unique rois (" << roi_nb << "): " << duration.count() << "ms\n";
    }
    // Get unique rois
    IntegerVector unique_classes (roi_nb);
    prev_i = 0;
    int roi_cnt = 0;
    for (int i=0;i!=col_roi.size();i++) {
        const int roi_i = col_roi[class_indexes[i]];
        if (roi_i == prev_i) {
            continue;
        } else {
            unique_classes[roi_cnt] = roi_i;
            roi_cnt++;
            prev_i = roi_i;
        }
    }
    if (debugL) {
        timeC = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(timeC - timeB);
        Rcpp::Rcout << "Getting unique rois: " << duration.count() << "s\n";
    }
    // Iterate over each rois
    std::size_t startv_i = 0;
    std::map<std::string, std::vector<double>> results = SummaryObj(unique_classes.size());
    if (debugL) {
        timeD = high_resolution_clock::now();
    }
    for (int cl_it = 0; cl_it != unique_classes.size();  ++cl_it ) {
        // Get indexes for this groups
        int classi = unique_classes[cl_it];
        IntegerVector class_vind;
        std::size_t last_index = 0;
        for (std::size_t indi=startv_i; indi!=class_indexes.size(); ++indi) {
            int cur_class = col_roi[class_indexes[indi]];
            if (cur_class == classi) {
                class_vind.push_back(static_cast<int>(class_indexes[indi]));
                last_index = indi+1;
            } else {
                break;
            }
        }
        // Change startv_i for next interation
        startv_i = last_index;
        // Get data
        UpdateObj(
            results,
            cl_it,
            classi,
            col_rt[class_vind],
            col_i[class_vind],
            col_mz[class_vind],
            col_mzrt[class_vind],
            col_rtsc[class_vind]
        );
    }
    if (debugL) {
        timeE = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(timeE - timeD);
        Rcpp::Rcout << "Processed all rois: " << duration.count() << "s\n";
    }
    // Convert to Data.Frame
    Rcpp::DataFrame output;
    for (auto const& resi : results) {
        output.push_back( resi.second, resi.first );
    }
    return output;
}