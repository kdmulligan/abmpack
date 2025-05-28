#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <numeric>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;
using namespace arma;

//' @name sample_day_mvts_cpp_seed
//' @title Determine number of patient mvts for the day.
//' @description function to randomly sample patient daily movement patterns
//' @param los Vector of patient los for the visit,
//' @param cur_room_type vector of patients' current room type (icu/non).
//' @param seed value of seed to be passed into the rcpp.
//' @return returns character vector with the number of "icu,non" movements
//' for the patients' upcoming day.

// [[Rcpp::export]]
CharacterVector sample_day_mvts_cpp_seed(IntegerVector los, IntegerVector cur_room_type, unsigned int seed) {
  int n = los.size();
  // number of mvts for the day "icu,non"
  CharacterVector mvt_opt_char = {"0,0", "0,1", "0,2", "0,3", "1,0", "1,1", "1,2", "2,0", "2,1", "3,0"};
  IntegerVector mvt_opt = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  // Create a vector to store assigned rooms
  NumericVector pat_mvts(n);
  CharacterVector pat_mvts_char(n);

  // Initialize random number generator with seed
  Environment base = Environment::namespace_env("base");
  Function set_seed = base["set.seed"];
  set_seed(seed);

  for(int i = 0; i < n; i++) {
    // determine number of movements in the day
    if(los[i] <= 30) {
      if(cur_room_type[i] == 0) {
        // short LOS, NON room start
        NumericVector probs1 = {0.7438, 0.1469, 0.0535, 0.0338, 0.0145, 0.0042, 0.0014, 0.0012, 0.0003, 0.0003};
        pat_mvts[i] = as<IntegerVector>(Rcpp::sample(mvt_opt, 1, true, probs1))[0];
      } else {
        // short LOS, ICU room start
        NumericVector probs2 = {0.8012, 0.0964, 0.0137, 0.0016, 0.0616, 0.0103, 0.0097, 0.0047, 0.0005, 0.0004};
        pat_mvts[i] = as<IntegerVector>(Rcpp::sample(mvt_opt, 1, true, probs2))[0];
      }
    } else {
      if(cur_room_type[i] == 0) {
        // long LOS, NON room start
        NumericVector probs3 = {0.9122, 0.0534, 0.0180, 0.0062, 0.0076, 0.0012, 0.0003, 0.0006, 0.0002, 0.0002};
        pat_mvts[i] = as<IntegerVector>(Rcpp::sample(mvt_opt, 1, true, probs3))[0];
      } else {
        // long LOS, ICU room start
        NumericVector probs4 = {0.9378, 0.0121, 0.0014, 0.0002, 0.0370, 0.0079, 0.0009, 0.0023, 0.0003, 0.0002};
        pat_mvts[i] = as<IntegerVector>(Rcpp::sample(mvt_opt, 1, true, probs4))[0];
      }
    }
    pat_mvts_char[i] = mvt_opt_char[pat_mvts[i] - 1];
  }

  return pat_mvts_char;
}


