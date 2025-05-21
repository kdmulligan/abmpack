#include <RcppArmadillo.h>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::export]]
IntegerVector combine_vectors(IntegerVector x, IntegerVector y) {
  // Create a new vector to store the combined result
  IntegerVector result(x.size() + y.size());
  // Copy elements from x to result
  std::copy(x.begin(), x.end(), result.begin());
  // Copy elements from y to result
  std::copy(y.begin(), y.end(), result.begin() + x.size());
  return result;
}

// [[Rcpp::export]]
Rcpp::IntegerVector sort_cpp(Rcpp::IntegerVector x) {
  std::sort(x.begin(), x.end());
  return x;
}

// [[Rcpp::export]]
Rcpp::DataFrame get_day_mvts_cpp(IntegerVector ids, IntegerVector icu_num, IntegerVector non_num, IntegerVector end_loc) {
  int n_pat = ids.size();
  int n_row = sum(icu_num) + sum(non_num);

  IntegerVector room_type(n_row);
  IntegerVector times(n_row);
  arma::mat mat_comb(0, 3);

  // iterate over all facilities admitting patients
  for(int i = 0; i < n_pat; i++) {
    int n_pat_moves = non_num[i] + icu_num[i];

    // for(int j = 0; j < n_pat_moves; j++) {
    // get middle patient moves
    IntegerVector middle_loc(n_pat_moves - 1);
    if(end_loc[i] == 1) {
      middle_loc = combine_vectors(rep(1, icu_num[i] - 1), rep(0, non_num[i]));
    } else {
      middle_loc = combine_vectors(rep(1, icu_num[i]), rep(0, non_num[i] - 1));
    }

    // sequence of movements for the day
    IntegerVector seq_mvt(n_pat_moves);
    if(icu_num[i] == 2 || non_num[i] == 2) {
      IntegerVector mid_order = Rcpp::sample(middle_loc, n_pat_moves - 1, false);
      seq_mvt = combine_vectors(mid_order, IntegerVector::create(end_loc[i]));
    } else {
      seq_mvt = combine_vectors(middle_loc, IntegerVector::create(end_loc[i]));
    }

    // times of movements
    IntegerVector time_mvts = Rcpp::sample(IntegerVector::create(1, 2, 3, 4), n_pat_moves, false);
    time_mvts = sort_cpp(time_mvts);

    arma::mat mat_i(n_pat_moves, 3);
    mat_i.col(0) = arma::vec(n_pat_moves, fill::value(ids[i]));
    mat_i.col(1) = as<arma::vec>(seq_mvt);
    mat_i.col(2) = as<arma::vec>(time_mvts);

    mat_comb = join_vert(mat_comb, mat_i);
  }
  // }

  // Create a DataFrame from the combined matrix
  return DataFrame::create(
    Named("patid") = mat_comb.col(0),
    Named("room_type") = mat_comb.col(1),
    Named("time_block") = mat_comb.col(2)
  );
}
