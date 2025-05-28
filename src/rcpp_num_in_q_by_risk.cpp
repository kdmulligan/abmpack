#include <Rcpp.h>
#include <numeric>
#include <iterator>
using namespace Rcpp;
using namespace std;

//' @name num_in_q_by_risk
//' @title Adjust admit patient list for the number of patients going to the queue.
//' @description Check how many patients from each facility need to go to the queue.
//' @param risks_beds Data frame with columns: hospital id, number of patients
//' at each risk level 1-4, number of icu beds, number of non-icu beds.
//' @return List with length equal to the number of facilities with the new
//' number of patients at each risk level after sending the proper number of patients
//' to the queue.

// [[Rcpp::export]]
IntegerVector adjust_for_queue(IntegerVector num_pat, int q_n) {
  // int n = x.size();
  int q_fr = 0;
  // output vector
  IntegerVector num_pat_adj(4);
  num_pat_adj = num_pat;
  // which risk levels is queue from?
  if (q_n <= num_pat[0]) {
    q_fr = 0;
  } else if (q_n <= num_pat[0] + num_pat[1]) {
    q_fr = 1;
  } else if (q_n <= num_pat[0] + num_pat[1] + num_pat[2]) {
    q_fr = 2;
  } else if (q_n <= num_pat[0] + num_pat[1] + num_pat[2] + num_pat[3]) {
    q_fr = 3;
  }
  // return q_fr;
  int sum_v = std::accumulate(num_pat.begin(), num_pat.begin() + q_fr + 1 , 0.0);
  // for(int i = 0; i < 4; i++) {
  int j = 0;
  while(j < 4) {
    if(j < q_fr) {
      num_pat_adj[j] = 0;
    } else if (j == q_fr) {
      // int sum_v = std::accumulate(num_pat.begin(), num_pat.begin() + q_fr , 0);
      num_pat_adj[j] = sum_v - q_n;
    }
    j++;
  }
  return num_pat_adj;
}


// // [[Rcpp::export]]
// // Function to sample n elements from a vector without replacement
// IntegerVector sub_vec(IntegerVector vec_1, IntegerVector vec_2) {
//   IntegerVector result(vec_1.size());
//   result = vec_1 - vec_2;
//   return result;
// }

// [[Rcpp::export]]
List num_in_q_by_risk(DataFrame risks_beds) {
  NumericVector id = risks_beds[0];
  NumericVector risk_1 = risks_beds[1];
  NumericVector risk_2 = risks_beds[2];
  NumericVector risk_3 = risks_beds[3];
  NumericVector risk_4 = risks_beds[4];
  NumericVector n_icu = risks_beds[5];
  NumericVector n_non = risks_beds[6];
  int n = id.size();
  List out(n);

  for(int i = 0; i < n; i++) {

    IntegerVector pat_i = IntegerVector::create(risk_1[i], risk_2[i], risk_3[i], risk_4[i]) ;

    int tot_beds = n_icu[i] + n_non[i];
    int tot_ad_pat = risk_1[i] + risk_2[i] + risk_3[i] + risk_4[i];
    int q_num = 0;
    if(tot_beds < tot_ad_pat) {
      q_num = tot_ad_pat - tot_beds;
    }
    // out[i] = adjust_for_queue(pat_i, q_num);
    IntegerVector temp = adjust_for_queue(pat_i, q_num);
    out[i] = IntegerVector::create(risk_1[i] - temp[0], risk_2[i] - temp[1], risk_3[i] - temp[2], risk_4[i] - temp[3]);
  }
  out.attr("names") = id;
  return out;
}

