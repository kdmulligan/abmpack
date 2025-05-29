#include <Rcpp.h>
#include <numeric>
#include <iterator>
#include <list>

using namespace Rcpp;
using namespace std;

//' @name adjust_for_queue_l
//' @title Adjust admit queue list
//' @description Check how many patients from each facility need to go to the queue.
//' @param patients List of patients at each facility
//' @param q_vec number in queue total
//' @return adjusted queue list.
 // [[Rcpp::export]]
IntegerVector adjust_for_queue_i(IntegerVector num_pat, int q_n) {
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
// [[Rcpp::export]]
List adjust_for_queue_l(List patients, IntegerVector q_vec) {
  int n = q_vec.size();

  // output vector
  List num_pat_adj_l = List::create();

  for(int i = 0; i < n; i++) {
    IntegerVector pat_i = patients[i];
    int q_i = q_vec[i];
    IntegerVector temp = adjust_for_queue_i(pat_i, q_i);
    num_pat_adj_l[i] = temp;
  }
  return num_pat_adj_l;

}
