// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// adjust_for_queue_l
List adjust_for_queue_l(List patients, IntegerVector q_vec);
RcppExport SEXP _abmpack_adjust_for_queue_l(SEXP patientsSEXP, SEXP q_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type patients(patientsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type q_vec(q_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(adjust_for_queue_l(patients, q_vec));
    return rcpp_result_gen;
END_RCPP
}
// assign_rooms_cpp_seed
Rcpp::DataFrame assign_rooms_cpp_seed(DataFrame pat_risks, SEXP icu, SEXP non, unsigned int seed);
RcppExport SEXP _abmpack_assign_rooms_cpp_seed(SEXP pat_risksSEXP, SEXP icuSEXP, SEXP nonSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type pat_risks(pat_risksSEXP);
    Rcpp::traits::input_parameter< SEXP >::type icu(icuSEXP);
    Rcpp::traits::input_parameter< SEXP >::type non(nonSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(assign_rooms_cpp_seed(pat_risks, icu, non, seed));
    return rcpp_result_gen;
END_RCPP
}
// combine_vectors
IntegerVector combine_vectors(IntegerVector x, IntegerVector y);
RcppExport SEXP _abmpack_combine_vectors(SEXP xSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(combine_vectors(x, y));
    return rcpp_result_gen;
END_RCPP
}
// sort_cpp
Rcpp::IntegerVector sort_cpp(Rcpp::IntegerVector x);
RcppExport SEXP _abmpack_sort_cpp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(sort_cpp(x));
    return rcpp_result_gen;
END_RCPP
}
// get_day_mvts_cpp
Rcpp::DataFrame get_day_mvts_cpp(IntegerVector ids, IntegerVector icu_num, IntegerVector non_num, IntegerVector end_loc);
RcppExport SEXP _abmpack_get_day_mvts_cpp(SEXP idsSEXP, SEXP icu_numSEXP, SEXP non_numSEXP, SEXP end_locSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ids(idsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type icu_num(icu_numSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type non_num(non_numSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type end_loc(end_locSEXP);
    rcpp_result_gen = Rcpp::wrap(get_day_mvts_cpp(ids, icu_num, non_num, end_loc));
    return rcpp_result_gen;
END_RCPP
}
// move_rooms_cpp_seed
Rcpp::DataFrame move_rooms_cpp_seed(DataFrame pat_rm_type, SEXP icu, SEXP non, unsigned int seed);
RcppExport SEXP _abmpack_move_rooms_cpp_seed(SEXP pat_rm_typeSEXP, SEXP icuSEXP, SEXP nonSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type pat_rm_type(pat_rm_typeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type icu(icuSEXP);
    Rcpp::traits::input_parameter< SEXP >::type non(nonSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(move_rooms_cpp_seed(pat_rm_type, icu, non, seed));
    return rcpp_result_gen;
END_RCPP
}
// adjust_for_queue
IntegerVector adjust_for_queue(IntegerVector num_pat, int q_n);
RcppExport SEXP _abmpack_adjust_for_queue(SEXP num_patSEXP, SEXP q_nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type num_pat(num_patSEXP);
    Rcpp::traits::input_parameter< int >::type q_n(q_nSEXP);
    rcpp_result_gen = Rcpp::wrap(adjust_for_queue(num_pat, q_n));
    return rcpp_result_gen;
END_RCPP
}
// num_in_q_by_risk
List num_in_q_by_risk(DataFrame risks_beds);
RcppExport SEXP _abmpack_num_in_q_by_risk(SEXP risks_bedsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type risks_beds(risks_bedsSEXP);
    rcpp_result_gen = Rcpp::wrap(num_in_q_by_risk(risks_beds));
    return rcpp_result_gen;
END_RCPP
}
// sample_day_mvts_cpp_seed
CharacterVector sample_day_mvts_cpp_seed(IntegerVector los, IntegerVector cur_room_type, unsigned int seed);
RcppExport SEXP _abmpack_sample_day_mvts_cpp_seed(SEXP losSEXP, SEXP cur_room_typeSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type los(losSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cur_room_type(cur_room_typeSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_day_mvts_cpp_seed(los, cur_room_type, seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_abmpack_adjust_for_queue_l", (DL_FUNC) &_abmpack_adjust_for_queue_l, 2},
    {"_abmpack_assign_rooms_cpp_seed", (DL_FUNC) &_abmpack_assign_rooms_cpp_seed, 4},
    {"_abmpack_combine_vectors", (DL_FUNC) &_abmpack_combine_vectors, 2},
    {"_abmpack_sort_cpp", (DL_FUNC) &_abmpack_sort_cpp, 1},
    {"_abmpack_get_day_mvts_cpp", (DL_FUNC) &_abmpack_get_day_mvts_cpp, 4},
    {"_abmpack_move_rooms_cpp_seed", (DL_FUNC) &_abmpack_move_rooms_cpp_seed, 4},
    {"_abmpack_adjust_for_queue", (DL_FUNC) &_abmpack_adjust_for_queue, 2},
    {"_abmpack_num_in_q_by_risk", (DL_FUNC) &_abmpack_num_in_q_by_risk, 1},
    {"_abmpack_sample_day_mvts_cpp_seed", (DL_FUNC) &_abmpack_sample_day_mvts_cpp_seed, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_abmpack(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
