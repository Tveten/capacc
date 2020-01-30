// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// hash_map_test
void hash_map_test();
RcppExport SEXP _mvcapaCor_hash_map_test() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    hash_map_test();
    return R_NilValue;
END_RCPP
}
// optimise_savings
Rcpp::List optimise_savings(const arma::mat& Q, const arma::vec& b, const Rcpp::List& penalty_components, const Rcpp::List& extended_nbs_list);
RcppExport SEXP _mvcapaCor_optimise_savings(SEXP QSEXP, SEXP bSEXP, SEXP penalty_componentsSEXP, SEXP extended_nbs_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type penalty_components(penalty_componentsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type extended_nbs_list(extended_nbs_listSEXP);
    rcpp_result_gen = Rcpp::wrap(optimise_savings(Q, b, penalty_components, extended_nbs_list));
    return rcpp_result_gen;
END_RCPP
}
// optimise_savings_list
Rcpp::List optimise_savings_list(const arma::mat& Q, const arma::vec& b, const Rcpp::List& penalty_components, const Rcpp::List& extended_nbs_list);
RcppExport SEXP _mvcapaCor_optimise_savings_list(SEXP QSEXP, SEXP bSEXP, SEXP penalty_componentsSEXP, SEXP extended_nbs_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type penalty_components(penalty_componentsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type extended_nbs_list(extended_nbs_listSEXP);
    rcpp_result_gen = Rcpp::wrap(optimise_savings_list(Q, b, penalty_components, extended_nbs_list));
    return rcpp_result_gen;
END_RCPP
}
// list_test
void list_test(Rcpp::List L);
RcppExport SEXP _mvcapaCor_list_test(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type L(LSEXP);
    list_test(L);
    return R_NilValue;
END_RCPP
}
// seq_int
std::vector<int> seq_int(const int& start, const int& end);
RcppExport SEXP _mvcapaCor_seq_int(SEXP startSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type start(startSEXP);
    Rcpp::traits::input_parameter< const int& >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_int(start, end));
    return rcpp_result_gen;
END_RCPP
}
// test_map
std::unordered_map<std::string, double> test_map();
RcppExport SEXP _mvcapaCor_test_map() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(test_map());
    return rcpp_result_gen;
END_RCPP
}
// test_ubt
void test_ubt(const int& n_levels);
RcppExport SEXP _mvcapaCor_test_ubt(SEXP n_levelsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n_levels(n_levelsSEXP);
    test_ubt(n_levels);
    return R_NilValue;
END_RCPP
}
// which_max_test
void which_max_test();
RcppExport SEXP _mvcapaCor_which_max_test() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    which_max_test();
    return R_NilValue;
END_RCPP
}
// test_rev_vec
void test_rev_vec();
RcppExport SEXP _mvcapaCor_test_rev_vec() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    test_rev_vec();
    return R_NilValue;
END_RCPP
}
// test_accu
double test_accu();
RcppExport SEXP _mvcapaCor_test_accu() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(test_accu());
    return rcpp_result_gen;
END_RCPP
}
// optimise_savings_vec_old
Rcpp::List optimise_savings_vec_old(const arma::mat& Q, const arma::vec& b, const Rcpp::List& penalty_components, const Rcpp::List& extended_nbs_list);
RcppExport SEXP _mvcapaCor_optimise_savings_vec_old(SEXP QSEXP, SEXP bSEXP, SEXP penalty_componentsSEXP, SEXP extended_nbs_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type penalty_components(penalty_componentsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type extended_nbs_list(extended_nbs_listSEXP);
    rcpp_result_gen = Rcpp::wrap(optimise_savings_vec_old(Q, b, penalty_components, extended_nbs_list));
    return rcpp_result_gen;
END_RCPP
}
// test_grow_tree
void test_grow_tree(const arma::mat& Q, const arma::vec& b, const Rcpp::List& penalty_components, const Rcpp::List& extended_nbs_list);
RcppExport SEXP _mvcapaCor_test_grow_tree(SEXP QSEXP, SEXP bSEXP, SEXP penalty_componentsSEXP, SEXP extended_nbs_listSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type penalty_components(penalty_componentsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type extended_nbs_list(extended_nbs_listSEXP);
    test_grow_tree(Q, b, penalty_components, extended_nbs_list);
    return R_NilValue;
END_RCPP
}
// optimise_savings_old
Rcpp::List optimise_savings_old(const arma::mat& Q, const arma::vec& b, const Rcpp::List& penalty_components, const Rcpp::List& extended_nbs_list);
RcppExport SEXP _mvcapaCor_optimise_savings_old(SEXP QSEXP, SEXP bSEXP, SEXP penalty_componentsSEXP, SEXP extended_nbs_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type penalty_components(penalty_componentsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type extended_nbs_list(extended_nbs_listSEXP);
    rcpp_result_gen = Rcpp::wrap(optimise_savings_old(Q, b, penalty_components, extended_nbs_list));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mvcapaCor_hash_map_test", (DL_FUNC) &_mvcapaCor_hash_map_test, 0},
    {"_mvcapaCor_optimise_savings", (DL_FUNC) &_mvcapaCor_optimise_savings, 4},
    {"_mvcapaCor_optimise_savings_list", (DL_FUNC) &_mvcapaCor_optimise_savings_list, 4},
    {"_mvcapaCor_list_test", (DL_FUNC) &_mvcapaCor_list_test, 1},
    {"_mvcapaCor_seq_int", (DL_FUNC) &_mvcapaCor_seq_int, 2},
    {"_mvcapaCor_test_map", (DL_FUNC) &_mvcapaCor_test_map, 0},
    {"_mvcapaCor_test_ubt", (DL_FUNC) &_mvcapaCor_test_ubt, 1},
    {"_mvcapaCor_which_max_test", (DL_FUNC) &_mvcapaCor_which_max_test, 0},
    {"_mvcapaCor_test_rev_vec", (DL_FUNC) &_mvcapaCor_test_rev_vec, 0},
    {"_mvcapaCor_test_accu", (DL_FUNC) &_mvcapaCor_test_accu, 0},
    {"_mvcapaCor_optimise_savings_vec_old", (DL_FUNC) &_mvcapaCor_optimise_savings_vec_old, 4},
    {"_mvcapaCor_test_grow_tree", (DL_FUNC) &_mvcapaCor_test_grow_tree, 4},
    {"_mvcapaCor_optimise_savings_old", (DL_FUNC) &_mvcapaCor_optimise_savings_old, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_mvcapaCor(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}