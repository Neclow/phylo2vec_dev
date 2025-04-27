#include <Rcpp.h>

#include "../base/to_newick.hpp"
#include "../base/to_vector.hpp"
#include "../ops/vector.hpp"

using namespace Rcpp;

// [[Rcpp::plugins(cpp17)]]

// [[Rcpp::export]]
std::string to_newick(const PhyloVec &v) { return toNewick(v); }

// [[Rcpp::export]]
PhyloVec to_vector(const std::string &newick) {
    // Rcpp doesn't like std::string_view
    std::string newick_str{newick};
    return toVector(newick_str);
}

// [[Rcpp::export]]
PhyloVec to_vector_no_parents(const std::string &newick) {
    std::string newick_str{newick};
    return toVectorNoParents(newick_str);
}

// [[Rcpp::export(name = "check")]]
PhyloVec check_p2v(const PhyloVec &v) { check_v(v); }

// [[Rcpp::export]]
void sample_p2v(const size_t &numLeaves, bool ordered) { sample(numLeaves, ordered); }