#include <pybind11/pybind11.h>

#include "base/to_newick.hpp"
#include "base/to_vector.hpp"

namespace py = pybind11;

PYBIND11_MODULE(phylo2vec, m) {
    m.doc() = "Phylo2Vec: a vector representation for binary trees";

    m.add("to_newick", &toNewick,
          "Recover a rooted tree(in Newick format) from a Phylo2Vec v");

    m.add("to_vector", &toVector,
          "Convert a Newick string with parent labels to a vector");

    m.add("to_vector_no_parents", &toVectorNoParents,
          "Convert a Newick string without parent labels to a vector");
}