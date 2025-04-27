// #include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../base/to_newick.hpp"
#include "../base/to_vector.hpp"
#include "../ops/vector.hpp"

namespace py = pybind11;

PYBIND11_MODULE(phylo2vec, m) {
    m.doc() = "Phylo2Vec: a vector representation for binary trees";

    m.def("to_newick", &toNewick, "Recover a rooted tree(in Newick format) from a Phylo2Vec v");

    m.def("to_vector", &toVector, "Convert a Newick string with parent labels to a vector");

    m.def("to_vector_no_parents", &toVectorNoParents,
          "Convert a Newick string without parent labels to a vector");

    m.def("sample", &sample, "Sample a random Phylo2Vec v for n leaves");

    m.def("check_v", &check_v, "Check that Phylo2Vec v is correct");
}