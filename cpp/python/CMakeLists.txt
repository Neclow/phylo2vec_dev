add_subdirectory(extern/pybind11)
pybind11_add_module(phylo2vec phylo2vec.cpp)
target_link_libraries(phylo2vec PUBLIC phylo2vec_cpp)
install(TARGETS phylo2vec DESTINATION .)