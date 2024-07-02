#include "base/core.hpp"
#include "base/to_newick.hpp"
#include "utils/random.hpp"

int main(int argc, char *argv[]) {
    PhyloVec v = sample(200, false);

    std::string newick = toNewick(v);

    return 0;
}