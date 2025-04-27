#include "pairwise.hpp"

#include <sstream>
#include <stdexcept>

#include "../base/to_newick.hpp"
#include "../ops/vector.hpp"

Matrix copheneticDistances(const PhyloVec &v, bool unrooted) {
    const size_t numLeaves = v.size() + 1;

    Ancestry anc = getAncestry(v);

    if (unrooted) {
        anc[anc.size() - 1][anc.size() - 1] -= 1;
    }

    const size_t numNodes = 2 * numLeaves - 1;

    // Initialize with zeros
    // std::vector<int> fullZeros(numNodes, 0);
    Matrix fullD(numNodes, std::vector<float>(numNodes, 0));
    // for (size_t i = 0; i < numNodes; ++i) {
    //     std::vector<int> row(numNodes, 0);
    //     D.push_back(row);
    // }

    std::vector<int> allVisited;

    for (size_t i = 0; i < numLeaves - 1; ++i) {
        auto &[c1, c2, p] = anc[numLeaves - i - 2];

        for (int j = 0; j < allVisited.size() - 1; ++j) {
            size_t visited = allVisited[j];
            // new distance from c1/c2 to visited is
            // 1 + the distance from p to visited
            float distFromVisited = fullD[p][visited] + 1.0;

            // c1 to visited
            fullD[c1][visited] = distFromVisited;
            fullD[visited][c1] = distFromVisited;

            // c2 to visited
            fullD[c2][visited] = distFromVisited;
            fullD[visited][c2] = distFromVisited;
        }

        // c1 to c2: path length = 2 --> c1 -- p --c2
        fullD[c1][c2] = 2.0;
        fullD[c2][c1] = 2.0;
        // c1 to parent: path length = 1 --> c1 -- p
        fullD[c1][p] = 1.0;
        fullD[p][c1] = 1.0;
        // c2 to parent: path length = 1 --> c2 -- p
        fullD[c2][p] = 1.0;
        fullD[p][c2] = 1.0;

        // all_visited.extend([c1, c2, p])
        allVisited.push_back(c1);
        allVisited.push_back(c2);
        allVisited.push_back(p);
    }

    // leafD = fullD[:numLeaves, :numLeaves]
    Matrix leafD(numLeaves, std::vector<float>(numLeaves, 0));

    for (size_t i = 0; i < numLeaves; ++i) {
        leafD[i] = std::vector<float>(fullD[i].begin(), fullD[i].begin() + numLeaves);
    }

    return leafD;
}

Matrix pairwiseDistances(const PhyloVec &v, std::string_view metric, bool unrooted) {
    if (metric == "cophenetic") {
        return copheneticDistances(v, unrooted);
    } else {
        std::ostringstream oss;
        oss << "Invalid metric name: " << metric;
        throw std::invalid_argument(oss.str());
    }
}