#include "../base/core.hpp"

PhyloVec sample(const size_t &numLeaves, bool ordered = false);
void check_v(const PhyloVec &v);
void reorder(PhyloVec &v, Leaf2Taxon &mapping, std::string_view method);
void reorderBirthDeath(Ancestry &abort, Leaf2Taxon &mapping, bool reorderInternal = true,
                       bool shuffleCols = true);
void rerootAtRandom(PhyloVec &v);
unsigned int removeLeaf(PhyloVec &v, unsigned int leaf);
void addLeaf(PhyloVec &v, unsigned int leaf, unsigned int pos);
std::pair<size_t, size_t> findCoordsOfFirstLeaf(const Ancestry &ancestry, int leaf);
std::vector<std::vector<unsigned int>> getAncestryPaths(const PhyloVec &v);
int getCommonAncestor(const PhyloVec &v, unsigned int node1, unsigned int node2);