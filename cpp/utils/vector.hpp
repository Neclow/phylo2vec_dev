#include "../base/core.hpp"

void reorder(PhyloVec &v, Leaf2Taxon &mapping, std::string_view method);
void reorderBirthDeath(Ancestry &abort, Leaf2Taxon &mapping,
                       bool reorderInternal = true, bool shuffleCols = true);
void rerootAtRandom(PhyloVec &v);
unsigned int removeLeaf(PhyloVec &v, const unsigned int &leaf);
void addLeaf(PhyloVec &v, const unsigned int &leaf, const unsigned int &pos);
std::pair<size_t, size_t> findCoordsOfFirstLeaf(const Ancestry &ancestry,
                                                const int &leaf);