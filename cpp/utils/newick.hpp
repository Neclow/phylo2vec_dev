#include <string>

#include "../base/core.hpp"

void removeParentLabels(std::string &newick);
void removeBranchLengthAnnotations(std::string &newick);
void applyLeafMapping(std::string &newick, const Leaf2Taxon &mapping);
Leaf2Taxon createLeafMapping(std::string &newick);
int getNumLeaves(const std::string &newick);