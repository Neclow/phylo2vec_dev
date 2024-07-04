#include <string>
#include <string_view>

#include "../base/core.hpp"

struct Converter {
    Leaf2Taxon mapping;
    std::string intNewick;
};

const std::string startDelimiters = "(,";
const std::string endDelimiters = ",);";

void removeParentLabels(std::string &newick);
void removeBranchLengthAnnotations(std::string &newick);
std::string toStringNewick(Converter converter);
Converter toIntNewick(std::string_view strNewick);
int getNumLeaves(std::string_view newick);