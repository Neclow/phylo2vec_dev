#include <string>
#include <string_view>

#include "../base/core.hpp"

const std::string startDelimiters = "(,";
const std::string endDelimiters = ",);";

void removeAnnotations(std::string &newick, const char delimiter,
                       const bool keepDelimiter);
void removeParentLabels(std::string &newick);
void removeBranchAnnotations(std::string &newick);
std::string toStringNewick(Converter converter);
Converter toIntNewick(std::string_view strNewick);
int getNumLeaves(std::string_view newick);