#include "newick.hpp"

void removeParentLabels(std::string &newick) {
    int open_idx = -1;
    const std::string delimiters = ",);";
    for (size_t i = 0; i < newick.size(); i++) {
        if (delimiters.find(newick[i]) != std::string::npos && open_idx != -1) {
            newick.erase(open_idx, i - open_idx);
            i -= i - open_idx;
            open_idx = -1;
        }
        if (newick[i] == ')') {
            open_idx = i + 1;
        }
    }
}

void removeBranchLengthAnnotations(std::string &newick) {
    int open_idx = -1;
    const std::string delimiters = ",);";
    for (size_t i = 0; i < newick.size(); i++) {
        if (delimiters.find(newick[i]) != std::string::npos && open_idx != -1) {
            newick.erase(open_idx, i - open_idx);
            i -= i - open_idx;
            open_idx = -1;
        }
        if (newick[i] == ':') {
            open_idx = i;
        }
    }
}
