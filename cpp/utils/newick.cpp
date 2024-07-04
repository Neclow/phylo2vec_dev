#include "newick.hpp"

void removeParentLabels(std::string &newick) {
    int open_idx = -1;
    const std::string delimiters = ",);";
    for (size_t i = 0; i < newick.size(); ++i) {
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
    for (size_t i = 0; i < newick.size(); ++i) {
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

Converter toIntNewick(std::string_view strNewick) {
    Converter result;

    int open_idx = -1;
    for (size_t i = 0; i < strNewick.size(); ++i) {
        // char is an end delimiter
        if (endDelimiters.find(strNewick[i]) != std::string::npos &&
            open_idx != -1) {

            // substring between start and end delimiter
            std::string_view subnewick =
                strNewick.substr(open_idx, i - open_idx);

            // Add the next int to the int Newick
            result.intNewick += std::to_string(result.mapping.size());

            // Fill mapping
            result.mapping.push_back(subnewick);

            // Reset
            open_idx = -1;
        }
        // current char is a start delimiter
        if (startDelimiters.find(strNewick[i]) != std::string::npos) {
            open_idx = i + 1;
            // Add to the int newick
            result.intNewick += strNewick[i];
        }

        // open_idx == -1 --> frame is not open --> add the current char
        if (open_idx == -1) {
            result.intNewick += strNewick[i];
        }
    }

    return result;
}

std::string toStringNewick(Converter converter) {
    std::string strNewick;

    int open_idx = -1;

    size_t j = 0;

    for (size_t i = 0; i < converter.intNewick.size(); ++i) {
        // char is an end delimiter
        if (endDelimiters.find(converter.intNewick[i]) != std::string::npos &&
            open_idx != -1) {

            // Add the next int to the int Newick
            strNewick += converter.mapping[j];

            ++j;

            // Reset
            open_idx = -1;
        }
        // current char is a start delimiter
        if (startDelimiters.find(converter.intNewick[i]) != std::string::npos) {
            open_idx = i + 1;
            // Add to the int newick
            strNewick += converter.intNewick[i];
        }

        // open_idx == -1 --> frame is not open --> add the current char
        if (open_idx == -1) {
            strNewick += converter.intNewick[i];
        }
    }

    return strNewick;
}