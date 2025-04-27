#include "newick.hpp"

void removeAnnotations(std::string &newick, const char delimiter,
                       int keepDelimiter) {
    int openIdx = -1;
    for (size_t i = 0; i < newick.size(); ++i) {
        char c = newick[i];
        // c is an end delimiter and the reading frame is open
        // erase what is between ')' and c
        // (i.e., the annotation of interest)
        if (endDelimiters.find(c) != std::string::npos && openIdx != -1) {
            size_t length = i - openIdx;
            newick.erase(openIdx, length);
            // Roll back i to the index before the substring
            i -= length;
            // close the reading frame
            openIdx = -1;
        }
        // c is a delimiter --> open the reading frame
        if (c == delimiter) {
            openIdx = i + keepDelimiter;
        }
    }
}

void removeParentLabels(std::string &newick) {
    const char parentDelimiter = ')';
    removeAnnotations(newick, parentDelimiter, 1);
}

void removeBranchAnnotations(std::string &newick) {
    const char branchDelimiter = ':';
    removeAnnotations(newick, branchDelimiter, 0);
}

Converter toIntNewick(std::string_view strNewick) {
    // Converter = 2 elements: mapping and intNewick
    Converter result;

    int openIdx = -1;
    for (size_t i = 0; i < strNewick.size(); ++i) {
        char c = strNewick[i];
        // char is an end delimiter
        if (endDelimiters.find(c) != std::string::npos && openIdx != -1) {

            // substring between start and end delimiter (= a taxon)
            std::string_view taxon = strNewick.substr(openIdx, i - openIdx);

            // Add the taxon to the mapping
            result.mapping.push_back(taxon);

            // Instead of a taxon, add the next int to the int Newick
            result.intNewick += std::to_string(result.mapping.size());

            // Reset
            openIdx = -1;
        }
        // current char is a start delimiter --> add to the int newick
        if (startDelimiters.find(c) != std::string::npos) {
            openIdx = i + 1;
            result.intNewick += c;
        }

        // open_idx == -1 --> frame is not open --> add the current char
        if (openIdx == -1) {
            result.intNewick += c;
        }
    }

    return result;
}

std::string toStringNewick(Converter converter) {
    // Output
    std::string strNewick;

    int openIdx = -1;

    // Number of leaves/taxa added
    size_t added = 0;

    for (size_t i = 0; i < converter.intNewick.size(); ++i) {
        char c = converter.intNewick[i];
        // char is an end delimiter
        if (endDelimiters.find(c) != std::string::npos && openIdx != -1) {

            // Add the next string to the string Newick
            strNewick += converter.mapping[added];
            ++added;

            // Reset
            openIdx = -1;
        }
        // current char is a start delimiter --> add to the int newick
        if (startDelimiters.find(c) != std::string::npos) {
            openIdx = i + 1;
            strNewick += c;
        }

        // open_idx == -1 --> frame is not open --> add the current char
        if (openIdx == -1) {
            strNewick += c;
        }
    }

    return strNewick;
}