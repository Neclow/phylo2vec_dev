#include "validation.hpp"

#include <sstream>
#include <stdexcept>

void check_v(const PhyloVec &v) {
    // check that v is valid: 0 <= v[i] <= 2i

    const size_t k = v.size();
    for (size_t i = 0; i < k; ++i) {
        if (v[i] > 2 * i) {
            std::ostringstream oss;
            oss << "Invalid value at index " << i
                << ": v[i] should be less than 2i, found " << v[i] << ".";
            throw std::out_of_range(oss.str());
        }
    }
}