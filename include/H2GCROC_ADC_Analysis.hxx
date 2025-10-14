#include "H2GCROC_Common.hxx"

// the average of the smallest 2 samples among the first 3 samples
inline void pedestal_subtraction_2minoffirst3(std::vector<int>& data, double& pedestal) {
    if (data.size() < 3) {
        pedestal = 0.0;
        return;
    }
    std::vector<int> first3(data.begin(), data.begin() + 3);
    std::sort(first3.begin(), first3.end());
    pedestal = (double(first3[0]) + double(first3[1])) / 2.0;
    for (auto& val : data) {
        val -= pedestal;
        if (val < 0) {
            val = 0;
        }
    }
}

