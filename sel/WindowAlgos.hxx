//
// Created by ecorona on 3/24/16.
//

#ifndef POSSEL_WINDOWALGOS_HXX
#define POSSEL_WINDOWALGOS_HXX

#include <vector>
#include <algorithm>
namespace win{

    auto getIndexBehindOrEqual(std::vector<int>& pos, int chrPos)
    {
        auto it = std::upper_bound(pos.begin(), pos.end(), chrPos); //---*|---- or //----|-----
        long startI = std::distance(pos.begin(), it);
        return startI == 0 ? startI : startI-1;
    }

    auto getIndexInfront(std::vector<int>& pos, int chrPos)
    {
        auto it = std::upper_bound(pos.begin(), pos.end(), chrPos); //---*|---- or //----|-----
        long startI = std::distance(pos.begin(), it);
        return startI == static_cast<int>(pos.size()) ? startI-1 : startI;
    }
}
#endif //POSSEL_WINDOWALGOS_HXX
