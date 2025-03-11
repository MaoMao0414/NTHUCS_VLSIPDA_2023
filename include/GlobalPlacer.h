#ifndef GLOBALPLACER_H
#define GLOBALPLACER_H

#include "Wrapper.hpp"
#include <vector>
#include <utility>
#include <unordered_set>

class GlobalPlacer
{
public:
    GlobalPlacer(wrapper::Placement &placement);

    void overlapPlace();
    void randomPlace(); // An example of random placement implemented by TA
    void place();
    void place2();
    void SAplace(int seedVal);     // use simulated annealing approach to do global placement
    double computeOverlap();
    double computeCost();
    double computeHPWLCost(int moduleIdx);
    double computeTwoModuleOverlap(int module1, int module2);   // 計算兩個module間的overlap量
    double computeModuleOverlap(int moduleIdx);       // 計算單一module的overlap量
    double computeModuleOverlap(int moduleIdx, std::unordered_set<int>& uset_allOverlapModuleIdx);
    void updateBest();

private:
    wrapper::Placement &_placement;
    double _initial_HPWL = -1;
    double _initial_overlap = -1;
    std::vector<std::pair<double, double> > _vec_movableCell_to_bestResult;
};

#endif // GLOBALPLACER_H
