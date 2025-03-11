#ifndef EXAMPLEFUNCTION_H
#define EXAMPLEFUNCTION_H

#include "NumericalOptimizerInterface.h"
#include "Wrapper.hpp"
#include <vector>
#include <map>
#include <set>

class ExampleFunction : public NumericalOptimizerInterface
{
public:
    ExampleFunction(wrapper::Placement &placement);

    void evaluateFG(const vector<double> &x, double &f, vector<double> &g);
    void evaluateF(const vector<double> &x, double &f);
    unsigned dimension();
    void evaluateF_WL(const vector<double> &x, double &f); 
    void evaluateF_WL2(const vector<double> &x, double &f); 
    void evaluateG_WL(const vector<double> &x, vector<double> &g);
    void evaluateF_Density(const vector<double> &x, double &f); 
    void evaluateG_Density(const vector<double> &x, vector<double> &g);

    int _num_blocks;
    double _smooth_parameter;
    std::set<int> _set_movable_cell_idx;
    std::vector<std::pair<int, double>> _vec_pin_x_diff;      // pin i連到的module要(+/-)多少值才能到pin的座標
    std::vector<std::pair<int, double>> _vec_pin_y_diff;  
    std::map<int, std::vector<int>> _map_netIdx_to_vecPinIdx;
    std::map<int ,int> _map_originalIdx_to_efIdx;
    std::vector<std::pair<double, double>> _vec_bin_mid_pt_cord;
    double _wb, _hb;
    double _target_density;
    wrapper::Placement &_placement;

    int _cur_module = -1;    // 當前正在優化的module
    int _beta;
    int _alpha;
};

#endif // EXAMPLEFUNCTION_H
