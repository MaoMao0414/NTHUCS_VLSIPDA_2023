#include "ExampleFunction.h"
#include "bits/stdc++.h"
#include <iomanip>

// minimize 3*x^2 + 2*x*y + 2*y^2 + 7

ExampleFunction::ExampleFunction(wrapper::Placement &placement) : _placement(placement)
{

}

void ExampleFunction::evaluateF_WL(const vector<double> &x, double &f)
{
    double lse_wl = 0;
    // cell_i x 座標: i*2; y 座標: i*2 + 1 
    for(auto piv : _map_netIdx_to_vecPinIdx)
    {
        vector<int> vec_pin_idx = piv.second;

        /*****求max{xi}部分**********/
        double x_max = -DBL_MAX;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_x_diff[pin_idx].first;
            double bias = _vec_pin_x_diff[pin_idx].second;
            
            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2] + bias;

            x_max = max(x_max, cord);
        }

        double exp_sum_x = 0;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_x_diff[pin_idx].first;
            double bias = _vec_pin_x_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2] + bias;

            exp_sum_x += exp((cord - x_max)/_smooth_parameter);
        }

        lse_wl = lse_wl + x_max + log(exp_sum_x)*_smooth_parameter;
        /******************************/

        /*****求max{-xi}部分**********/
        double x_max_neg = -DBL_MAX;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_x_diff[pin_idx].first;
            double bias = _vec_pin_x_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2] + bias;

            x_max_neg = max(x_max_neg, -cord);
        }

        double exp_sum_xneg = 0;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_x_diff[pin_idx].first;
            double bias = _vec_pin_x_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2] + bias;

            exp_sum_xneg += exp((-cord - x_max_neg)/_smooth_parameter);
        }

        lse_wl = lse_wl + x_max_neg + log(exp_sum_xneg)*_smooth_parameter;
        /******************************/

        /*****求max{yi}部分**********/
        double y_max = -DBL_MAX;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_y_diff[pin_idx].first;
            double bias = _vec_pin_y_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2 + 1] + bias;

            y_max = max(y_max, cord);
        }

        double exp_sum_y = 0;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_y_diff[pin_idx].first;
            double bias = _vec_pin_y_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2 + 1] + bias;

            //cout << cord << " " << y_max << endl;
            //cout << exp(cord - y_max) << endl;
            exp_sum_y += exp((cord - y_max)/_smooth_parameter);
        }

        //cout << exp_sum_y << endl;
        lse_wl = lse_wl + y_max + log(exp_sum_y)*_smooth_parameter;
        /******************************/

        /*****求max{-yi}部分**********/
        double y_max_neg = -DBL_MAX;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_y_diff[pin_idx].first;
            double bias = _vec_pin_y_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2 + 1] + bias;

            y_max_neg = max(y_max_neg, -cord);
        }

        //cout << y_max_neg << endl;
        double exp_sum_y_neg = 0;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_y_diff[pin_idx].first;
            double bias = _vec_pin_y_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2 + 1] + bias;

            exp_sum_y_neg += exp((-cord - y_max_neg)/_smooth_parameter);
        }

        //cout << exp_sum_y_neg << endl;
        lse_wl = lse_wl + y_max_neg + log(exp_sum_y_neg)*_smooth_parameter;
        /******************************/
    }

    f += _alpha*lse_wl;
}

void ExampleFunction::evaluateG_WL(const vector<double> &x, vector<double> &g)
{
    for(auto piv : _map_netIdx_to_vecPinIdx)
    {
        vector<int> vec_pin_idx = piv.second;

        /*****計算x坐標的max和exp_sum**********/
        double x_max = -DBL_MAX;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_x_diff[pin_idx].first;
            double bias = _vec_pin_x_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2] + bias;

            x_max = max(x_max, cord);
        }

        double exp_sum_x = 0;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_x_diff[pin_idx].first;
            double bias = _vec_pin_x_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2] + bias;

            exp_sum_x += exp((cord - x_max)/_smooth_parameter);
        }

        /*****計算x坐標的max_neg和exp_sum_neg**********/
        double x_max_neg = -DBL_MAX;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_x_diff[pin_idx].first;
            double bias = _vec_pin_x_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2] + bias;

            x_max_neg = max(x_max_neg, -cord);
        }

        double exp_sum_xneg = 0;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_x_diff[pin_idx].first;
            double bias = _vec_pin_x_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2] + bias;

            exp_sum_xneg += exp((-cord - x_max_neg)/_smooth_parameter);
        }

        /*****計算y坐標的max和exp_sum**********/
        double y_max = -DBL_MAX;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_y_diff[pin_idx].first;
            double bias = _vec_pin_y_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2 + 1] + bias;

            y_max = max(y_max, cord);
        }

        double exp_sum_y = 0;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_y_diff[pin_idx].first;
            double bias = _vec_pin_y_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2 + 1] + bias;

            exp_sum_y += exp((cord - y_max)/_smooth_parameter);
        }

        /*****計算y坐標的max_neg和exp_sum_neg**********/
        double y_max_neg = -DBL_MAX;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_y_diff[pin_idx].first;
            double bias = _vec_pin_y_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2 + 1] + bias;

            y_max_neg = max(y_max_neg, -cord);
        }

        double exp_sum_yneg = 0;
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_y_diff[pin_idx].first;
            double bias = _vec_pin_y_diff[pin_idx].second;

            double cord = (connect_module_idx == -1) ? bias : x[connect_module_idx*2 + 1] + bias;

            exp_sum_yneg += exp((-cord - y_max_neg)/_smooth_parameter);
        }

        // 計算x坐標的梯度
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_x_diff[pin_idx].first;
            if (connect_module_idx != -1)
            {
                double bias = _vec_pin_x_diff[pin_idx].second;
                double cord = x[connect_module_idx*2] + bias;

                // 對於x_max的梯度
                double grad_x_max = exp((cord - x_max)/_smooth_parameter) / exp_sum_x;

                // 對於x_max_neg的梯度
                double grad_x_max_neg = -exp((-cord - x_max_neg)/_smooth_parameter) / exp_sum_xneg;

                // 將這兩部分的梯度加到總梯度上
                g[connect_module_idx*2] += (grad_x_max + grad_x_max_neg) * _smooth_parameter * _alpha;
            }
        }

        // 計算y坐標的梯度
        for(int pin_idx : vec_pin_idx)
        {
            int connect_module_idx = _vec_pin_y_diff[pin_idx].first;
            if (connect_module_idx != -1)
            {
                double bias = _vec_pin_y_diff[pin_idx].second;
                double cord = x[connect_module_idx*2 + 1] + bias;

                // 對於y_max的梯度
                double grad_y_max = exp((cord - y_max)/_smooth_parameter) / exp_sum_y;

                // 對於y_max_neg的梯度
                double grad_y_max_neg = -exp((-cord - y_max_neg)/_smooth_parameter) / exp_sum_yneg;

                // 將這兩部分的梯度加到總梯度上
                g[connect_module_idx*2 + 1] += (grad_y_max + grad_y_max_neg) * _smooth_parameter * _alpha;
            }
        }
    }

}

void ExampleFunction::evaluateF_Density(const vector<double> &x, double &f)
{
    double totalDensity = 0;

    unordered_map<int, vector<pair<int, double>>> umap_cellInBin_coverArea;     // cell i 在 bin i 中的覆蓋量

    for(int binIdx = 0; binIdx < _vec_bin_mid_pt_cord.size(); binIdx++)
    {
        double bin_x = _vec_bin_mid_pt_cord[binIdx].first;
        double bin_y = _vec_bin_mid_pt_cord[binIdx].second;

        for(int moduleIdx = 0; moduleIdx < _placement.numModules(); moduleIdx++)
        {
            //cout << moduleIdx << endl;
            double module_x, module_y;
            
            if(_placement.module(moduleIdx).isFixed() || (_cur_module != -1 && moduleIdx != _cur_module))
            {
                module_x = _placement.module(moduleIdx).centerX();
                module_y = _placement.module(moduleIdx).centerY();                
            }
            else
            {
                int ef_idx = _map_originalIdx_to_efIdx[moduleIdx];

                module_x = x[ef_idx*2];
                module_y = x[ef_idx*2 + 1];
            }
            
            double wi = _placement.module(moduleIdx).width();
            double hi = _placement.module(moduleIdx).height();

            // 先來算x方向的重疊量
            double x_overlap_area = 0;
            
            double dx = abs(module_x - bin_x);
            double a_x = 4 / ((_wb + wi) * (2*_wb + wi));
            double b_x = 4 / (_wb * (2*_wb + wi));

            if(0 <= dx && dx <= (_wb/2) + (wi/2))
                x_overlap_area = 1 - a_x * dx * dx;
            else if((_wb/2) + (wi/2) <= dx && dx <= _wb + (wi/2))
                x_overlap_area = b_x * (dx - _wb - wi/2) * (dx - _wb - wi/2);
            else
                x_overlap_area = 0;

            // 再來算x方向的重疊量
            double y_overlap_area = 0;
            
            double dy = abs(module_y - bin_y);
            double a_y = 4 / ((_hb + hi) * (2*_hb + hi));
            double b_y = 4 / (_hb * (2*_hb + hi));

            if(0 <= dy && dy <= (_hb/2) + (hi/2))
                y_overlap_area = 1 - a_y * dy * dy;
            else if((_hb/2) + (hi/2) <= dy && dy <= _hb + (hi/2))
                y_overlap_area = b_y * (dy - _hb - hi/2) * (dy - _hb - hi/2);
            else
                y_overlap_area = 0;
        
            umap_cellInBin_coverArea[moduleIdx].push_back({binIdx, x_overlap_area*y_overlap_area});
        }     
    }

    // 計算c
    vector<double> vec_c(_placement.numModules());
    for(int moduleIdx = 0; moduleIdx < vec_c.size(); moduleIdx++)
    {
        double totalArea = 0;
        for(auto pid : umap_cellInBin_coverArea[moduleIdx])
            totalArea += pid.second;
        
        vec_c[moduleIdx] = _placement.module(moduleIdx).area() / totalArea;
    }

    vector<double> vec_binDensity(_vec_bin_mid_pt_cord.size());
    for(int moduleIdx = 0; moduleIdx < vec_c.size(); moduleIdx++)
    {
        for(auto pid : umap_cellInBin_coverArea[moduleIdx])
            vec_binDensity[pid.first] += vec_c[moduleIdx] * pid.second;
    }

    for(double den : vec_binDensity)
        totalDensity += (den - _target_density) * (den - _target_density);

    // cout << "den : " << endl;
    // for(double den : vec_binDensity)
    //     cout << den << endl;


    f += _beta*totalDensity;
}

void ExampleFunction::evaluateG_Density(const vector<double> &x, vector<double> &g) {

    // Temporary data structure to store the overlap areas for each cell in each bin
    unordered_map<int, vector<pair<int, double>>> umap_cellInBin_coverArea;

    for(int binIdx = 0; binIdx < _vec_bin_mid_pt_cord.size(); binIdx++)
    {
        double bin_x = _vec_bin_mid_pt_cord[binIdx].first;
        double bin_y = _vec_bin_mid_pt_cord[binIdx].second;

        for(int moduleIdx = 0; moduleIdx < _placement.numModules(); moduleIdx++)
        {
            //cout << moduleIdx << endl;
            double module_x, module_y;

            if(_placement.module(moduleIdx).isFixed() || (_cur_module != -1 && moduleIdx != _cur_module))
            {
                module_x = _placement.module(moduleIdx).centerX();
                module_y = _placement.module(moduleIdx).centerY();                
            }
            else
            {
                int ef_idx = _map_originalIdx_to_efIdx[moduleIdx];

                module_x = x[ef_idx*2];
                module_y = x[ef_idx*2 + 1];
            }
            
            double wi = _placement.module(moduleIdx).width();
            double hi = _placement.module(moduleIdx).height();

            // 先來算x方向的重疊量
            double x_overlap_area = 0;
            
            double dx = abs(module_x - bin_x);
            double a_x = 4 / ((_wb + wi) * (2*_wb + wi));
            double b_x = 4 / (_wb * (2*_wb + wi));

            if(0 <= dx && dx <= (_wb/2) + (wi/2))
                x_overlap_area = 1 - a_x * dx * dx;
            else if((_wb/2) + (wi/2) <= dx && dx <= _wb + (wi/2))
                x_overlap_area = b_x * (dx - _wb - wi/2) * (dx - _wb - wi/2);
            else
                x_overlap_area = 0;

            // 再來算x方向的重疊量
            double y_overlap_area = 0;
            
            double dy = abs(module_y - bin_y);
            double a_y = 4 / ((_hb + hi) * (2*_hb + hi));
            double b_y = 4 / (_hb * (2*_hb + hi));

            if(0 <= dy && dy <= (_hb/2) + (hi/2))
                y_overlap_area = 1 - a_y * dy * dy;
            else if((_hb/2) + (hi/2) <= dy && dy <= _hb + (hi/2))
                y_overlap_area = b_y * (dy - _hb - hi/2) * (dy - _hb - hi/2);
            else
                y_overlap_area = 0;
        
            umap_cellInBin_coverArea[moduleIdx].push_back({binIdx, x_overlap_area*y_overlap_area});
        }     
    }

    // 計算c
    vector<double> vec_c(_placement.numModules());
    for(int moduleIdx = 0; moduleIdx < vec_c.size(); moduleIdx++)
    {
        double totalArea = 0;
        for(auto pid : umap_cellInBin_coverArea[moduleIdx])
            totalArea += pid.second;
        
        vec_c[moduleIdx] = _placement.module(moduleIdx).area() / totalArea;
    }

    // Now, calculate the gradient
    for (int moduleIdx = 0; moduleIdx < _placement.numModules(); moduleIdx++) {
        // Skip if module is fixed
        if(_cur_module != -1 && moduleIdx != _cur_module)
            continue;
        
        if (!_placement.module(moduleIdx).isFixed()) {
            int ef_idx = _map_originalIdx_to_efIdx[moduleIdx];
            double module_x = x[ef_idx * 2];
            double module_y = x[ef_idx * 2 + 1];
            double wi = _placement.module(moduleIdx).width();
            double hi = _placement.module(moduleIdx).height();

            for (int binIdx = 0; binIdx < _vec_bin_mid_pt_cord.size(); binIdx++) {
                double bin_x = _vec_bin_mid_pt_cord[binIdx].first;
                double bin_y = _vec_bin_mid_pt_cord[binIdx].second;
                double dx = module_x - bin_x;
                double dy = module_y - bin_y;
                double abs_dx = abs(dx);
                double abs_dy = abs(dy);

                // Compute constants 'a' and 'b' for both x and y
                double a_x = 4 / ((_wb + wi) * (2 * _wb + wi));
                double b_x = 4 / (_wb * (2 * _wb + wi));
                double a_y = 4 / ((_hb + hi) * (2 * _hb + hi));
                double b_y = 4 / (_hb * (2 * _hb + hi));

                // Calculate the gradient for x and y coordinates
                double grad_x = 0.0;
                double grad_y = 0.0;

                // Calculate the gradient for the x coordinate
                if (0 <= abs_dx && abs_dx <= (_wb/2) + (wi/2)) {
                    grad_x = -2 * a_x * dx;
                } else if ((_wb/2) + (wi/2) < abs_dx && abs_dx <= _wb + (wi/2)) {
                    grad_x = 2 * b_x * (dx - _wb - wi/2) * (dx >= 0 ? 1 : -1);
                }

                // Calculate the gradient for the y coordinate
                if (0 <= abs_dy && abs_dy <= (_hb/2) + (hi/2)) {
                    grad_y = -2 * a_y * dy;
                } else if ((_hb/2) + (hi/2) < abs_dy && abs_dy <= _hb + (hi/2)) {
                    grad_y = 2 * b_y * (dy - _hb - hi/2) * (dy >= 0 ? 1 : -1);
                }

                // Multiply by the normalization factor and the current density difference
                double density_diff = vec_c[moduleIdx] * umap_cellInBin_coverArea[moduleIdx][binIdx].second - _target_density;
                g[ef_idx * 2] += 2 * density_diff * grad_x * _beta;
                g[ef_idx * 2 + 1] += 2 * density_diff * grad_y * _beta;
            }
        }
    }

}

void ExampleFunction::evaluateFG(const vector<double> &x, double &f, vector<double> &g)
{
    f = 0;
    for(int i = 0; i < g.size(); i++)
        g[i] = 0;

    evaluateF_WL(x, f);
    evaluateF_Density(x, f);
    
    evaluateG_WL(x, g);
    evaluateG_Density(x, g);
}


void ExampleFunction::evaluateF(const vector<double> &x, double &f)
{   
    f = 0;
    evaluateF_WL(x, f);
    evaluateF_Density(x, f);

    //printf("@@@@@@@@@@@@@@@@@@@@ final objective(wirelength part): %.0f @@@@@@@@@@@@@@@@@@@@n", f);
}

unsigned ExampleFunction::dimension()
{
    return _num_blocks*2; // num_blocks*2
    // each two dimension represent the X and Y dimensions of each block
}
