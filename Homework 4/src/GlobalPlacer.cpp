#include "GlobalPlacer.h"
#include "ExampleFunction.h"
#include "NumericalOptimizer.h"

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <chrono>
#include "bits/stdc++.h"
using namespace std;

struct Interval {
    double start, end;
    int moduleIdx;

    Interval(double s, double e, int idx) : start(s), end(e), moduleIdx(idx) {}
};

struct TreeNode {
    Interval* interval;
    double maxEnd;
    TreeNode *left, *right, *parent;

    TreeNode(Interval* interval) 
        : interval(interval), maxEnd(interval->end), left(nullptr), right(nullptr), parent(nullptr) {}
};

class IntervalTree {
private:
    TreeNode* root;
    std::map<int, TreeNode*> moduleIdxToNode;

    TreeNode* insert(TreeNode* node, Interval* interval) {
        if (!node) {
            node = new TreeNode(interval);
            moduleIdxToNode[interval->moduleIdx] = node;
            return node;
        }

        double start = interval->start;
        if (start < node->interval->start) {
            node->left = insert(node->left, interval);
            node->left->parent = node;
        } else {
            node->right = insert(node->right, interval);
            node->right->parent = node;
        }

        node->maxEnd = std::max(node->maxEnd, interval->end);
        return node;
    }

    void queryOverlaps(TreeNode* node, Interval& interval, std::unordered_set<int>& overlapIndices) {
        if (!node) return;

        if (node->interval->start < interval.end && node->interval->end > interval.start) {
            overlapIndices.insert(node->interval->moduleIdx);
        }

        if (node->left && node->left->maxEnd > interval.start) {
            queryOverlaps(node->left, interval, overlapIndices);
        }

        if (node->right) {
            queryOverlaps(node->right, interval, overlapIndices);
        }
    }

public:
    IntervalTree() : root(nullptr) {}

    void insertInterval(double start, double end, int moduleIdx) {
        Interval* interval = new Interval(start, end, moduleIdx);
        root = insert(root, interval);
    }

    void modifyInterval(int moduleIdx, double newStart, double newEnd) {
        if (moduleIdxToNode.find(moduleIdx) == moduleIdxToNode.end()) {
            std::cout << "ModuleIdx not found" << std::endl;
            return;
        }

        TreeNode* node = moduleIdxToNode[moduleIdx];
        node->interval->start = newStart;
        node->interval->end = newEnd;

        while (node != nullptr) {
            node->maxEnd = std::max(node->interval->end, 
                                    std::max(node->left ? node->left->maxEnd : 0, 
                                             node->right ? node->right->maxEnd : 0));
            node = node->parent;
        }
    }

    void queryOverlaps(double start, double end, std::unordered_set<int>& overlapIndices) {
        Interval queryInterval(start, end, -1);

        queryOverlaps(root, queryInterval, overlapIndices);
    }
};

GlobalPlacer::GlobalPlacer(wrapper::Placement &placement)
    : _placement(placement)
{
}

void GlobalPlacer::overlapPlace()
{
    // 初始重疊比率
    double initialOverlapRatioWidth = 0.05;
    double initialOverlapRatioHeight = 0.05;

    // 當需要更多空間時增加的重疊比率
    double additionalOverlapRatio = 0.05;

    double currentOverlapRatioWidth = initialOverlapRatioWidth;
    double currentOverlapRatioHeight = initialOverlapRatioHeight;

    double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
    double coreHeight = _placement.boundryTop() - _placement.boundryBottom();

    double currentX = _placement.boundryLeft();
    double currentY = _placement.boundryBottom();
    double maxHeightInRow = 0;

    for (size_t i = 0; i < _placement.numModules(); ++i)
    {
        double width = _placement.module(i).width();
        double height = _placement.module(i).height();

        // 考慮固定模塊
        if (_placement.module(i).isFixed())
        {
            // 更新當前行的最大高度
            maxHeightInRow = max(maxHeightInRow, height);
            continue;
        }

        double effectiveWidth = width * (1 - currentOverlapRatioWidth);
        double effectiveHeight = height * (1 - currentOverlapRatioHeight);

        if (currentX + effectiveWidth > _placement.boundryRight())
        {
            currentX = _placement.boundryLeft();
            currentY += maxHeightInRow;
            maxHeightInRow = 0;
            currentOverlapRatioWidth += additionalOverlapRatio; // 增加重疊比率
            currentOverlapRatioHeight += additionalOverlapRatio; // 增加重疊比率
        }

        if (currentY + effectiveHeight > _placement.boundryTop())
        {
            currentY = _placement.boundryBottom(); // 從頂部重新開始
            currentOverlapRatioWidth += additionalOverlapRatio; // 增加重疊比率
            currentOverlapRatioHeight += additionalOverlapRatio; // 增加重疊比率
        }

        _placement.module(i).setPosition(currentX, currentY);
        currentX += effectiveWidth;
        maxHeightInRow = std::max(maxHeightInRow, effectiveHeight);
    }

    cout << "end of initial placement" << endl;
}

void GlobalPlacer::randomPlace()
{
    double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
    double coreHeight = _placement.boundryTop() - _placement.boundryBottom();
    for (size_t i = 0; i < _placement.numModules(); ++i)
    {
        if (_placement.module(i).isFixed())
            continue;

        double width = _placement.module(i).width();
        double height = _placement.module(i).height();
        double x = rand() % static_cast<int>(coreWidth - width) + _placement.boundryLeft();
        double y = rand() % static_cast<int>(coreHeight - height) + _placement.boundryBottom();
        _placement.module(i).setPosition(x, y);
    }
}

double GlobalPlacer::computeOverlap()
{
    double overlap_area = 0;

    for (size_t i = 0; i < _placement.numModules(); i++)
    {
        double module_1_width = _placement.module(i).width();
        double module_1_height = _placement.module(i).height();

        double module_1_center_x = _placement.module(i).centerX();
        double module_1_center_y = _placement.module(i).centerY();

        for (size_t j = i+1; j < _placement.numModules(); j++)
        {
            double module_2_width = _placement.module(j).width();
            double module_2_height = _placement.module(j).height();

            double module_2_center_x = _placement.module(j).centerX();
            double module_2_center_y = _placement.module(j).centerY();

            double x_overlap = max(0.0, min(module_1_center_x + module_1_width/2, module_2_center_x + module_2_width/2) - max(module_1_center_x - module_1_width/2, module_2_center_x - module_2_width/2));
            double y_overlap = max(0.0, min(module_1_center_y + module_1_height/2, module_2_center_y + module_2_height/2) - max(module_1_center_y - module_1_height/2, module_2_center_y - module_2_height/2));

            overlap_area = overlap_area + x_overlap*y_overlap;
            //cout << overlap_area << endl;
        }
    }

    return overlap_area;
}

double GlobalPlacer::computeTwoModuleOverlap(int module1, int module2)
{
    double module_1_width = _placement.module(module1).width();
    double module_1_height = _placement.module(module1).height();

    double module_1_center_x = _placement.module(module1).centerX();
    double module_1_center_y = _placement.module(module1).centerY();   

    double module_2_width = _placement.module(module2).width();
    double module_2_height = _placement.module(module2).height();

    double module_2_center_x = _placement.module(module2).centerX();
    double module_2_center_y = _placement.module(module2).centerY();

    double x_overlap = max(0.0, min(module_1_center_x + module_1_width/2, module_2_center_x + module_2_width/2) - max(module_1_center_x - module_1_width/2, module_2_center_x - module_2_width/2));
    double y_overlap = max(0.0, min(module_1_center_y + module_1_height/2, module_2_center_y + module_2_height/2) - max(module_1_center_y - module_1_height/2, module_2_center_y - module_2_height/2));

    return x_overlap*y_overlap;
}

double GlobalPlacer::computeModuleOverlap(int moduleIdx)
{
    double totalOverlapArea = 0;

    for(int calModuleIdx = 0; calModuleIdx < _placement.numModules(); calModuleIdx++)
    {
        if(calModuleIdx == moduleIdx)
            continue;
        
        totalOverlapArea += computeTwoModuleOverlap(moduleIdx, calModuleIdx);
    }

    return totalOverlapArea;
}

double GlobalPlacer::computeModuleOverlap(int moduleIdx, std::unordered_set<int>& uset_allOverlapModuleIdx)
{
    double totalOverlapArea = 0;

    for(auto m_idx : uset_allOverlapModuleIdx)
    {
        if(m_idx == moduleIdx)
            continue;

        totalOverlapArea += computeTwoModuleOverlap(moduleIdx, m_idx);
    }

    return totalOverlapArea;
}

double GlobalPlacer::computeCost()
{
    double current_HPWL = _placement.computeHpwl();
    double current_overlap = computeOverlap();

    if(_initial_HPWL == -1)
        _initial_HPWL = current_HPWL;
    
    if(_initial_overlap == -1)
        _initial_overlap = current_overlap;
    
    // TODO:
    // cost function = aC1 + bC2
    double HPWL_weight = 1;         // hyperparameter
    double overlap_weight = 0;      // hyperparameter

    return HPWL_weight*(current_HPWL/_initial_HPWL) + overlap_weight*(current_overlap/_initial_overlap)*(current_overlap/_initial_overlap);
}

double GlobalPlacer::computeHPWLCost(int moduleIdx)
{
    double hpwlCost = 0;
    
    // 這個i是迴圈counter和pin的index不同
    for(int i = 0; i < _placement.module(moduleIdx).numPins(); i++)
    {
        // 計算這些pin連到的net的HPWL
        int netIdx = _placement.module(moduleIdx).pin(i).netId();
        
        // 每個net的rectengular box
        double leftBound = DBL_MAX;
        double rightBound = -DBL_MAX;

        double bottomBound = DBL_MAX;
        double topBound = -DBL_MAX;

        // j 同理
        for(int j = 0; j < _placement.net(netIdx).numPins(); j++)
        {
            double pinCordX = _placement.net(netIdx).pin(j).x();
            double pinCordY = _placement.net(netIdx).pin(j).y();   

            leftBound = min(leftBound, pinCordX);
            rightBound = max(rightBound, pinCordX);

            bottomBound = min(bottomBound, pinCordY);
            topBound = max(topBound, pinCordY);
        }

        hpwlCost += (rightBound - leftBound) + (topBound - bottomBound);
    }

    return hpwlCost;
}

void GlobalPlacer::updateBest()
{
    _vec_movableCell_to_bestResult.clear();

    for (size_t i = 0; i < _placement.numModules(); ++i)
    {
        double center_x = _placement.module(i).centerX();
        double center_y = _placement.module(i).centerY();

        _vec_movableCell_to_bestResult.push_back(make_pair(center_x, center_y));
    }

    return;
}

void GlobalPlacer::SAplace(int seedVal)
{
    srand(seedVal);

    auto start_time = chrono::high_resolution_clock::now();

    //overlapPlace();
    randomPlace();
    
    // 建立兩顆interval tree用來記錄每個module x, y軸的區間
    // IntervalTree tree_xDir, tree_yDir;
    // for(size_t moduleIdx = 0; moduleIdx < _placement.numModules(); moduleIdx++)
    // {
    //     double moduleLeftx = _placement.module(moduleIdx).x();
    //     double moduleLowy =  _placement.module(moduleIdx).y();

    //     double moduleRightx = moduleLeftx + _placement.module(moduleIdx).width();
    //     double moduleTopy = moduleLowy + _placement.module(moduleIdx).height();

    //     tree_xDir.insertInterval(moduleLeftx, moduleRightx, moduleIdx);
    //     tree_yDir.insertInterval(moduleLowy, moduleTopy, moduleIdx);
    // }

    double initialWL = 0;
    //double initialOverlap = computeOverlap();
    double initialOverlap = 0;
    if(initialOverlap == 0)
        initialOverlap = 1;

    double currentWL = initialWL;
    double currentOverlap = initialOverlap;

    double coreWidth = _placement.boundryRight() - _placement.boundryLeft();
    double coreHeight = _placement.boundryTop() - _placement.boundryBottom();

    double totalModuleArea = 0;
    for(int i = 0; i < _placement.numModules(); i++)
        totalModuleArea += _placement.module(i).area();
    
    double maxHPWL = coreWidth + coreHeight;
    
    // annealing schedule
    double initialTemperature = DBL_MAX;
    double decreaseRate = 0.995;
    double terminate_temperature = 1e-150;
    double fixBoundRate = 0.99995;

    double greedyBound = 1e-4;    // TODO                                                                                   
    double fixBound = 1e-6;  

    double N = 128;     // 每個溫度產生的解數目 // TODO

    // cost function weight
    // 調餐
    // TODO
    double weightHPWL = 1; 
    double weightOverlap = 0.5;


    // double bestCost = computeCost();
    // updateBest();

    int op1_ac = 0;

    // cost function
    // 用wirelength "變化量" + overlap "變化量"來看 
    double lastACrate = 1;
    double curTemperature = initialTemperature;
    while(curTemperature >= terminate_temperature)
    {
        //cout << "temp: " << curTemperature << " | decreaseRate: " << decreaseRate;

        double totalAttempt = 0, reject = 0, upHill = 0;

        while(upHill <= N && totalAttempt <= 2*N)
        {
            auto end_time_parser = chrono::high_resolution_clock::now();
            auto duration_parser = chrono::duration_cast<chrono::milliseconds>(end_time_parser-start_time).count();

            if(duration_parser > 590000)
            {
                return;
            }

            // neighborhood structure
            double perturbProb = (double)rand()/(RAND_MAX+1.0);
            if(perturbProb < 0.001) // TODO (0.8)
            {
                // option 1 : 把一個cell move到其他位置
                
                // 先挑選一個幸運movable cell
                int cellIdx = -1;
                do {cellIdx = rand() % _placement.numModules();} 
                while(_placement.module(cellIdx).isFixed());

                // 紀錄一下他的原位置
                double old_center_x = _placement.module(cellIdx).centerX();
                double old_center_y = _placement.module(cellIdx).centerY();

                // 再來幫他挑一個幸運位置
                double width = _placement.module(cellIdx).width();
                double height = _placement.module(cellIdx).height();
                double new_center_x = (rand() % static_cast<int>(coreWidth - width) + _placement.boundryLeft()) + width/2;
                double new_center_y = (rand() % static_cast<int>(coreHeight - height) + _placement.boundryBottom()) + height/2;
                
                // 先算module舊的WL和overlap量
                double moduleOldWL = 0, moduleOldOverlap = 0;

                // HPWL
                moduleOldWL = computeHPWLCost(cellIdx);

                // overlap
                unordered_set<int> uset_ovmodule;
                double old_leftX = _placement.module(cellIdx).x();
                double old_rightX = _placement.module(cellIdx).x() + width;

                double old_botY = _placement.module(cellIdx).y();
                double old_topY = _placement.module(cellIdx).y() + height;
                
                //tree_xDir.queryOverlaps(old_leftX, old_rightX, uset_ovmodule);
                //tree_yDir.queryOverlaps(old_botY, old_topY, uset_ovmodule);

                //moduleOldOverlap = computeModuleOverlap(cellIdx, uset_ovmodule);

                double curCost = weightHPWL*(moduleOldWL/maxHPWL) + weightOverlap*(moduleOldOverlap/totalModuleArea);

                // 把這個cell丟到新位置
                _placement.module(cellIdx).setCenterPosition(new_center_x, new_center_y);

                // 算module移動後的WL和overlap量
                double moduleNewWL = 0, moduleNewOverlap = 0;

                // HPWL
                moduleNewWL = computeHPWLCost(cellIdx);

                // overlap
                uset_ovmodule.clear();
                double new_leftX = _placement.module(cellIdx).x();
                double new_rightX = _placement.module(cellIdx).x() + width;

                double new_botY = _placement.module(cellIdx).y();
                double new_topY = _placement.module(cellIdx).y() + height;

                //tree_xDir.queryOverlaps(new_leftX, new_rightX, uset_ovmodule);
                //tree_yDir.queryOverlaps(new_botY, new_topY, uset_ovmodule);

                //moduleNewOverlap = computeModuleOverlap(cellIdx, uset_ovmodule);

                // 算更新後的總wirelength和overlap量
                currentWL = currentWL + (moduleNewWL - moduleOldWL);
                currentOverlap = currentOverlap + (moduleNewOverlap - moduleOldOverlap);
                
                // 接著算他的cost
                //double newCost = weightHPWL*((moduleNewWL - moduleOldWL)/maxHPWL) + weightOverlap*((moduleNewOverlap - moduleOldOverlap)/totalModuleArea);
                double newCost = weightHPWL*(moduleNewWL/maxHPWL) + weightOverlap*(moduleNewOverlap/totalModuleArea);

                if(newCost < curCost)   // down hill
                {
                    // modify 區間樹
                    //tree_xDir.modifyInterval(cellIdx, new_leftX, new_rightX);
                    //tree_yDir.modifyInterval(cellIdx, new_botY, new_topY);
                    op1_ac++;
                }
                else 
                {
                    double changeProbability = exp(-(newCost-curCost)/curTemperature);
                    double p = (double)rand()/(RAND_MAX+1.0);

                    //cout << "op1 p: " << p << endl;

                    if(p < changeProbability && curTemperature < greedyBound)   // accept     // TODO
                    {
                        // modify 區間樹
                      //  tree_xDir.modifyInterval(cellIdx, new_leftX, new_rightX);
                        //tree_yDir.modifyInterval(cellIdx, new_botY, new_topY);
                        op1_ac++;
                    }
                    else    // reject
                    {
                        _placement.module(cellIdx).setCenterPosition(old_center_x, old_center_y);                        
                        reject++;
                    }

                    upHill++;
                }
            }
            else
            {
                // option 2 : 交換兩個cell位置

                // 先挑選一個幸運movable cell
                int cellIdx_1 = -1;
                do {cellIdx_1 = rand() % _placement.numModules();} 
                while(_placement.module(cellIdx_1).isFixed());

                // 再挑一個不一樣的
                int cellIdx_2 = -1;
                do {cellIdx_2 = rand() % _placement.numModules();} 
                while(_placement.module(cellIdx_2).isFixed() || cellIdx_2 == cellIdx_1);

                // 紀錄一下他們的原位置
                double old_center_x_1 = _placement.module(cellIdx_1).centerX();
                double old_center_y_1 = _placement.module(cellIdx_1).centerY();

                double old_center_x_2 = _placement.module(cellIdx_2).centerX();
                double old_center_y_2 = _placement.module(cellIdx_2).centerY();

                double oldWL = 0, oldOverlap = 0;

                // 計算交換前的HPWL和overlap
                oldWL = computeHPWLCost(cellIdx_1) + computeHPWLCost(cellIdx_2);

                // overlap
                unordered_set<int> uset_ovmodule_1;
                double old_leftX_1 = _placement.module(cellIdx_1).x();
                double old_rightX_1 = _placement.module(cellIdx_1).x() + _placement.module(cellIdx_1).width();

                double old_botY_1 = _placement.module(cellIdx_1).y();
                double old_topY_1 = _placement.module(cellIdx_1).y() + _placement.module(cellIdx_1).height();

                //tree_xDir.queryOverlaps(old_leftX_1, old_rightX_1, uset_ovmodule_1);
                //tree_yDir.queryOverlaps(old_botY_1, old_topY_1, uset_ovmodule_1);

                unordered_set<int> uset_ovmodule_2;
                double old_leftX_2 = _placement.module(cellIdx_2).x();
                double old_rightX_2 = _placement.module(cellIdx_2).x() + _placement.module(cellIdx_2).width();

                double old_botY_2 = _placement.module(cellIdx_2).y();
                double old_topY_2 = _placement.module(cellIdx_2).y() + _placement.module(cellIdx_2).height();;

                //tree_xDir.queryOverlaps(old_leftX_2, old_rightX_2, uset_ovmodule_2);
                //tree_yDir.queryOverlaps(old_botY_2, old_topY_2, uset_ovmodule_2);

                //oldOverlap = computeModuleOverlap(cellIdx_1, uset_ovmodule_1) + computeModuleOverlap(cellIdx_2, uset_ovmodule_2) - computeTwoModuleOverlap(cellIdx_1, cellIdx_2);

                double curCost = weightHPWL*(oldWL/maxHPWL) + weightOverlap*(oldOverlap/totalModuleArea);

                // 交換他們
                _placement.module(cellIdx_1).setCenterPosition(old_center_x_2, old_center_y_2);
                _placement.module(cellIdx_2).setCenterPosition(old_center_x_1, old_center_y_1);

                // 計算交換後的HPWL和overlap

                double newWL = 0, newOverlap = 0;

                newWL = computeHPWLCost(cellIdx_1) + computeHPWLCost(cellIdx_2);

                // overlap
                uset_ovmodule_1.clear();
                double new_leftX_1 = _placement.module(cellIdx_1).x();
                double new_rightX_1 = _placement.module(cellIdx_1).x() + _placement.module(cellIdx_1).width();

                double new_botY_1 = _placement.module(cellIdx_1).y();
                double new_topY_1 = _placement.module(cellIdx_1).y() + _placement.module(cellIdx_1).height();

                //tree_xDir.queryOverlaps(new_leftX_1, new_rightX_1, uset_ovmodule_1);
                //tree_yDir.queryOverlaps(new_botY_1, new_topY_1, uset_ovmodule_1);

                uset_ovmodule_2.clear();
                double new_leftX_2 = _placement.module(cellIdx_2).x();
                double new_rightX_2 = _placement.module(cellIdx_2).x() + _placement.module(cellIdx_2).width();

                double new_botY_2 = _placement.module(cellIdx_2).y();
                double new_topY_2 = _placement.module(cellIdx_2).y() + _placement.module(cellIdx_2).height();;

               // tree_xDir.queryOverlaps(new_leftX_2, new_rightX_2, uset_ovmodule_2);
                //tree_yDir.queryOverlaps(new_botY_2, new_topY_2, uset_ovmodule_2);

                //newOverlap = computeModuleOverlap(cellIdx_1, uset_ovmodule_1) + computeModuleOverlap(cellIdx_2, uset_ovmodule_2) - computeTwoModuleOverlap(cellIdx_1, cellIdx_2);


                // 算更新後的總wirelength和overlap量
                currentWL = currentWL + (newWL - oldWL);
                currentOverlap = currentOverlap + (newOverlap - oldOverlap);

                // 接著算他的cost
                //double newCost = weightHPWL*(currentWL/initialWL) + weightOverlap*(currentOverlap/initialOverlap);
                //double newCost = weightHPWL*((newWL - oldWL)/maxHPWL) + weightOverlap*((newOverlap - oldOverlap)/totalModuleArea);

                double newCost = weightHPWL*(newWL/maxHPWL) + weightOverlap*(newOverlap/totalModuleArea);

                if(newCost <= curCost)   // down hill
                {

                    // // 看看best solution能不能update
                    // if(newCost < bestCost)
                    // {
                    //     bestCost = newCost;
                    //     updateBest();
                    // }
                    //tree_xDir.modifyInterval(cellIdx_1, new_leftX_1, new_rightX_1);
                    //tree_yDir.modifyInterval(cellIdx_1, new_botY_1, new_topY_1);

                    //tree_xDir.modifyInterval(cellIdx_2, new_leftX_2, new_rightX_2);
                    //tree_yDir.modifyInterval(cellIdx_2, new_botY_2, new_topY_2);
                }
                else 
                {
                    double changeProbability = exp(-(newCost-curCost)/curTemperature);
                    double p = (double)rand()/(RAND_MAX+1.0);

                    //cout << "op2 p: " << p << endl;

                    if(p < changeProbability && curTemperature < greedyBound)   // accept // TODO
                    {
                        //cout << "accept with bad solution" << endl;
                        //tree_xDir.modifyInterval(cellIdx_1, new_leftX_1, new_rightX_1);
                        //tree_yDir.modifyInterval(cellIdx_1, new_botY_1, new_topY_1);

                        //tree_xDir.modifyInterval(cellIdx_2, new_leftX_2, new_rightX_2);
                        //tree_yDir.modifyInterval(cellIdx_2, new_botY_2, new_topY_2);
                    }
                    else    // reject
                    {
                        _placement.module(cellIdx_1).setCenterPosition(old_center_x_1, old_center_y_1);           
                        _placement.module(cellIdx_2).setCenterPosition(old_center_x_2, old_center_y_2);                  
                        reject++;
                    }

                    upHill++;
                }
            }

            totalAttempt++;
        }

        // based on VPR dynamic change decreaseRate by acceptance rate
        double acceptanceRate = 1 - (reject/(totalAttempt));
        
        //cout << " | acceptRate: " << acceptanceRate << endl;

        curTemperature = curTemperature*decreaseRate;

        if(curTemperature < fixBound)
            decreaseRate = fixBoundRate;

        // if(acceptanceRate > 0.95)
        //     decreaseRate = 0.8;
        // else if(0.8 < acceptanceRate && acceptanceRate <= 0.95)
        //     decreaseRate = 0.95;  
        // else if(0.15 < acceptanceRate && acceptanceRate <= 0.8)
        //     decreaseRate = 0.995; 
        // else if(acceptanceRate <= 0.15)
        //     decreaseRate = 0.8; 
    }

    // // 最後，根據best result更新cell位置
    // for (size_t i = 0; i < _placement.numModules(); ++i)
    // {
    //     double center_x = _vec_movableCell_to_bestResult[i].first;
    //     double center_y = _vec_movableCell_to_bestResult[i].second;

    //    _placement.module(i).setCenterPosition(center_x, center_y);
    // }

    // set module coordinate
    for (size_t i = 0; i < _placement.numModules(); i++)
    {
        if(_placement.module(i).isFixed())
            continue;

        double new_center_x = _placement.module(i).centerX();
        double new_center_y = _placement.module(i).centerY();

        double module_width = _placement.module(i).width();
        double module_height = _placement.module(i).height();

        if(new_center_x - (module_width/2) < _placement.boundryLeft())
        {
            //cout << "refine!" << endl;
            new_center_x = _placement.boundryLeft() + (module_width/2) + 1;
        }

        if(new_center_x + (module_width/2) > _placement.boundryRight())
        {
            //cout << "refine!" << endl;
            new_center_x = _placement.boundryRight() - (module_width/2) - 1;            
        }

        if(new_center_y - (module_height/2) < _placement.boundryBottom())
        {
            //cout << "refine!" << endl;
            new_center_y = _placement.boundryBottom() + (module_height/2) + 1;
        }

        if(new_center_y + (module_height/2) > _placement.boundryTop())
        {
            //cout << "refine!" << endl;
            new_center_y = _placement.boundryTop() - (module_height/2) - 1;            
        }

        //cout << new_center_x << " " << new_center_y << endl;
        _placement.module(i).setCenterPosition(new_center_x, new_center_y);
    }

    cout << "use op1 aceept cnt: " << op1_ac << endl;
    return;
}

//// #############################################################################################

void GlobalPlacer::place2()
{
    ///////////////////////////////////////////////////////////////////
    // The following example is only for analytical methods.
    // if you use other methods, you can skip and delete it directly.
    //////////////////////////////////////////////////////////////////

    cout << "cell information: " << "t" << _placement.boundryTop() << " b" << _placement.boundryBottom() << " l" << _placement.boundryLeft() << " r" << _placement.boundryRight() << endl; 
    cout << "module_num" << _placement.numModules() << " " << "net_num" << _placement.numNets() <<  " " << "pin_num" << _placement.numPins() << endl;

    // 把整個boundary切成m*m個grid(bin)
    int m = ceil(log2(sqrt(_placement.numModules())));
    int grid_dimension = min(m, 1024);

    int grid_dimension_x = 4, grid_dimension_y = 4;
    cout << "m= " << grid_dimension << endl;
    cout << "n= " << grid_dimension_y << endl;

    double grid_width = (_placement.boundryRight() - _placement.boundryLeft())/grid_dimension_x;
    double grid_height = (_placement.boundryTop() - _placement.boundryBottom())/grid_dimension_y;

    cout << grid_width << "w h" << grid_height << endl;

    // start cutting bin
    vector<pair<double, double>> vec_tempV;
    for(double bottomleft_x = _placement.boundryLeft(); bottomleft_x + grid_width <= _placement.boundryRight(); bottomleft_x += grid_width)
    {
        for(double bottomleft_y = _placement.boundryBottom(); bottomleft_y + grid_height <= _placement.boundryTop(); bottomleft_y += grid_height)
        {
            double gird_mid_x = bottomleft_x + (grid_width)/2, gird_mid_y = bottomleft_y + (grid_height)/2;

            vec_tempV.push_back({gird_mid_x, gird_mid_y});
        }
    }

    cout << "after cut:" << vec_tempV.size() << endl;

    // 計算target density
   
    double chipArea = (_placement.boundryRight() - _placement.boundryLeft())*(_placement.boundryTop() - _placement.boundryBottom());
    double cellTotalArea = 0;
    for(size_t i = 0; i < _placement.numModules(); i++)
        cellTotalArea += _placement.module(i).area();
    double targetD = cellTotalArea/(vec_tempV.size());

    cout << "target density: " << targetD << endl;

    // An example of random placement implemented by TA.
    // If you want to use it, please uncomment the folllwing 1 line.
    //randomPlace();

    vector<pair<double, double>> vec_old, vec_new;
    for (size_t i = 0; i < _placement.numModules(); ++i)
    {
        if(!_placement.module(i).isFixed())
            vec_old.push_back({_placement.module(i).centerX(), _placement.module(i).centerY()});
    }


    double init_res = _placement.computeHpwl();

    printf("initial HPWL: %.0f\n", init_res);

    int maxIteration = 5;
    for(int iter = 0; iter < maxIteration; iter++)
    {
        // _placement.numModules()
        for(int curModuleIdx = 0; curModuleIdx < 100; curModuleIdx++)
        {
            cout << "cc:" << curModuleIdx << endl;
            if(_placement.module(curModuleIdx).isFixed())
                continue;

            ExampleFunction ef(_placement); // require to define the object function and gradient function
            ef._smooth_parameter = 0.01;    // user defined

            ef._wb = grid_width;
            ef._hb = grid_height;

            ef._num_blocks = 1;
            ef._cur_module = curModuleIdx;

            ef._vec_pin_x_diff = vector<pair<int, double>>(_placement.numPins());
            ef._vec_pin_y_diff = vector<pair<int, double>>(_placement.numPins()); 

            ef._vec_bin_mid_pt_cord = vector<std::pair<double, double>>(vec_tempV.begin(), vec_tempV.end());
            ef._target_density = targetD;

            vector<double> vec_module_center(2);
            vec_module_center[0] = _placement.module(curModuleIdx).centerX();
            vec_module_center[1] = _placement.module(curModuleIdx).centerY();

            for(int pin_idx = 0; pin_idx < _placement.numPins(); pin_idx++)
            {
                double pin_cord_x = _placement.pin(pin_idx).x();
                double pin_cord_y = _placement.pin(pin_idx).y();

                int connect_net_idx = _placement.pin(pin_idx).netId();
                int connect_module_idx = _placement.pin(pin_idx).moduleId();

                double module_center_x = _placement.module(connect_module_idx).centerX();
                double module_center_y = _placement.module(connect_module_idx).centerY();

                // module固定位置，代表pin的位置不會再改變了
                if (connect_module_idx != curModuleIdx)
                {
                    ef._vec_pin_x_diff[pin_idx] = {-1, pin_cord_x};
                    ef._vec_pin_y_diff[pin_idx] = {-1, pin_cord_y};
                }
                else
                {
                    ef._vec_pin_x_diff[pin_idx] = {ef._map_originalIdx_to_efIdx[0], pin_cord_x-module_center_x};
                    ef._vec_pin_y_diff[pin_idx] = {ef._map_originalIdx_to_efIdx[0], pin_cord_y-module_center_y};
                }

                ef._map_netIdx_to_vecPinIdx[connect_net_idx].push_back(pin_idx);
            }
            
            NumericalOptimizer no(ef);
            no.setX(vec_module_center);  // set initial solution
            no.setNumIteration(25); // user-specified parameter
            no.setStepSizeBound(10); // user-specified parameter
            no.solve();             // Conjugate Gradient solver

            
            // set module coordinate
            double new_center_x = no.x(0);
            double new_center_y = no.x(1);

            double module_width = _placement.module(curModuleIdx).width();
            double module_height = _placement.module(curModuleIdx).height();

            if(new_center_x - (module_width/2) < _placement.boundryLeft())
            {
                cout << "refine!" << endl;
                new_center_x = _placement.boundryLeft() + (module_width/2) + 1;
            }

            if(new_center_x + (module_width/2) > _placement.boundryRight())
            {
                cout << "refine!" << endl;
                new_center_x = _placement.boundryRight() - (module_width/2) - 1;            
            }

            if(new_center_y - (module_height/2) < _placement.boundryBottom())
            {
                cout << "refine!" << endl;
                new_center_y = _placement.boundryBottom() + (module_height/2) + 1;
            }

            if(new_center_y + (module_height/2) > _placement.boundryTop())
            {
                cout << "refine!" << endl;
                new_center_y = _placement.boundryTop() - (module_height/2) - 1;            
            }

            vec_new.push_back({new_center_x, new_center_y});
            _placement.module(curModuleIdx).setCenterPosition(new_center_x, new_center_y);
        }

        for(int i = 0; i < vec_old.size(); i++)
        {
            cout << "old cord: " << vec_old[i].first << " " << vec_old[i].second;
            cout << " | new cord: " << vec_new[i].first << " " << vec_new[i].second << endl;
        }
    }



    // double tt;
    // ef.evaluateF(vec_module_center, tt);

	// cout << "Current solution:\n";
	// for (unsigned i = 0; i < no.dimension(); i++)
	// {
	// 	cout << "x[" << i << "] = " << no.x(i) << "\n";
	// }
	// cout << "Objective: " << no.objective() << "\n";
    


    printf("initial HPWL: %.0f\n", init_res);

    /* @@@ TODO
     * 1. Understand above example and modify ExampleFunction.cpp to implement the analytical placement
     * 2. You can choose LSE or WA as the wirelength model, the former is easier to calculate the gradient
     * 3. For the bin density model, you could refer to the lecture notes
     * 4. You should first calculate the form of wirelength model and bin density model and the forms of their gradients ON YOUR OWN
     * 5. Replace the value of f in evaluateF() by the form like "f = alpha*WL() + beta*BinDensity()"
     * 6. Replace the form of g[] in evaluateG() by the form like "g = grad(WL()) + grad(BinDensity())"
     * 7. Set the initial vector x in place(), set step size, set #iteration, and call the solver like above example
     * */
}

//// #############################################################################################

void GlobalPlacer::place()
{
    ///////////////////////////////////////////////////////////////////
    // The following example is only for analytical methods.
    // if you use other methods, you can skip and delete it directly.
    //////////////////////////////////////////////////////////////////

    ExampleFunction ef(_placement); // require to define the object function and gradient function
    ef._smooth_parameter = 0.001;    // user defined

    for(size_t i = 0; i < _placement.numModules(); i++)
    {
        if (_placement.module(i).isFixed())
            continue;

        ef._set_movable_cell_idx.insert(i);
    }
    ef._num_blocks = ef._set_movable_cell_idx.size();

    ef._vec_pin_x_diff = vector<pair<int, double>>(_placement.numPins());
    ef._vec_pin_y_diff = vector<pair<int, double>>(_placement.numPins());    

    ////////////////////////////////////////////////////////////////

    cout << "cell information: " << "t" << _placement.boundryTop() << " b" << _placement.boundryBottom() << " l" << _placement.boundryLeft() << " r" << _placement.boundryRight() << endl; 
    cout << "module_num" << _placement.numModules() << " " << "net_num" << _placement.numNets() <<  " " << "pin_num" << _placement.numPins() << " " << "movable_num" << ef._num_blocks << endl;

    // 把整個boundary切成m*m個grid(bin)
    int m = ceil(log2(sqrt(_placement.numModules())));
    int grid_dimension = min(m, 1024);

    int grid_dimension_x = 14, grid_dimension_y = 14;
    cout << "m= " << grid_dimension << endl;
    cout << "n= " << grid_dimension_y << endl;

    double grid_width = (_placement.boundryRight() - _placement.boundryLeft())/grid_dimension_x;
    double grid_height = (_placement.boundryTop() - _placement.boundryBottom())/grid_dimension_y;

    ef._wb = grid_width;
    ef._hb = grid_height;

    cout << grid_width << "w h" << grid_height << endl;

    // start cutting bin
    for(double bottomleft_x = _placement.boundryLeft(); bottomleft_x + grid_width <= _placement.boundryRight(); bottomleft_x += grid_width)
    {
        for(double bottomleft_y = _placement.boundryBottom(); bottomleft_y + grid_height <= _placement.boundryTop(); bottomleft_y += grid_height)
        {
            double gird_mid_x = bottomleft_x + (grid_width)/2, gird_mid_y = bottomleft_y + (grid_height)/2;

            ef._vec_bin_mid_pt_cord.push_back({gird_mid_x, gird_mid_y});
        }
    }

    cout << "after cut:" << ef._vec_bin_mid_pt_cord.size() << endl;

    // 計算target density
    double chipArea = (_placement.boundryRight() - _placement.boundryLeft())*(_placement.boundryTop() - _placement.boundryBottom());
    double cellTotalArea = 0;
    for(size_t i = 0; i < _placement.numModules(); i++)
        cellTotalArea += _placement.module(i).area();
    ef._target_density = cellTotalArea / chipArea;

    cout << "target density: " << ef._target_density << endl;

    // An example of random placement implemented by TA.
    // If you want to use it, please uncomment the folllwing 1 line.
    randomPlace();

    double init_res = _placement.computeHpwl();

    printf("initial HPWL: %.0f\n", init_res);

    // 為了求HPWL，需要算出每個pin的x-y座標且用module x-y center座標function表示
    // step1. 求出pin的x-y座標function
    int ef_idx = 0;
    for(int module_idx = 0; module_idx < _placement.numModules(); module_idx++)
    {
        double cur_module_center_x = _placement.module(module_idx).centerX();
        double cur_module_center_y = _placement.module(module_idx).centerY();

        //cout << module_idx << " " << cur_module_center_x << " " << cur_module_center_y << endl;
        if (_placement.module(module_idx).isFixed())
            continue;

        ef._map_originalIdx_to_efIdx[module_idx] = ef_idx++; 
    }

    for(int pin_idx = 0; pin_idx < _placement.numPins(); pin_idx++)
    {
        double pin_cord_x = _placement.pin(pin_idx).x();
        double pin_cord_y = _placement.pin(pin_idx).y();

        int connect_net_idx = _placement.pin(pin_idx).netId();
        int connect_module_idx = _placement.pin(pin_idx).moduleId();

        double module_center_x = _placement.module(connect_module_idx).centerX();
        double module_center_y = _placement.module(connect_module_idx).centerY();

        // module固定位置，代表pin的位置不會再改變了
        if (_placement.module(connect_module_idx).isFixed())
        {
            ef._vec_pin_x_diff[pin_idx] = {-1, pin_cord_x};
            ef._vec_pin_y_diff[pin_idx] = {-1, pin_cord_y};
        }
        else
        {
            ef._vec_pin_x_diff[pin_idx] = {ef._map_originalIdx_to_efIdx[connect_module_idx], pin_cord_x-module_center_x};
            ef._vec_pin_y_diff[pin_idx] = {ef._map_originalIdx_to_efIdx[connect_module_idx], pin_cord_y-module_center_y};
        }

        ef._map_netIdx_to_vecPinIdx[connect_net_idx].push_back(pin_idx);
    }

    int maxIter = 5;
    for(int iter = 0; iter < maxIter; iter++)
    {
        vector<double> vec_module_center;

        for(int module_idx = 0; module_idx < _placement.numModules(); module_idx++)
        {
            double cur_module_center_x = _placement.module(module_idx).centerX();
            double cur_module_center_y = _placement.module(module_idx).centerY();

            //cout << module_idx << " " << cur_module_center_x << " " << cur_module_center_y << endl;
            if (_placement.module(module_idx).isFixed())
                continue;

            vec_module_center.push_back(cur_module_center_x);
            vec_module_center.push_back(cur_module_center_y);
        }

        ef._alpha = 1;
        ef._beta = iter*6000;

        NumericalOptimizer no(ef);
        no.setX(vec_module_center);  // set initial solution
        no.setNumIteration(30); // user-specified parameter
        no.setStepSizeBound((_placement.boundryRight() - _placement.boundryLeft())*7); // user-specified parameter
        no.solve();             // Conjugate Gradient solver

        // set module coordinate
        for (size_t i = 0; i < _placement.numModules(); i++)
        {
            if (_placement.module(i).isFixed())
                continue;

            int ef_idx = ef._map_originalIdx_to_efIdx[i];

            double new_center_x = no.x(ef_idx*2);
            double new_center_y = no.x(ef_idx*2 + 1);

            double module_width = _placement.module(i).width();
            double module_height = _placement.module(i).height();

            if(new_center_x - (module_width/2) < _placement.boundryLeft())
            {
                cout << "refine!" << endl;
                new_center_x = _placement.boundryLeft() + (module_width/2) + 1;
            }

            if(new_center_x + (module_width/2) > _placement.boundryRight())
            {
                cout << "refine!" << endl;
                new_center_x = _placement.boundryRight() - (module_width/2) - 1;            
            }

            if(new_center_y - (module_height/2) < _placement.boundryBottom())
            {
                cout << "refine!" << endl;
                new_center_y = _placement.boundryBottom() + (module_height/2) + 1;
            }

            if(new_center_y + (module_height/2) > _placement.boundryTop())
            {
                cout << "refine!" << endl;
                new_center_y = _placement.boundryTop() - (module_height/2) - 1;            
            }

            cout << new_center_x << " " << new_center_y << endl;
            _placement.module(i).setCenterPosition(new_center_x, new_center_y);
        }
    }


    // double tt;
    // ef.evaluateF(vec_module_center, tt);

	// cout << "Current solution:\n";
	// for (unsigned i = 0; i < no.dimension(); i++)
	// {
	// 	cout << "x[" << i << "] = " << no.x(i) << "\n";
	// }
	// cout << "Objective: " << no.objective() << "\n";


    printf("initial HPWL: %.0f\n", init_res);

    /* @@@ TODO
     * 1. Understand above example and modify ExampleFunction.cpp to implement the analytical placement
     * 2. You can choose LSE or WA as the wirelength model, the former is easier to calculate the gradient
     * 3. For the bin density model, you could refer to the lecture notes
     * 4. You should first calculate the form of wirelength model and bin density model and the forms of their gradients ON YOUR OWN
     * 5. Replace the value of f in evaluateF() by the form like "f = alpha*WL() + beta*BinDensity()"
     * 6. Replace the form of g[] in evaluateG() by the form like "g = grad(WL()) + grad(BinDensity())"
     * 7. Set the initial vector x in place(), set step size, set #iteration, and call the solver like above example
     * */
}

