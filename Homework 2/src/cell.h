#ifndef CELL_H_DEF
#define CELL_H_DEF

struct Net;

struct Cell
{
    std::string name;
    long long arr_cellArea[2];     // 0是放A的area, 1是放B的area
    std::vector<Net*> vec_connectNet;
    bool isLock;
    bool curDie;     // 位在哪個die : -1表未分配，0為dieA, 1為dieB
    long long gain;
    float cellAratio;
};

struct Net
{
    std::string name;
    std::vector<Cell*> vec_connectCell;
    long long arr_netDistribution[2];     // 紀錄位於A, B兩個die上的cell數量
};

#endif