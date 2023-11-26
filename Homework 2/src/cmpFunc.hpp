#include "cell.h"

// 1, 2, 3, 4, 5
static bool comparefun1(Cell* lhs, Cell* rhs)
{
    return lhs->arr_cellArea[0]-lhs->arr_cellArea[1] > rhs->arr_cellArea[0]-rhs->arr_cellArea[1];     // dieB面積越大，dieA面積越小的排前面
}

// 1, 2, 3, 4, 5
static bool comparefun2(Cell* lhs, Cell* rhs)
{
    return lhs->arr_cellArea[0]-lhs->arr_cellArea[1] < rhs->arr_cellArea[0]-rhs->arr_cellArea[1];     // dieB面積越大，dieA面積越小的排後面
}

static bool comparefun3(Cell* lhs, Cell* rhs)
{
    return lhs->cellAratio < rhs->cellAratio;     // 用die A佔的ratio由小到大排
}

static bool comparefun4(Cell* lhs, Cell* rhs)
{
    return lhs->cellAratio > rhs->cellAratio;     // 用die A佔的ratio由大到小排
}

// 1, 2, 3, 5
static bool comparefun5(Cell* lhs, Cell* rhs)
{
    return lhs->arr_cellArea[0] < rhs->arr_cellArea[0];     // 用die A的cell由小到大排
}

// 1, 2, 3
static bool comparefun6(Cell* lhs, Cell* rhs)
{
    return lhs->arr_cellArea[0] > rhs->arr_cellArea[0];     // 用die A的cell由大到小排
}

// 1, 2, 3, 5
static bool comparefun7(Cell* lhs, Cell* rhs)
{
    return lhs->arr_cellArea[1] < rhs->arr_cellArea[1];     // 用die B的cell由小到大排
}

// 1, 2, 3
static bool comparefun8(Cell* lhs, Cell* rhs)
{
    return lhs->arr_cellArea[1] > rhs->arr_cellArea[1];     // 用die B的cell由大到小排
}