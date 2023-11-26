#ifndef WIRE_H_DEF
#define WIRE_H_DEF

#include "cell.h"

struct Net
{
    std::string name;
    std::vector<Cell*> vec_connectCell;
};

#endif