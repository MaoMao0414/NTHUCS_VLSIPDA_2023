#include "module.h"

bool cmpIncreaseXCord(Module* lhs, Module* rhs)
{
    return lhs->cordX < rhs->cordX;
}

bool cmpDecreaseArea(Module* lhs, Module* rhs)
{
    return lhs->minArea > rhs->minArea;
}

