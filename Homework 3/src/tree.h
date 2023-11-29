#ifndef TREE_H_DEF
#define TREE_H_DEF

#include "module.h"

struct Node
{
    int moduleIdx;
    Node* left, *right;
};

#endif