#ifndef TREE_NODE_H
#define TREE_NODE_H
#include "particle.h"

class tree_node
{
public:
    tree_node();
    double center[3];    // geom. center
    double len;

    double cm[3];         // center of mass
    double mass;

    tree_node *suns[8];    //subnodes
    particle* p;

    bool empty;
};

#endif // TREE_NODE_H
