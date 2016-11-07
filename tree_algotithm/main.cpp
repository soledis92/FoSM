#include <iostream>
#include "tree.h"

using namespace std;

int main()
{
    int max_points = 5000;
    int max_nodes = 5* max_points;
    double opening_threshold = 0.8;
    double eps = 0.001;
    tree* mytree = new tree(max_nodes, max_points, opening_threshold, eps);
    mytree->simulation();
    delete mytree;
    return 0;
}
