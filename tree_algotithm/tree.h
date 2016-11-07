#ifndef TREE_H
#define TREE_H

#include "tree_node.h"
#include "particle.h"
#include <vector>

class tree
{
public:
    tree(int max_nodes, int max_points, double opening_threshold, double eps);
    void simulation();
private:
    int max_nodes = 5000;
    int max_points = 5* max_nodes;
    double opening_threshold;
    double eps;
    //store data in array:
    std::vector <tree_node> data_tree;
    std::vector <particle> star;
    // variables:
    int count_nodes = 0;
    //methods:
    tree_node* get_empty_node();
    int get_subnode_index(tree_node* current, particle* p);
    void insert_particle(tree_node* current, particle* pnew);
    void calc_multipole_moments(tree_node* current);
    double get_opening_angle(tree_node* current, double pos[3]);
    void walk_tree(tree_node* current, double pos[3], double acc[3]);
};

#endif // TREE_H
