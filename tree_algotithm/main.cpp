#include <iostream>
#include "tree.h"
#include <vector>

using namespace std;

int main()
{
    double eps = 0.001;
/*
 * following mesh was used
*/
    double i = 0.2;
    double j = 0.4;
    double k = 0.8;
    int a = 5000;
    int b = 10000;
    int c = 20000;
    int d = 40000;

    /*
    tree* mytree = new tree(5*(a), a, i, eps);
    mytree->simulation();
    delete mytree;
    tree* mytree2 = new tree(5*(b), b, i, eps);
    mytree2->simulation();
    delete mytree2;
    tree* mytree3 = new tree(5*(c), c, i, eps);
    mytree3->simulation();
    delete mytree3;
    tree* mytree4 = new tree(5*(d), d, i, eps);
    mytree4->simulation();
    delete mytree4;

    tree* mytree5 = new tree(5*(a), a, j, eps);
    mytree5->simulation();
    delete mytree5;
    tree* mytree6 = new tree(5*(b), b, j, eps);
    mytree6->simulation();
    delete mytree6;
    tree* mytree7 = new tree(5*(c), c, j, eps);
    mytree7->simulation();
    delete mytree7;
    tree* mytree8 = new tree(5*(d), d, j, eps);
    mytree8->simulation();
    delete mytree8;
*/
    tree* mytree9 = new tree(5*(a), a, k, eps);
    mytree9->simulation();
    delete mytree9;
    tree* mytree10 = new tree(5*(b), b, k, eps);
    mytree10->simulation();
    delete mytree10;
    tree* mytree11 = new tree(5*(c), c, k, eps);
    mytree11->simulation();
    delete mytree11;
    tree* mytree12 = new tree(5*(d), d, k, eps);
    mytree12->simulation();
    delete mytree12;
    return 0;
}
