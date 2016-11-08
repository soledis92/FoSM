#include "tree.h"
#include <iostream>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>


tree::tree(int max_nodes, int max_points, double opening_threshold, double eps)
{
    this->max_nodes = max_nodes;
    this->max_points = max_points;
    this->eps = eps;
    this->opening_threshold = opening_threshold;
    this->data_tree.resize(max_nodes);
    this->star.resize(max_points);
}


tree_node* tree::get_empty_node()
{
    tree_node *no;
    static int count_nodes = 0;

    if(count_nodes < this->max_nodes)
      {
        no = &data_tree[count_nodes++];        //first free entry of tree
        memset(no, 0, sizeof(tree_node));      //set node "no" to 0
      }
    else
      {
        std::cout << "sorry, we are out of tree nodes." << std::endl;
        exit(1);
      }

    return no;
}

int tree::get_subnode_index(tree_node *current, particle *p)
{
    int index = 0;

    if(p->pos[0] > current->center[0])
    {
        index += 4;
    }
    if(p->pos[1] > current->center[1])
    {
        index += 2;
    }
    if(p->pos[2] > current->center[2])
    {
        index += 1;
    }
    return index;
}

void tree::insert_particle(tree_node *current, particle *pnew)
{
    tree_node *sun;
    int i, j, k, n, p_subnode, pnew_subnode;

    if(current->p)
    {
        /* The node contains a particle already.
         Need to create a new set of 8 subnodes, and then move this particle to one of them */

        for(i = 0, n = 0; i<2; i++)
            for(j=0; j<2; j++)
                for(k=0; k<2; k++)
                {
                    sun = get_empty_node();
                    current->suns[n++] = sun;
                    sun->len = 0.5 * current->len;
                    sun->center[0] = current->center[0] + 0.25 * (2*i-1) * current->len;
                    sun->center[1] = current->center[1] + 0.25 * (2*j-1) * current->len;
                    sun->center[2] = current->center[2] + 0.25 * (2*k-1) * current->len;
                }

        /* determine in which subnode the old particle sits */
        p_subnode = get_subnode_index(current, current->p);

        /* now move the particle to this subnode */
        current->suns[p_subnode]->p = current->p;
        current->p = NULL;

        /* determine in which subnode the new particle sits */
        pnew_subnode = get_subnode_index(current, pnew);

        /* now try to insert the new particle there */
        insert_particle(current->suns[pnew_subnode], pnew);
    }
    else
    {
        /* check in which subnode the new particle would fall */
        pnew_subnode = get_subnode_index(current, pnew);

        /* if the corresponding subnode exists, we try to insert the particle there,
         otherwise we know there are no subnodes in the node, so we can put the particle into the current node */
        if(current->suns[pnew_subnode])
            insert_particle(current->suns[pnew_subnode], pnew);
        else
            current->p = pnew;
    }
}

void tree::calc_multipole_moments(tree_node *current)
{
    int n, j;

    if(current->suns[0])   /* do we have subnodes? */
    {
        /* yes, so let's first calculate their multipole moments */
        for(n = 0; n < 8; n++)
            calc_multipole_moments(current->suns[n]);

        /* initialize the node multipole moments to zero */
        current->mass  = 0;
        for(j = 0; j < 3; j++)
            current->cm[j] = 0;

        /* now calculate the moment of the current node from those of its suns */
        // calculate monopole: Q = sum(m_i)
        for (n=0; n < 8; n++)
        {
            current->mass += current->suns[n]->mass;
        }
        // calculate center of mass:
        for (n=0; n < 8; n++)
        {
            for(j=0; j<3; j++)
            {
                current->cm[j] += current->suns[n]->mass * current->suns[n]->cm[j] / current->mass;
            }
        }
    }
    else
    {
        if(current->p)  /* do we at least have a particle? */
        {
            /* yes, so let's copy this particle to the multipole moments of the node */

            current->mass = current->p->mass;
            for(j = 0; j < 3; j++)
                current->cm[j] = current->p->pos[j];
        }
        else
        {
            /* nothing in here at all; let's initialize the multipole moments to zero */
            current->mass  = 0;
            for(j = 0; j < 3; j++)
                current->cm[j] = 0;
        }
    }
}

double tree::get_opening_angle(tree_node *current, double pos[])
{
    int j;
    double r2 = 0;

    for(j = 0; j < 3; j++)
    {
        r2 += (current->cm[j] - pos[j]) * (current->cm[j] - pos[j]);
    }

    return current->len / (sqrt(r2) + 1.0e-35);       //"+ 1.0e-35" avoids divison by 0
}

void tree::walk_tree(tree_node *current, double pos[], double acc[])
{
    int n;
    double theta;

    if(current->mass)   /* only do something if there is mass in this branch of the tree (i.e. if it is not empty) */
    {
        theta = get_opening_angle(current, pos);

        /* if the node is seen under a small enough angle or contains a single particle,
           * we take its multipole expansion, and we're done for this branch
           */
        if(theta < opening_threshold || current->p)
        {
            calc_multipole_moments(current);
            acc[0] += current->mass / sqrt(current->cm[0]*current->cm[0] + current->cm[1]*current->cm[1] + current->cm[2]*current->cm[2]);
            acc[1] += current->mass / sqrt(current->cm[0]*current->cm[0] + current->cm[1]*current->cm[1] + current->cm[2]*current->cm[2]);
            acc[2] += current->mass / sqrt(current->cm[0]*current->cm[0] + current->cm[1]*current->cm[1] + current->cm[2]*current->cm[2]);
        }
        else
        {
            /* otherwise we have to open the node and look at all daughter nodes in turn */

            if(current->suns[0])             /* make sure that we actually have subnodes */
                for(n=0; n<8; n++)
                    walk_tree(current->suns[n], pos, acc);
        }
    }
}

void tree::simulation()
{
    tree_node *root;
    int i, j, N;
    double t0, t1;

    N = this->max_points;

    srand48(42);   /* set a random number seed */

    /* create a random particle set, uniformly distributed in a box */
    for(i=0; i < N; i++)
    {
        star.at(i).mass = 1.0 / N;       // total mass = 1

        for(j=0; j<3; j++)
            star.at(i).pos[j] = drand48();
    }

    /* create an empty root node for the tree */
    root = get_empty_node();


    /* set the dimension and position of the root node */
    root->len = 1.0;
    for(j=0; j<3; j++)
        root->center[j] = 0.5;

    /* insert the particles into the tree */
    for(i=0; i < N; i++)
        insert_particle(root, &star[i]);

    /* calculate the multipole moments */
    calc_multipole_moments(root);


    /* set a timer */
    t0 = (double) clock();

    /* now calculate the accelerations with the tree */
    for(i = 0; i < N; i++)
    {
        for(j = 0; j < 3; j++)
            star[i].acc_tree[j] = 0;

        walk_tree(root, star[i].pos, star[i].acc_tree);
    }

    t1 = (double) clock();
    std::cout << "force calculation with tree took: "
              << '\n'
              << (t1 - t0) / CLOCKS_PER_SEC << " sec" << std::endl;


    t0 = (double) clock();

    /* now calculate the accelerations with direct summation, for comparison */
    for(i = 0; i < N; i++)
    {
        for(j=0; j < N; j++)
        {
            if(i==j){continue;}
            star[i].acc_exact[0] += star[j].mass / sqrt(star[j].pos[0]*star[j].pos[0] + star[j].pos[1]*star[j].pos[1] + star[j].pos[2]*star[j].pos[2]);
            star[i].acc_exact[1] += star[j].mass / sqrt(star[j].pos[0]*star[j].pos[0] + star[j].pos[1]*star[j].pos[1] + star[j].pos[2]*star[j].pos[2]);
            star[i].acc_exact[2] += star[j].mass / sqrt(star[j].pos[0]*star[j].pos[0] + star[j].pos[1]*star[j].pos[1] + star[j].pos[2]*star[j].pos[2]);

        }
    }

    t1 = (double) clock();
    std::cout << "calculation with direct summation took: "
              << '\n'
              << (t1 - t0) / CLOCKS_PER_SEC << " sec" << std::endl;

    /* now do the calculation of the mean relative error
     */

    double err_sum = 0;
    double err[3] = {0};
    for(i = 0; i < N; i++)
    {
        for (int i = 0; i< 3; i++)
        {
            err[i] += fabs(star[i].acc_tree[i]-star[i].acc_exact[i]) / star[i].acc_exact[i];
        }
    }
    for(i = 0; i< 3; i++)
    {
        err[i] /= N;
        err_sum += err[i];
    }

    std::cout << " following parameters: " << "N = " << this->max_nodes << ", theta = " << this->opening_threshold << std::endl;
    std::cout << "Average relative error: "
              << err_sum
              << '\n'
              << "err_vec: a(x, y, z) = "
              << err[0] << ", " << err[1] << ", " << err[2] << '\n' << std::endl;
    /*
    for(i=0; i<3;i++)
    std::cout << "exact_a_"<< i << ": " << star[1].acc_exact[i] << " vs. tree_a_" << i << ": " << star[1].acc_tree[i] << std::endl;
    */
}
