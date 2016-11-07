#ifndef PARTICLE_H
#define PARTICLE_H

class particle
{
public:
    particle();
    double pos[3];
    double mass;
    double acc_tree[3]; //acceleration of particle calculated from tree
    double acc_exact[3]; //acceleration of particle calculated by exact summation
};

#endif // PARTICLE_H
