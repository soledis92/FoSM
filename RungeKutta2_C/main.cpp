#include <iostream>
#include <cmath>
#include <vector>
#include "rk2.hpp"
#include <stdio.h>
#include <fstream>
using namespace std;


//Constants

const double m1 = 0.5;
const double m2 = 1.0;
const double l1 = 2.0;
const double l2 = 1.0;

const double g = 1;

const double dt = 0.05;
const double T = 100.0;

// define the differential equation which is has to be solved

void dgl(double* q, double t, double* dq_dt, int dim){

    //first coordinate:  dphi1_dt
    //dq_dt[0] = (q[2] / (pow(l1,2) * (m1 + m2)) - (q[3] * cos(q[0] - q[1])) / (l1 * l2 * (m1 +m2)))  *
    //        (1. / (1. - (m1 * pow(cos(q[0] - q[1]), 2) / (m1 + m2) )));
    dq_dt[0] = (q[2] - (l1/l2) * cos(q[0]-q[1]) * q[3])  /  (pow(l1, 2) * (m1+m2) * (1 - pow(cos(q[0]-q[1]),2)));

    // second coordinate: dphi2_dt
    dq_dt[1] = ((q[2]) - m2*l1*l2* dq_dt[0] * cos(q[0] - q[1])) / (m2 * pow(l2,2));



    // third coordinate: dq1_dt  (NICHT ZU VERWECHSELN MIT DEM VEKTOR dq_dt HIER !!!!! dq1_dt = dq_dt[2] )
    dq_dt[2] = - m2*l1*l2* dq_dt[0] * dq_dt[1] * sin(q[0] - q[1]) - (m1+m2)*g*l1 * sin(q[0]);



    // fourth coordinate: dq2_dt (NICHT ZU VERWECHSELN MIT DEM VEKTOR dq_dt HIER !!!!! dq1_dt = dq_dt[3] )
    dq_dt[3] = m2*l1*l2 * dq_dt[0] * dq_dt[1] * sin(q[0] - q[1]) - m2*g*l2 * sin(q[1]);


}



int main(){

    // Initial values

    double phi1 = 50.0 / 180.0 * M_PI;
    double phi2 = -120 / 180.0 * M_PI;
    double dphi1_dt = 0.0;
    double dphi2_dt = 0.0;
    double q1 = 0.0;
    double q2 = 0.0;
    int n = 4;

    // Initial value vector y0

    double* y0 = (double*)malloc(n*sizeof(double));
    y0[0] = phi1;
    y0[1] = phi2;
    y0[2] = q1;
    y0[3] = q2;

    double t = 0;

    fstream file;
    file.open("FoSM2_3b.dat", ios::out);


    while(t < T) {
        file << t <<" "<< y0[0] <<" "<< y0[1]<<" "<<y0[2]<<" "<<y0[3]<<endl;
        rk2(dgl, y0, t, dt, n);
        t += dt;
    }

    file.close();

    return 0;
}
