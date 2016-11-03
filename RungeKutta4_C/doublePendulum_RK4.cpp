#include <iostream>
#include <cmath>
#include <vector>
#include "rk4.hpp"
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
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
     //       (1. / (1. - (m1 * pow(cos(q[0] - q[1]), 2) / (m1 + m2) )));
    dq_dt[0] = (q[2] - (l1/l2) * cos(q[0]-q[1]) * q[3])  /  (pow(l1, 2) * m1 + pow(l1, 2) * m2 * (1 - pow(cos(q[0]-q[1]),2)));

    // second coordinate: dphi2_dt
    dq_dt[1] = ((q[3]) - m2*l1*l2* dq_dt[0] * cos(q[0] - q[1])) / (m2 * pow(l2,2));



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

    double* y = (double*)malloc(n*sizeof(double));
    y[0] = phi1;
    y[1] = phi2;
    y[2] = q1;
    y[3] = q2;

    double* phi_punkt = (double*)malloc(n*sizeof(double));
    phi_punkt[0] = 0.0;
    phi_punkt[1] = 0.0;
    phi_punkt[2] = 0.0;  // not used but ne
    phi_punkt[3] = 0.0;


    double t = 0;

    fstream file;
    file.open("FoSM2_3d.dat", ios::out);

    double x1;
    double x2;
    double y1;
    double y2;
    double phi_punkt1;
    double phi_punkt2;
    double E_kin;
    double E_pot;

    while(t < T) {

        x1 = l1*sin(y[0]);
        y1 = -l1*cos(y[0]);
        x2 = x1 + l2*sin(y[1]);
        y2 = y1 - l2*cos(y[1]);
        phi_punkt1 = phi_punkt[0];
        phi_punkt2 = phi_punkt[1];

        // calculate kinetic (T) energy and potential (V) energy of the system
        E_kin = 0.5 * ((m1+m2)*l1*l1*phi_punkt1*phi_punkt1 + 2*m2*l1*l2*cos(y[0]-y[1])*phi_punkt1*phi_punkt2 +
                m2*l2*l2*phi_punkt2*phi_punkt2);
        E_pot = -g*((m1+m2)*l1*cos(y[0]) + m2*l2*cos(y[1]));

        file << t <<" "<< x1 <<" "<< y1 <<" "<< x2 <<" "<< y2 <<" "<< E_kin <<" "<< E_pot <<endl;

        rk4(dgl, y, t, dt, n);
        dgl(y, t, phi_punkt, n); //call the dgl function to get phi_punkt 1&2 out of it.

        t += dt;
    }

    file.close();

    return 0;
}
