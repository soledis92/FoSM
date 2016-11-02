#include <iostream>
#include <cmath>
#include <vector>
#include "rk2.hh"
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

void dgl(vector<double>& q, double t, vector<double>& dq_dt, int dim){

    //first coordinate:  dphi1_dt
    dq_dt[0] = (q[2] / (pow(l1,2) * (m1 + m2)) - (q[3] * cos(q[1] - q[2])) / (l1 * l2 * (m1 +m2)))  *
            (1 / (1 - (m1 * pow(cos(q[1] - q[2]), 2) / (m1 + m2) )));

    // second coordinate: dphi2_dt
    dq_dt[1] = ((q[2]) - m2*l1*l2* dq_dt[0] * cos(q[1] - q[2])) / (m2 * pow(l2,2));

    // third coordinate: dq1_dt (NICHT ZU VERWECHSELN MIT DEM VEKTOR dq_dt HIER !!!!!
    dq_dt[2] = - m2*l1*l2* dq_dt[0] * dq_dt[1] * sin(q[1] - q[2]) - (m1+m2)*g*l1 * sin(q[0]);
}



int main(){



    // Initial values

    double phi1 = 50.0 / 180.0 * M_PI;
    double phi2 = -120 / 180.0 * M_PI;
    double dphi1_dt = 0.0;
    double dphi2_dt = 0.0;
    double q1 = 0.0;
    double q2 = 0.0;

    // Initial value vector y0

    vector<double> y0 (4, 0.0);
    y0[0] = phi1;
    y0[1] = phi2;
    y0[2] = q1;
    y0[3] = q2;

    double t = 0;

    for(int i = 0; i < T; i++) {
        rk2(dgl, y0, t, dt, y0.size());
        t += dt;

        cout << A[0] <<endl;
    }



    return 0;
}
