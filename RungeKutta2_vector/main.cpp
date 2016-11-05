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

void dgl2(vector<double>& q, double t, vector<double>& dq_dt, int dim){
    dq_dt[0] = q[0];
    //dq_dt[1] = 0;
    //dq_dt[2] = 0;
    //dq_dt[3] = 0;


}


// define the differential equation which is has to be solved

void dgl(vector<double>& q, double t, vector<double>& dq_dt, int dim){

    //first coordinate:  dphi1_dt
    dq_dt[0] = (q[2] / (pow(l1,2) * (m1 + m2)) - (q[3] * cos(q[0] - q[1])) / (l1 * l2 * (m1 +m2)))  *
            (1. / (1. - (m1 * pow(cos(q[0] - q[1]), 2) / (m1 + m2) )));
    cout << "dq_dt:  "<< dq_dt[0]<<endl;

    // second coordinate: dphi2_dt
    dq_dt[1] = ((q[2]) - m2*l1*l2* dq_dt[0] * cos(q[0] - q[1])) / (m2 * pow(l2,2));
    cout << dq_dt[1]<<endl;


    // third coordinate: dq1_dt  (NICHT ZU VERWECHSELN MIT DEM VEKTOR dq_dt HIER !!!!! dq1_dt = dq_dt[2] )
    dq_dt[2] = - m2*l1*l2* dq_dt[0] * dq_dt[1] * sin(q[0] - q[1]) - (m1+m2)*g*l1 * sin(q[0]);

    cout << dq_dt[2]<<endl;

    // fourth coordinate: dq2_dt (NICHT ZU VERWECHSELN MIT DEM VEKTOR dq_dt HIER !!!!! dq1_dt = dq_dt[3] )
    dq_dt[3] = m2*l1*l2 * dq_dt[0] * dq_dt[1] * sin(q[0] - q[1]) - m2*g*l2 * sin(q[1]);
    cout << dq_dt[3]<<endl;

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

    vector<double> y0 (1, 0.0);
    y0[0] = phi1;
    //y0[1] = phi2;
    //y0[2] = q1;
    //y0[3] = q2;

    double t = 0;


    for(int i = 0; i < 21; i++) {
        rk2(dgl2, y0, t, dt, 1);
        t += dt;

        cout << y0[0] <<endl;

    }

    vector< vector<double> > r;
    r.resize(3, vector<double> (3, 0.0));

    vector<double> b = r[0];
    cout<<"B: "<< b.size()<<endl;



    cout << "hallo world"<< endl;

    return 0;
}
