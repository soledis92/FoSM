#pragma once
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

//de(double*,double,double*,int) = differential equation with arguments (q-vector , t-time, dy_dt-time derivative of y, n-dimensions)


void rk2( void(*de)(vector<double>&, double, vector<double>&, int), vector<double>& q, double t, double dt, int n ) {
    int i, k, p, PM = 2;
    static const double c2 = 1./2, c1 = 1.;

    vector< vector<double> > a;
    a.resize(PM, vector<double> (PM, 0.0));

    //init a  Matrix with RK-coefficients
    a[0][0] = c2; a[0][1] = 0.;
    a[1][0] = 0.; a[1][1] = c1;

    vector<double> b (PM, 0.0);
    //init b Vector with RK-coefficients
    b[0] = 0.; b[1] = c1;

    vector< vector<double> > y;
    y.resize(PM+1, vector<double> (PM+1, 0.0));

    vector< vector<double> > z;
    z.resize(PM, vector<double> (PM, 0.0));

    for(i = 0; i < n; i++)
        y[0][i] = q[i];
    for(p = 1; p <= PM; p++) {
        de(y[p-1], t+b[p-1]*dt, z[p-1], n);
        for(i = 0; i < n; i++)
            y[p][i] = q[i];
        for(k = 0; k < p; k++)
            for(i = 0; i < n; i++)
                y[p][i] = y[p][i] + dt*a[p-1][k]*z[k][i];
    }
    for(i = 0; i < n; i++)
        q[i] = y[PM][i];
}
