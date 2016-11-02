#include <iostream>
#include <cmath>

void f_D_2(double ** ptr, int row, int col){
    for(int i = 0; i < row; i++){
        free(ptr[i]);
    }
    free(ptr);
}

void f_D_1(double * ptr, int row){
    free(ptr);
}


void rk2(void(*de)(double*,double,double*,int), double* q, double t, double dt, int n) {
    int i, k, p, PM = 2;
    static const double c2 = 1./2, c1 = 1.;
    
    double** a = (double**)malloc(PM*sizeof(double*));
    for(i = 0; i < PM; i++)
    a[i] = (double*)malloc(PM*sizeof(double));
    
    //init a
    a[0][0] = c2; a[0][1] = 0.;
    a[1][0] = 0.; a[1][1] = c1;
    
    double* b = (double*)malloc(PM*sizeof(double));
    //init b
    b[0] = 0.; b[1] = c1;
    
    double** y = (double**)malloc((PM+1)*sizeof(double*));
    for(i = 0; i < PM+1; i++)
    y[i] = (double*)malloc(n*sizeof(double));
    
    double** z = (double**)malloc(PM*sizeof(double*));
    for(i = 0; i < PM; i++)
    z[i] = (double*)malloc(n*sizeof(double));
    
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
    
    f_D_2(a,PM,PM); f_D_1(b,PM); f_D_2(y,PM+1,n); f_D_2(z,PM,n);
}

