#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <stdint.h>
using namespace std;

int main() {

    ifstream fin("numbers.dat");

    int32_t ii;

    double d1 = 0, d2;
    double sumup;

    fin >> ii;
    fin >> d1;
    fin >> d2;

    cout << ii  << endl;
    cout << d1  << endl;
    cout << d2  << endl;


    for (int i = 1; i <= 1000000; i++) {
        fin >> d1;
        sumup +=d1;
        //cout << d1<< endl;
        d1 = 0;
    }

    //cout << sumup<< endl;

}