#include <iostream>

using namespace std;

int main() {

    //---------------------------------------------------
    // Machine Epsilon
    //-------------------------------------------

    float a = 1;

    float e = 1;

    float h = 1.1;

    while ((a + e) > a) {
        e = e / h;
    }

    cout << (e) << endl;


    double b = 1;

    double eb = 1;

    double hb = 1.1;

    while ((b + eb) > b) {
        eb = eb / hb;
    }
    cout << (eb) << endl;

    long double c = 1;

    long double ec = 1;

    long double hc = 1.1;

    while ((c + ec) > c){
        ec = ec/hc;
    }

    cout<< (ec) << endl;

    //---------------------------------------------------
    // Pitfalls of floating point arithmetic
    //-------------------------------------------

    double aa = 1.0e17;
    double bb = -1.0e17;
    double cc = 1.0;
    double x = (aa + bb) + cc;
    double y = aa + (bb + cc);


    cout << "x = " << x << ", " << "y = " << y << endl;

    //---------------------------------------------------
    // Pitfalls of floating point representation
    //-------------------------------------------

    // a)

    float xx = 0.01;
    double yy = xx;
    double zz = 0.01;

    cout << "x = " << xx << ", " << "y = " << yy << ", z = " << zz << endl;

    int i = xx * 10000;
    int j = yy * 10000;
    int k = zz * 10000;

    printf("%d %d %d\n", i, j, k);


    return 0;
}
