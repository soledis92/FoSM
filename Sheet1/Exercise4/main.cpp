#include <iostream>

int main(){

    float x = 1;

    float e = 1;

    float h = 1.1;

    while ((x + e) > x){
        e = e/h;

    }
    std::cout<< (e) << std::endl;


    float y = 255;

    float e2 = 1;

    float h2 = 1.1;

    while ((y + e2) > y){
        e2 = e2/h2;
    }
    std::cout<< (e2) << std::endl;


    double xx = 1;

    double ee = 1;

    double hh = 1.1;

    while ((xx + ee) > xx){
        ee = ee/hh;

    }
    std::cout<< (ee) << std::endl;


    double yy = 255;

    double ee2 = 1;

    double hh2 = 1.1;

    while ((yy + ee2) > yy){
        ee2 = ee2/hh2;
    }
    std::cout<< (ee2) << std::endl;


    return 0;
}