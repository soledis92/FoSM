#include <iostream>
#include <bitset>

using namespace std;
int main()
{
    float x = 0.01;
    double f = x;
    cout<<bitset<sizeof f*8>(*(long unsigned int*)(&f))<<endl;


    float f2 = 0.01;
    cout<<bitset<sizeof f2*8>(*(long unsigned int*)(&f2))<<endl;
    return 0;
}
