#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <fftw3.h>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;
typedef boost::array<double, 3> vec3;
typedef boost::array<vec3, 3> mat3;

size_t test(int i)
{
    if (i == 2)
    {
        return 2;
    }
    else
    {
        return -1;
    }
}

int main()
{
    if (test(2))
    {
        cout << test(1) << endl;
        cout << test(2) << endl;
    }
}