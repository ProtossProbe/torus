#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <numeric>
#include <boost/array.hpp>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <cstdlib>
#include <quadpackpp/workspace.hpp>

using namespace Eigen;
using namespace std;
/* logarithmic.cpp
 */

typedef double Real;

#define ABS(x) (((x) < Real(0)) ? -(x) : (x))

Real integrand(Real x, Real *alpha)
{
    return pow(x, *alpha) * log(1 / x);
}

Real exact(Real alpha)
{
    return pow(alpha + 1, -2);
}

int main(int argc, char **argv)
{


    return 0;
}
