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
    size_t limit = 128;
    size_t m_deg = 10; // sets (2m+1)-point Gauss-Kronrod

    Workspace<Real> Work(limit, m_deg);

    // Set parameters and function class...
    Real alpha = (argc > 1) ? Real(atof(argv[1])) : Real(1);

    Function<Real, Real> F(integrand, &alpha);

    // Set default quadrature tolerance from machine epsilon...
    Real epsabs = Work.get_eps() * 1e2;
    Real epsrel = Real(0);
    std::cout << "# epsabs : " << epsabs << std::endl;

    int status = 0;
    Real result, abserr;

    try
    {
        status = Work.qag(F, Real(0), Real(1), epsabs, epsrel, result, abserr);
    }
    catch (const char *reason)
    {
        std::cerr << reason << std::endl;
        return status;
    }

    std::cout << "result: " << result;
    std::cout << " error: " << ABS(result - exact(alpha));
    std::cout << std::endl;

    std::cout << Real(1) << endl;

    return 0;
}
