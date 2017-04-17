#ifndef _ELLIPTIC_INTETRAL_HPP_
#define _ELLIPTIC_INTETRAL_HPP_

#include <vector>
#include <boost/array.hpp>
#include "../torus.hpp"

typedef std::vector<boost::array<double, 2>> container;
typedef boost::array<double, 4> array_4;
typedef boost::array<double, 2> array_2;

namespace Elliptic_Integral
{
double ceik(double m);
double ceie(double m);
double ceib(double m);
double ceid(double m);
double ceis(double m);
container zonal_toroidal_harmonics_scale(double u, size_t iter);
array_4 zonal_toroidal_harmonics_seed(double u);
}


#endif