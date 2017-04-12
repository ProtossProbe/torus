//
//  torus_gravity.h
//  
//
//  Created by Protoss Probe on 2017/04/09.
//  Copyright © 2016-2017年 probe. All rights reserved.
//

#ifndef _TORUS_HPP_
#define _TORUS_HPP_


#include <iostream>
#include <boost/array.hpp>
#include <string>
#include <cmath>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/heuman_lambda.hpp>
#include <gsl/gsl_integration.h>
#include "torus.hpp"

typedef boost::array<double, 2> pos;
typedef boost::array<double, 3> vel;
typedef boost::array<double, 5> state_type5;

struct f_params
{
    double r;
    double x3;
    double r0;
};

struct f_params2
{
    double r;
    double x3;
    double eta;
    double r0;
};

const double pi = 3.1415926535897932384626433832795028841971693993751;

#endif