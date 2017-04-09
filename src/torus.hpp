//
//  torus_gravity.h
//  
//
//  Created by Protoss Probe on 2017/04/09.
//  Copyright © 2016-2017年 probe. All rights reserved.
//

#ifndef _TORUS_HPP_
#define _TORUS_HPP_

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/heuman_lambda.hpp>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/tools/roots.hpp>

#include <gsl/gsl_integration.h>

typedef boost::array<double, 2> state_type2;
typedef boost::array<double, 5> state_type5;
typedef boost::array<double, 6> state_type6;
typedef boost::array<double, 21> state_type21;
typedef boost::array<double, 23> state_type23;

typedef runge_kutta_dopri5<double, double, double, double, vector_space_algebra, default_operations, never_resizer> stepper_type;

#endif