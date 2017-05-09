//
//  torus.hpp
//
//
//  Created by Protoss Probe on 2017/04/09.
//  Copyright © 2016-2017年 probe. All rights reserved.
//

#ifndef _POLY_GRAV_HPP_
#define _POLY_GRAV_HPP_

#include "fukushima/elliptic_integral.hpp"
#include "poly_grav.hpp"
#include "torus.hpp"
#include <boost/array.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/heuman_lambda.hpp>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <string>
#include <vector>

typedef Eigen::Matrix<double, 3, 1> pos_vec;
typedef Eigen::Matrix<size_t, 3, 1> connect_vec3;
typedef Eigen::Matrix<size_t, 4, 1> connect_vec4;
typedef std::vector<pos_vec> points_data;
typedef std::vector<connect_vec3> polygons3_data;
typedef std::vector<connect_vec4> polygons4_data;

class Torus;
class PolyGrav {
  public:
    PolyGrav();
    ~PolyGrav();
    size_t vert_n, edge_n, face_n;
    points_data points;
    polygons3_data polygons;
    void import_3d_obj(std::string dir);
};

#endif