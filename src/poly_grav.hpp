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

typedef boost::array<double, 3> pos3;
typedef std::vector<Eigen::Vector3d> points_data;
typedef std::vector<Eigen::Matrix<size_t, 3, 1>> polygons3_data;
typedef std::vector<Eigen::Matrix<size_t, 4, 1>> polygons4_data;

class Torus;
class PolyGrav {
  public:
    PolyGrav();
    ~PolyGrav();
    PolyGrav(std::string dir);
    std::string dir;
    size_t vert_n, edge_n, face_n;
    Eigen::Vector3d mc;
    Eigen::Matrix3d jj;
    points_data points;
    polygons3_data polygons;
    void init();

  private:
    void import_3d_obj(std::string dir);
    void export_3d_txt(std::string dir);
    void calexec(std::string dir);
    void import_info(std::string dir);
};

#endif