//
//  poly_grav.hpp
//
//
//  Created by Protoss Probe on 2017/05/09.
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
#include <eigen3/Eigen/Eigenvalues>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <string>
#include <vector>

typedef boost::array<double, 3> vec3;
typedef boost::array<vec3, 3> mat3;
typedef boost::array<size_t, 3> connect3;
typedef boost::array<size_t, 4> connect4;
typedef std::vector<vec3> points_data;
typedef std::vector<connect3> polygons3_data;
typedef std::vector<connect4> polygons4_data;

class Torus;
class PolyGrav {
  public:
    PolyGrav();
    ~PolyGrav();
    PolyGrav(std::string dir);
    std::string dir;
    size_t vert_n, edge_n, face_n;
    double co = 0.5;
    vec3 mc;
    vec3 abc;
    mat3 jj;
    mat3 rotmat;
    points_data points;
    polygons3_data polygons;
    void init();
    void principle_axes();
    void export_3d_txt(std::string dir, char acc);
    double potential(vec3 field_p);

  private:
    void import_3d_obj(std::string dir);
    void calexec(std::string dir);
    void import_info(std::string dir);
    double L_e(double a, double b, double e);
    double ccos(vec3 a, vec3 b);
    Eigen::Vector3d boost2eigen_vec(vec3 vec);
    Eigen::Matrix3d boost2eigen_mat(mat3 mat);
    vec3 eigen2boost_vec(Eigen::Vector3d vec);
    mat3 eigen2boost_mat(Eigen::Matrix3d mat);
    double norm(const vec3 &vec);
    double dot(const vec3 &vec1, const vec3 &vec2);
    vec3 cross(const vec3 &vec1, const vec3 &vec2);
    vec3 mul(const mat3 &mat, const vec3 &vec);
    mat3 outer(const vec3 &vec1, const vec3 &vec2);
};

#endif