//
//  poly_grav.cpp
//
//
//  Created by Protoss Probe on 2017/05/09.
//  Copyright © 2016-2017年 probe. All rights reserved.
//

#include <boost/array.hpp>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string>

#include "fukushima/elliptic_integral.hpp"
#include "poly_grav.hpp"
#include "torus.hpp"

using namespace std;
using namespace Eigen;

PolyGrav::PolyGrav() = default;
PolyGrav::~PolyGrav() = default;

void PolyGrav::import_3d_obj(string dir) {
    string data;
    pos_vec temp;
    connect_vec3 temp_c;
    ifstream objfile(dir);
    if (objfile.is_open()) {
        while (!objfile.eof()) {
            objfile >> data;
            if (data == "v") {
                for (size_t i = 0; i < 3; i++) {
                    objfile >> data;
                    temp(i) = stof(data);
                }
                points.push_back(temp);
            } else if (data == "f") {
                for (size_t i = 0; i < 3; i++) {
                    objfile >> data;
                    data.resize(data.size() - 2);
                    temp_c(i) = stoi(data) - 1;
                }
                polygons.push_back(temp_c);
            }
        }
        objfile.close();
    }
    vert_n = points.size();
    face_n = polygons.size();
    edge_n = vert_n + face_n - 2;
};
