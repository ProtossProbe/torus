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
PolyGrav::PolyGrav(string dir) : dir(dir){};

void PolyGrav::import_3d_obj(string dir) {
    // parse .obj file to get vertexs and polygon data
    string data;
    Vector3d temp;
    Matrix<size_t, 3, 1> temp_c;
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

void PolyGrav::export_3d_txt(string dir) {
    ofstream txtfile;
    txtfile.open(dir);
    txtfile << vert_n << endl;
    txtfile << endl;
    for (size_t i = 0; i < vert_n; i++) {
        for (size_t j = 0; j < 3; j++) {
            txtfile << float(points[i](j)) << ' ';
        }
        txtfile << endl;
    }
    txtfile << endl;
    txtfile << face_n << endl;
    txtfile << endl;
    for (size_t i = 0; i < face_n; i++) {
        txtfile << "3 ";
        for (size_t j = 0; j < 3; j++) {
            txtfile << polygons[i](j) << ' ';
        }
        txtfile << endl;
    }
}

void PolyGrav::import_info(string dir) {
    ifstream txtfile(dir);
    if (txtfile.is_open()) {
        while (!txtfile.eof()) {
            for (size_t i = 0; i < 3; i++) {
                txtfile >> mc(i);
            }
            for (size_t i = 0; i < 3; i++) {
                for (size_t j = 0; j < 3; j++) {
                    txtfile >> jj(i, j);
                }
            }
        }
    }
}

void PolyGrav::calexec(string dir) {
    string exe = "../bin/volInt ";
    exe += dir;
    const char *input = exe.c_str();
    // execl(executable, executable, input, (char *)NULL);

    system(input);
}

void PolyGrav::init() {
    string filename = dir;
    PolyGrav::import_3d_obj("../assets/" + filename + ".obj");
    PolyGrav::export_3d_txt("../assets/" + filename + ".txt");
    PolyGrav::calexec("../assets/" + filename + ".txt");
    PolyGrav::import_info("../assets/" + filename + "_info.txt");
    cout << "Initialization Completed!\n\n";
    cout << "Vertex Number: " << vert_n << "\n\n"
         << "Faces Number: " << face_n << "\n\n";
    cout << "Center of Mass: \n"
         << mc << "\n\n"
         << "Inertia Tensor: \n"
         << jj << "\n\n";
}
