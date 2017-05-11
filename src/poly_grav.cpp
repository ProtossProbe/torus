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
#include <eigen3/Eigen/Eigenvalues>
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

void PolyGrav::export_3d_txt(string dir, char acc = 'f') {
    ofstream txtfile;
    txtfile.open(dir);
    txtfile << vert_n << endl;
    txtfile << endl;
    for (size_t i = 0; i < vert_n; i++) {
        for (size_t j = 0; j < 3; j++) {
            if (acc == 'f') {
                txtfile << float(points[i](j)) << ' ';

            } else if (acc == 'd') {
                txtfile << setprecision(8) << points[i](j) << ' ';
            }
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
    // cout << "Initialization Completed!\n\n";
    // cout << "Vertex Number: " << vert_n << "\n\n"
    //      << "Faces Number: " << face_n << "\n\n";
    // cout << "Center of Mass: \n"
    //      << mc << "\n\n"
    //      << "Inertia Tensor: \n"
    //      << jj << "\n\n";
}

void PolyGrav::principle_axes() {
    EigenSolver<MatrixXd> es(jj);
    abc = es.eigenvalues().real();
    rotmat = es.eigenvectors().real();
    for (auto it = points.begin(); it != points.end(); ++it) {
        *it = *it - mc;
        *it = rotmat.transpose() * *it;
    }
}

double PolyGrav::L_e(double a, double b, double e) {
    return log((a + b + e) / (a + b - e));
}

double PolyGrav::S_j(double c1, double c2, double c3) {
    double c = (c3 - c1 * c2) / (sqrt(1 - c1 * c1) * sqrt(1 - c2 * c2));
    return acos(c);
}

double PolyGrav::ccos(Vector3d a, Vector3d b) {
    double result = a.dot(b);
    result /= a.norm();
    result /= b.norm();
    return result;
}

double PolyGrav::potential(Vector3d field_p) {
    double result = 0, E_term = 0, F_term = 0;
    for (auto polygon : polygons) {
        Vector3d normal;
        Matrix3d edges, edges_n, r, p, temp, E_e, F_f;
        p << points[polygon(0)], points[polygon(1)], points[polygon(2)];
        temp << field_p, field_p, field_p;
        r = p - temp;
        edges << p.col(1) - p.col(0), p.col(2) - p.col(1), p.col(0) - p.col(2);
        normal = edges.col(0).cross(edges.col(1));

        edges_n << edges.col(0).cross(normal), edges.col(1).cross(normal),
            edges.col(2).cross(normal);
        for (size_t i = 0; i < 3; i++) {
            E_e = normal * edges_n.col(i).transpose();
            E_term += r.col(i).dot(E_e * r.col(i)) *
                      PolyGrav::L_e(r.col(i).norm(), r.col((i + 1) % 3).norm(),
                                    edges.col(i).norm());
        }
        Matrix3d F_e = normal * normal.transpose();
        double c1 = PolyGrav::ccos(r.col(0), r.col(1));
        double c2 = PolyGrav::ccos(r.col(1), r.col(2));
        double c3 = PolyGrav::ccos(r.col(2), r.col(0));
        double S = PolyGrav::S_j(c1, c2, c3) + PolyGrav::S_j(c2, c3, c1) +
                   PolyGrav::S_j(c3, c1, c2);
        int sgn;
        sgn = (r.col(0).dot(normal) > 0) ? 1 : -1;
        double omega_f = (S - M_PI) * sgn;
        F_term += r.col(0).dot(F_e * r.col(0)) * omega_f;
        // cout << E_term << endl << F_term << endl;
    }
    result = co * (E_term - F_term);
    return result;
}