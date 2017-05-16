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

vec3 operator+(const vec3 &a1, const vec3 &a2) {
    vec3 a;
    for (size_t i = 0; i < 3; i++)
        a[i] = a1[i] + a2[i];
    return a;
}

vec3 operator-(const vec3 &a1, const vec3 &a2) {
    vec3 a;
    for (size_t i = 0; i < 3; i++)
        a[i] = a1[i] + a2[i];
    return a;
}

vec3 operator*(const vec3 &a1, const double &a2) {
    vec3 a;
    for (size_t i = 0; i < 3; i++)
        a[i] = a1[i] * a2;
    return a;
}

vec3 operator/(const vec3 &a1, const double &a2) {
    vec3 a;
    for (size_t i = 0; i < 3; i++)
        a[i] = a1[i] / a2;
    return a;
}

mat3 operator+(const mat3 &a1, const mat3 &a2) {
    mat3 a;
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            a[i][j] = a1[i][j] + a2[i][j];
        }
    }
    return a;
}

mat3 operator-(const mat3 &a1, const mat3 &a2) {
    mat3 a;
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            a[i][j] = a1[i][j] - a2[i][j];
        }
    }
    return a;
}

PolyGrav::PolyGrav() = default;
PolyGrav::~PolyGrav() = default;
PolyGrav::PolyGrav(string dir) : dir(dir){};

void PolyGrav::import_3d_obj(string dir) {
    // parse .obj file to get vertexs and polygon data
    string data;
    vec3 temp;
    connect3 temp_c;
    ifstream objfile(dir);
    if (objfile.is_open()) {
        while (!objfile.eof()) {
            objfile >> data;
            if (data == "v") {
                for (size_t i = 0; i < 3; i++) {
                    objfile >> data;
                    temp[i] = stof(data);
                }
                points.push_back(temp);
            } else if (data == "f") {
                for (size_t i = 0; i < 3; i++) {
                    objfile >> data;
                    data.resize(data.size() - 2);
                    temp_c[i] = stoi(data) - 1;
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
                txtfile << float(points[i][j]) << ' ';

            } else if (acc == 'd') {
                txtfile << setprecision(8) << points[i][j] << ' ';
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
            txtfile << polygons[i][j] << ' ';
        }
        txtfile << endl;
    }
}

void PolyGrav::import_info(string dir) {
    ifstream txtfile(dir);
    if (txtfile.is_open()) {
        while (!txtfile.eof()) {
            for (size_t i = 0; i < 3; i++) {
                txtfile >> mc[i];
            }
            for (size_t i = 0; i < 3; i++) {
                for (size_t j = 0; j < 3; j++) {
                    txtfile >> jj[i][j];
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
    Matrix3d temp_mat, temp_mat2, rotmat;
    temp_mat = boost2eigen_mat(jj);
    EigenSolver<MatrixXd> es(temp_mat);
    temp_mat2 = es.eigenvalues().real();
    abc = eigen2boost_vec(temp_mat2);
    rotmat = es.eigenvectors().real();
    for (auto it = points.begin(); it != points.end(); ++it) {
        *it = *it - mc;
        *it = eigen2boost_vec(rotmat.transpose() * boost2eigen_vec(*it));
    }
}

double PolyGrav::L_e(double a, double b, double e) {
    return log((a + b + e) / (a + b - e));
}

double PolyGrav::ccos(vec3 a, vec3 b) {
    double result = a.dot(b);
    result /= a.norm();
    result /= b.norm();
    return result;
}

double PolyGrav::potential(vec3 field_p) {
    double result = 0, E_term = 0, F_term = 0;
    int n = 50, i = 0;
    for (auto polygon : polygons) {
        vec3 normal;
        mat3 edges, edges_n, r, p, temp, E_e, F_f;
        p = {{points[polygon[0]], points[polygon[1]], points[polygon[2]]}};
        temp = {{field_p, field_p, field_p}};
        r = p - temp;
        edges = {{p[1] - p[0], p[2] - p[1], p[0] - p[2]}};
        for (size_t i = 0; i < 3; i++) {
            edges[i] = edges[i] / PolyGrav::norm(edges[i]);
        }
        normal = PolyGrav::cross(edges[0], edges[1]);

        edges_n = {{PolyGrav::cross(edges[0], normal),
                    PolyGrav::cross(edges[1], normal),
                    PolyGrav::cross(edges[2], normal)}};

        for (size_t i = 0; i < 3; i++) {
            E_e = PolyGrav::outer(normal, edges_n[i]);
            E_term += r.col(i).dot(E_e * r.col(i)) *
                      PolyGrav::L_e(r.col(i).norm(), r.col((i + 1) % 3).norm(),
                                    edges.col(i).norm());
        }
        mat3 F_e = normal * normal.transpose();
        double c1 = PolyGrav::ccos(r.col(0), r.col(1));
        double c2 = PolyGrav::ccos(r.col(1), r.col(2));
        double c3 = PolyGrav::ccos(r.col(2), r.col(0));
        double c4 = r.col(0).dot(r.col(1).cross(r.col(2))) /
                    (r.col(0).norm() * r.col(1).norm() * r.col(2).norm());
        double cc = (1 + c1) * (1 + c2) * (1 + c3);

        double cosomega = 1 - c4 * c4 / cc;
        double sinomega = (1 + c1 + c2 + c3) / cc * abs(c4);
        int sgn;
        sgn = (r.col(0).dot(normal) > 0) ? 1 : -1;
        double omega_f = 2 * atan2(1 - cosomega, sinomega) * sgn;
        F_term += r.col(0).dot(F_e * r.col(0)) * omega_f;

        if (i % n == 0) {
            // cout << E_term << endl;
            // cout << F_term << endl;
            // cout << r.col(0).dot(F_e * r.col(0)) << endl;
            // cout << cosomega << '\t' << sinomega << endl;
            // cout << omega_f << endl;
            // cout << endl;
        }

        // cout << setprecision(8) << E_term << endl << F_term << endl;
        i++;
    }
    result = co * (E_term - F_term);
    return result;
}

Vector3d PolyGrav::boost2eigen_vec(vec3 vec) {
    Vector3d output;
    for (size_t i = 0; i < 3; i++) {
        output(i) = vec[i];
    }
    return output;
}

Matrix3d PolyGrav::boost2eigen_mat(mat3 mat) {
    Matrix3d output;
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            output(i, j) = mat[i][j];
        }
    }
    return output;
}

vec3 PolyGrav::eigen2boost_vec(Vector3d vec) {
    vec3 output;
    for (size_t i = 0; i < 3; i++) {
        output[i] = vec(i);
    }
    return output;
}

mat3 PolyGrav::eigen2boost_mat(Matrix3d mat) {
    mat3 output;
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            output[i][j] = mat(i, j);
        }
    }
    return output;
}

double PolyGrav::norm(const vec3 &vec) {
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

double PolyGrav::dot(const vec3 &vec1, const vec3 &vec2) {
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

vec3 PolyGrav::cross(const vec3 &vec1, const vec3 &vec2) {
    vec3 output;
    output[0] = -vec1[2] * vec2[1] + vec1[1] * vec2[2];
    output[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    output[2] = -vec1[1] * vec2[0] + vec1[0] * vec2[1];
    return output;
}

vec3 PolyGrav::mul(const mat3 &mat, const vec3 &vec) {
    vec3 output;
    for (size_t i = 0; i < 3; i++) {
        output[i] =
            mat[i][0] * vec[0] + mat[i][1] * vec[1] + mat[i][2] * vec[2];
    }
    return output;
}

mat3 PolyGrav::outer(const vec3 &vec1, const vec3 &vec2) {
    mat3 output;
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            output[i][j] = vec1[i] * vec2[j];
        }
    }
}