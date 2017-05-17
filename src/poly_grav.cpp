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
        a[i] = a1[i] - a2[i];
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
    vec3 temp, normal, edge_len;
    mat3 p, edge, F_f;
    ten3 E_e;
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
                    temp_c[i] = stoi(data) - 1;
                }

                p = {{points[temp_c[0]], points[temp_c[1]], points[temp_c[2]]}};
                edge = {{p[1] - p[0], p[2] - p[1], p[0] - p[2]}};
                edge_len = {{PolyGrav::norm(edge[0]), PolyGrav::norm(edge[1]),
                             PolyGrav::norm(edge[2])}};
                for (size_t i = 0; i < 3; i++) {
                    edge[i] = edge[i] / edge_len[i];
                }

                normal = PolyGrav::cross(edge[0], edge[1]);
                normal = normal / PolyGrav::norm(normal);
                for (size_t i = 0; i < 3; i++) {
                    E_e[i] = PolyGrav::outer(normal,
                                             PolyGrav::cross(edge[i], normal));
                }
                F_f = PolyGrav::outer(normal, normal);

                polygons.push_back(temp_c);
                edges.push_back(edge);
                edges_len.push_back(edge_len);
                normals.push_back(normal);
                E_es.push_back(E_e);
                F_fs.push_back(F_f);
            }
        }
        objfile.close();
        vert_n = points.size();
        face_n = polygons.size();
        edge_n = vert_n + face_n - 2;
    };
}

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
    Matrix3d temp_mat;
    temp_mat = boost2eigen_mat(jj);
    EigenSolver<MatrixXd> es(temp_mat);
    abc = eigen2boost_vec(es.eigenvalues().real());
    mat3 rotmat = eigen2boost_mat(es.eigenvectors().real());
    for (auto it = points.begin(); it != points.end(); ++it) {
        *it = *it - mc;
        *it = PolyGrav::mul(PolyGrav::transpose(rotmat), *it);
    }
}

double PolyGrav::L_e(double a, double b, double e) {
    return log((a + b + e) / (a + b - e));
}

double PolyGrav::omega_f(mat3 r, vec3 r_len) {
    double y = PolyGrav::dot(r[0], PolyGrav::cross(r[1], r[2]));
    double x = 0;
    x += r_len[0] * r_len[1] * r_len[2];
    for (size_t i = 0; i < 3; i++) {
        x += r_len[i] * PolyGrav::dot(r[(i + 1) % 3], r[(i + 2) % 3]);
    }

    return 2 * atan2(y, x);
}

double PolyGrav::potential(vec3 field_p) {
    double result = 0, E_term = 0, F_term = 0;
    for (size_t n = 0; n < face_n; n++) {
        mat3 r, p, temp;
        vec3 r_len;
        auto polygon = polygons[n];
        auto edge_len = edges_len[n];
        auto E_e = E_es[n];
        auto F_f = F_fs[n];
        auto normal = normals[n];

        temp = {{field_p, field_p, field_p}};
        p = {{points[polygon[0]], points[polygon[1]], points[polygon[2]]}};
        r = p - temp;
        r_len = {
            {PolyGrav::norm(r[0]), PolyGrav::norm(r[1]), PolyGrav::norm(r[2])}};

        for (size_t i = 0; i < 3; i++) {
            E_term += PolyGrav::dot(r[i], PolyGrav::mul(E_e[i], r[i])) *
                      PolyGrav::L_e(r_len[i], r_len[(i + 1) % 3], edge_len[i]);
        }

        double omega_f = PolyGrav::omega_f(r, r_len);
        F_term += PolyGrav::dot(r[0], PolyGrav::mul(F_f, r[0])) * omega_f;
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
    return output;
}

mat3 PolyGrav::transpose(const mat3 &mat) {
    mat3 output;
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
            output[i][j] = mat[j][i];
        }
    }
    return output;
}