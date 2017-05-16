#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <fftw3.h>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;
typedef boost::array<double, 3> vec3;
typedef boost::array<vec3, 3> mat3;

template <class T> T operator+(const T &a1, const T &a2) {
    T a;
    for (typename T::size_type i = 0; i < a1.size(); i++)
        a[i] = a1[i] + a2[i];
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

vec3 mul(const mat3 &mat, const vec3 &vec) {
    vec3 output;
    for (size_t i = 0; i < 3; i++) {
        output[i] =
            mat[i][0] * vec[0] + mat[i][1] * vec[1] + mat[i][2] * vec[2];
    }
    return output;
}

int main() {
    vec3 a = {1, 2, 4};
    vec3 b = {3, 4, 5};
    vec3 c = a + b;
    cout << c[0] << endl;

    mat3 a1 = {{{1, 2, 4}, {2, 3, 4}, {3, 4, 5}}};
    mat3 b1 = {{{1, 2, 4}, {2, 3, 4}, {3, 4, 5}}};
    mat3 c1 = {{a, a, b}};

    vec3 test = mul(a1, a);
    for (auto ele : c1) {
        for (auto el : ele) {
            cout << el << endl;
        }
    }
}