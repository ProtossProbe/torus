//
//  torus.cpp
//
//
//  Created by Protoss Probe on 2017/04/09.
//  Copyright © 2016-2017年 probe. All rights reserved.
//

#include <boost/array.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <iostream>
#include <quadpackpp/workspace.hpp>
#include <string>

#define EIGEN_STACK_ALLOCATION_LIMIT 10000000
#include "fukushima/elliptic_integral.hpp"
#include "torus.hpp"
#include <eigen3/Eigen/Dense>
#include <fftw3.h>

using namespace std;
using namespace boost::math;
using namespace boost::numeric::odeint;
using namespace Eigen;
using namespace Elliptic_Integral;

class Torus
{
  public:
    Torus() = default;
    ~Torus() = default;
    Torus(double r0) : r0(r0) {}
    Torus(double r0, double R0) : r0(r0), R0(R0) {}
    const double r0 = 0.25, R0 = 1.0;
    const double R1u = R0 - r0, R2u = R0 + r0;
    const double au = sqrt(R0 * R0 - r0 * r0);
    const double ubu = R0 / r0;
    double ur = 3.6;

    static const size_t Ns = 512;
    typedef Matrix<double, Ns, 1> sample_vec;
    typedef Matrix<double, Ns, Ns> sample_mat;

    pos r2u(const pos rz)
    {
        return {utr(rz[0], rz[1], au), thetatr(rz[0], rz[1], au)};
    }

    pos u2r(const pos ut)
    {
        return {rto(ut[0], ut[1], au), zto(ut[0], ut[1], au)};
    }

    void sampling()
    {
        double theta, v;
        sample_vec t;
        for (size_t i = 1; i <= Ns; i++)
        {
            theta = (2 * i - 1) * pi / (2 * Ns);
            v = potential_ut(ur, theta) / sqrt(dto(ur, theta)) * (-au);
            thetalist(i - 1) = theta;
            vlist(i - 1) = v;
            for (size_t j = 0; j < Ns; j++)
            {
                t(j) = cos(j * theta);
            }
            tlist.row(i - 1) = t;
        }
        // cout << "theta: " << endl
        //      << thetalist << endl
        //      << "v: " << endl
        //      << vlist << endl
        //      << "t: " << endl
        //      << tlist << endl;
    }

    void dct()
    {
        clist = vlist / (int(Ns));
        double *p = &clist[0];
        fftw_plan plan =
            fftw_plan_r2r_1d(Ns, p, p, FFTW_REDFT10, FFTW_ESTIMATE);
        fftw_execute(plan);
        clist(0) /= 2;
        // cout << setprecision(15) << "c: " << endl
        //      << clist;
    }

    void legendre_ratio(const double u, const size_t n)
    {
        plist = Elliptic_Integral::ratio_p(u, ur, n);
    }

    double potential_harmonics_ut(const double u, const double theta,
                                  const size_t n)
    {
        legendre_ratio(u, n);
        double result = 0;
        for (size_t i = 0; i < n; i++)
        {
            result += basic_func(theta, i);
        }
        result *= (-1 / au) * sqrt(dto(u, theta));
        return result;
    }

    double basic_func(const double theta, const size_t i)
    {
        return clist(i) * cos(i * theta) * plist[i][0];
    }

    double potential_harmonics_rz(const double r, const double z,
                                  const size_t n)
    {
        pos ut = r2u({r, z});
        return potential_harmonics_ut(ut[0], ut[1], n);
    }

    double potential_ut(const double u, const double theta, char method = 'a')
    {
        pos rz = u2r({u, theta});
        return potential(rz[0], rz[1], method);
    }

    double potential(const double r, const double x3, char method = 'a')
    {
        // Two different ways to calculate the potential of a torus.
        // method 'a' use a single quadrature and has high accuracy and
        // efficiency. method 'b' use a double quadrature and its accuracy and
        // efficiency is low.

        double result;
        switch (method)
        {
        case 'a':
        {
            result = potential_func(r, x3);
            break;
        }
        case 'b':
        {
            result = potential_func2(r, x3);
            break;
        }
            // case 'c':
            // {
            //     result = potential_func3(r, x3);
            //     break;
            // }
        }
        return result;
    }

    void test()
    {
        double r = 1.25, z = 0.0;
        struct f_params params
        {
            r, z, r0
        };
        cout << "value: " << func_single_quad(asin(z / r0), &params) << endl;
    }

  private:
    sample_vec thetalist, vlist;
    sample_mat tlist;
    sample_vec clist;
    container plist;

    double ftr(const double r, const double z, const double a)
    {
        return sqrt(pow((r * r - a * a), 2) + 2 * (r * r + a * a) * z * z +
                    pow(z, 4));
    }

    double utr(const double r, const double z, const double a)
    {
        return (r * r + z * z + a * a) / ftr(r, z, a);
    }

    double thetatr(const double r, const double z, const double a)
    {
        return atan2(2 * a * z, r * r + z * z - a * a);
    }

    double dtr(const double r, const double z, const double a)
    {
        return 2 * a * a / ftr(r, z, a);
    }

    double dto(const double u, const double theta) { return u - cos(theta); }

    double vto(const double u) { return sqrt(u * u - 1); }

    double rto(const double u, const double theta, const double a)
    {
        return a * vto(u) / dto(u, theta);
    }

    double zto(const double u, const double theta, const double a)
    {
        return a * sin(theta) / dto(u, theta);
    }

    double potential_func(const double r, const double x3)
    {
        // use potential_single_quad to do the quadrature.
        double u_p;
        if (abs(r) < 1e-12)
        {
            double lam, k;
            lam = 2 * r0 * sqrt(1 + pow(x3, 2)) / (1 + pow(r0, 2) + pow(x3, 2));
            k = sqrt(2 * lam / (1 + lam));
            u_p = (4.0 / 3.0 / (pi * r0) * sqrt(4 - 2 * pow(k, 2)) / pow(k, 2) *
                   sqrt(1 + pow(r0, 2) / (1 + pow(x3, 2))) *
                   (ellint_2(k) -
                    2 * (1 - pow(k, 2)) / (2 - pow(k, 2)) * ellint_1(k)));
        }
        else
        {
            struct f_params params
            {
                r, x3, r0
            };
            Function<double, void> F(func_single_quad, &params);
            u_p = quadrature(F, 0, 2 * pi, 'b');
        }
        return u_p;
    }

    static double func_single_quad(double theta, void *params)
    {
        // Integrand for potential function.
        // variable: theta
        // parameters: r, x3, r0

        struct f_params *value = (struct f_params *)params;
        double r = value->r, x3 = value->x3, r0 = value->r0;
        r = abs(r);

        double R1, p, q, k, k_p, z, phi, result, ei1, ei2, ei1_phi, ei2_phi;
        int sgn = 0;
        R1 = 1 + r0 * cos(theta);
        z = x3 - r0 * sin(theta);
        p = sqrt(pow((R1 + r), 2) + z * z);
        q = sqrt(pow((R1 - r), 2) + z * z);
        k = 2 / p * sqrt(r * R1);

        k_p = sqrt(1 - k * k);
        phi = asin(abs(z) / q);

        ei1 = ellint_1(k);
        ei2 = ellint_2(k);
        ei1_phi = ellint_1(k_p, phi);
        ei2_phi = ellint_2(k_p, phi);

        if (R1 - r > 0)
        {
            sgn = 1;
        }
        else if (R1 - r < 0)
        {
            sgn = -1;
        }

        result =
            (-p * ei2 - (R1 * R1 - r * r) / p * ei1 +
             abs(z) * (pi / 2 + pi / 2 * sgn) -
             abs(z) * sgn * (ei2 * ei1_phi + ei1 * ei2_phi - ei1 * ei1_phi));

        return -result * cos(theta) / (pi * pi * r0);
    }

    double potential_func2(const double r, const double x3)
    {
        // use func_double_quad to do a double quadature
        struct f_params params
        {
            abs(r), x3, r0
        };
        Function<double, void> F(func_double_quad, &params);
        return quadrature(F, -r0, r0, 'c') / (pi * pi * r0 * r0);
    }

    static double func_double_quad(double eta, void *params)
    {
        struct f_params *value = (struct f_params *)params;
        struct f_params2 params2;
        params2.r = value->r;
        params2.x3 = value->x3;
        params2.eta = eta;

        double r0 = value->r0;
        params2.r0 = r0;

        Function<double, void> F(func_double_quad2, &params2);
        double lim = sqrt(r0 * r0 - eta * eta);
        return quadrature(F, -lim, lim, 'c');
    }

    static double func_double_quad2(double zeta, void *params)
    {
        struct f_params2 *value = (struct f_params2 *)params;
        double r = value->r, x3 = value->x3, eta = value->eta, r0 = value->r0;

        return (2 *
                ellint_1(sqrt(4 * r * (1 + eta) /
                              (pow((1 + eta + r), 2) + pow((x3 - zeta), 2)))) *
                sqrt(pow((1 + eta), 2) /
                     (pow((1 + eta + r), 2) + pow((x3 - zeta), 2))));
    }

    static double quadrature(Function<double, void> F, double low, double top,
                             char method = 'b')
    {
        double result, abserr, epsabs = 1e-10, epsrel = 0;
        switch (method)
        {
        case 'a':
        {
            // do the quadrature with quadpackpp.
            size_t limit = 128, m_deg = 10;
            Workspace<double> Work(limit, m_deg);
            int status = 0;
            try
            {
                status = Work.qag(F, low, top, epsabs, epsrel, result, abserr);
            }
            catch (const char *reason)
            {
                cerr << reason << endl;
                return status;
            }
            break;
        }
        case 'b':
        {
            // do the quadrature with gsl_integration.
            gsl_integration_workspace *w =
                gsl_integration_workspace_alloc(5000);
            gsl_function Func;
            Func.function = F.function_;
            Func.params = F.params_;

            struct f_params *value = (struct f_params *)Func.params;
            double r0 = value->r0, x3 = value->x3;

            if (abs(x3) <= r0)
            {
                double singular = asin(x3 / r0);
                double pts[] = {low, singular, top};
                size_t npts = 3;
                gsl_integration_qagp(&Func, pts, npts, epsabs, epsrel, 5000, w,
                                     &result, &abserr);
            }
            else
            {
                gsl_integration_qags(&Func, low, top, epsabs, epsrel, 5000, w,
                                     &result, &abserr);
            }
            gsl_integration_workspace_free(w);
            break;
        }
        case 'c':
        {
            // do the quadrature with gsl_integration but it's for double
            // integral.
            epsabs = 1e-8;
            gsl_integration_workspace *w =
                gsl_integration_workspace_alloc(5000);
            gsl_function Func;
            Func.function = F.function_;
            Func.params = F.params_;
            gsl_integration_qags(&Func, low, top, epsabs, epsrel, 5000, w,
                                 &result, &abserr);
            gsl_integration_workspace_free(w);
            break;
        }
        }
        // cout << "epsabs: " << epsabs << endl;
        // cout << setprecision(12) << "result: " << result << endl
        //      << "abserr: " << abserr << endl;
        return result;
        // return gauss_legendre(128, func_single_quad, &params, 0, 2 * pi);
    }

    double potential_func3() { return 0; }
};

class Particle
{
  public:
    Particle() = default;
    ~Particle() = default;
    Particle(pos position) : r(position[0]), z(position[1]) {}
    Particle(state_type5 state)
        : r(state[0]), z(state[1]), r_dot(state[2]), z_dot(state[3]),
          lam_dot(state[4]) {}
    Particle(state_type5 state, double time)
        : r(state[0]), z(state[1]), r_dot(state[2]), z_dot(state[3]),
          lam_dot(state[4]), t(time) {}

    double r = 0.0, z = 0.0, r_dot = 0.0, z_dot = 0.0, lam_dot = 0.0, t = 0.0;

    state_type5 convert2state() const { return {r, z, r_dot, z_dot, lam_dot}; }
    pos convert2pos() const { return {r, z}; }
    vel convert2vel() const { return {r_dot, z_dot, lam_dot}; }
};

int main()
{
    Torus torus;

    Particle par1({1.25, 0.2, 3, 4, 5});
    // double temp;
    // temp = torus.potential(par1.r, par1.z, 'a');
    // cout << setprecision(15) << temp << endl;
    // temp = torus.potential(par1.r, par1.z, 'b');
    // cout << setprecision(15) << temp << endl;
    // double m = 0.94;

    // container rr = ratio_p(1.5, 3.5, 20);
    // for (auto ele : rr)
    // {
    // pos pp = torus.r2u({1.25, 0.0});

    // cout << pp[0] << endl
    //  << pp[1] << endl;

    torus.sampling();
    torus.dct();

    double v1, v2;
    const clock_t start = clock();
    // for (size_t i = 0; i < 1000; i++)
    // {
    v1 = torus.potential_harmonics_rz(par1.r, par1.z, 60);
    v2 = torus.potential(par1.r, par1.z);
    // }

    cout << endl
         << "Cpu Time: "
         << static_cast<double>(clock() - start) / CLOCKS_PER_SEC << endl;
    cout << setprecision(15) << v1 << endl
         << v2 << endl;
}
