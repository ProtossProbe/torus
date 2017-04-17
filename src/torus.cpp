//
//  torus.cpp
//
//
//  Created by Protoss Probe on 2017/04/09.
//  Copyright © 2016-2017年 probe. All rights reserved.
//

#include <iostream>
#include <boost/array.hpp>
#include <string>
#include <cmath>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/heuman_lambda.hpp>
#include <boost/numeric/odeint.hpp>
#include <gsl/gsl_integration.h>
#include <quadpackpp/workspace.hpp>
#include "fukushima/elliptic_integral.hpp"
#include "torus.hpp"

using namespace std;
using namespace boost::math;
using namespace boost::numeric::odeint;
using namespace Elliptic_Integral;

class Torus
{
  public:
    Torus() = default;
    ~Torus() = default;
    Torus(double r0) : r0(r0) {}
    Torus(double r0, double R0) : r0(r0), R0(R0) {}
    const double r0 = 0.25, R0 = 1.0;

    double potential(const double r, const double x3, char method = 'a')
    {
	// Two different ways to calculate the potential of a torus.
	// method 'a' use a single quadrature and has high accuracy and efficiency.
	// method 'b' use a double quadrature and its accuracy and efficiency is low.

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
	cout
	    << "value: " << func_single_quad(asin(z / r0), &params) << endl;
    }

  private:
    double potential_func(const double r, const double x3)
    {
	// use potential_single_quad to do the quadrature.
	double u_p;
	if (abs(r) < 1e-12)
	{
	    double lam, k;
	    lam = 2 * r0 * sqrt(1 + pow(x3, 2)) / (1 + pow(r0, 2) + pow(x3, 2));
	    k = sqrt(2 * lam / (1 + lam));
	    u_p = (4.0 / 3.0 / (pi * r0) * sqrt(4 - 2 * pow(k, 2)) / pow(k, 2) * sqrt(1 + pow(r0, 2) / (1 + pow(x3, 2))) *
		   (ellint_2(k) - 2 * (1 - pow(k, 2)) / (2 - pow(k, 2)) * ellint_1(k)));
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

	// cout << setprecision(15) << "k: " << k << endl;

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

	result = (-p * ei2 - (R1 * R1 - r * r) / p * ei1 + abs(z) * (pi / 2 + pi / 2 * sgn) - abs(z) * sgn * (ei2 * ei1_phi + ei1 * ei2_phi - ei1 * ei1_phi));

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

	return (2 * ellint_1(sqrt(4 * r * (1 + eta) / (pow((1 + eta + r), 2) + pow((x3 - zeta), 2)))) * sqrt(pow((1 + eta), 2) / (pow((1 + eta + r), 2) + pow((x3 - zeta), 2))));
    }

    static double quadrature(Function<double, void> F, double low, double top, char method = 'a')
    {
	double result, abserr, epsabs = 1e-10, epsrel = 0;
	switch (method)
	{
	case 'a':
	{
		//do the quadrature with quadpackpp.
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
		//do the quadrature with gsl_integration.
	    gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);
	    gsl_function Func;
	    Func.function = F.function_;
	    Func.params = F.params_;

	    struct f_params *value = (struct f_params *)Func.params;
	    double r0 = value->r0, x3 = value->x3;
	    double singular = asin(x3 / r0);
	    double pts[] = {low, singular, top};
	    size_t npts = 3;

	    gsl_integration_qagp(&Func, pts, npts, epsabs, epsrel, 5000,
				 w, &result, &abserr);
	    gsl_integration_workspace_free(w);
	    break;
	}
	case 'c':
	{
		//do the quadrature with gsl_integration but it's for double integral.
	    epsabs = 1e-8;
	    gsl_integration_workspace *w = gsl_integration_workspace_alloc(5000);
	    gsl_function Func;
	    Func.function = F.function_;
	    Func.params = F.params_;
	    gsl_integration_qags(&Func, low, top, epsabs, epsrel, 5000,
				 w, &result, &abserr);
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
};

class Particle
{
  public:
    Particle() = default;
    ~Particle() = default;
    Particle(pos position) : r(position[0]), z(position[1]) {}
    Particle(state_type5 state) : r(state[0]), z(state[1]), r_dot(state[2]), z_dot(state[3]), lam_dot(state[4]) {}
    Particle(state_type5 state, double time) : r(state[0]), z(state[1]), r_dot(state[2]), z_dot(state[3]), lam_dot(state[4]), t(time) {}

    double r = 0.0, z = 0.0, r_dot = 0.0, z_dot = 0.0, lam_dot = 0.0, t = 0.0;

    state_type5 convert2state() const
    {
	return {r, z, r_dot, z_dot, lam_dot};
    }
    pos convert2pos() const
    {
	return {r, z};
    }
    vel convert2vel() const
    {
	return {r_dot, z_dot, lam_dot};
    }
};

int main()
{
    // Torus torus1;
    // Particle par1({1.25, 0.001, 3, 4, 5});
    // double temp = torus1.potential(par1.r, par1.z, 'a');
    // cout << setprecision(15) << temp << endl;
    // temp = torus1.potential(par1.r, par1.z, 'b');
    // cout << setprecision(15) << temp << endl;
    // double m = 0.94;
    // cout << setprecision(15) << ceik(m * m) - ellint_1(m) << endl;

	container cc = zonal_toroidal_harmonics(1.5,1024);
	cout << setprecision(20) << cc[1024][1] << endl;
}
