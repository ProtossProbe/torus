//
//  torus_run.cpp
//
//
//  Created by Protoss Probe on 2016/10/13.
//  Copyright © 2016年 probe. All rights reserved.
//

#include <iostream>
#include <libiomp/omp.h>
#include <math.h>
#include <time.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/heuman_lambda.hpp>
#include <gsl/gsl_integration.h>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/tools/roots.hpp>
#include "gauss_legendre.h"
// #include "torus_gravity.hpp"

using namespace std;
using namespace boost::math;
using namespace boost::numeric::odeint;

typedef boost::array<double, 2> state_type2;
typedef boost::array<double, 5> state_type5;
typedef boost::array<double, 6> state_type6;
typedef boost::array<double, 21> state_type21;
typedef boost::array<double, 23> state_type23;

typedef runge_kutta_dopri5<double, double, double, double, vector_space_algebra, default_operations, never_resizer> stepper_type;

const double r0 = 0.25; // 定义 torus 中小圆的半径 r0，默认大圆半径 R0=1
const double lam = 0;
//const double lam_value1 = 1.53765218258;
//const double lam_value2 = 1.69824839234;
const double pi = 3.141592653589793238463;
const string filename = "/Users/huangyukun/Desktop/test/test/torus_data.txt";
const string dotfile = "/Users/huangyukun/Desktop/test/test/dot_data.txt";
const string lamfile = "/Users/huangyukun/Desktop/test/test/lam_data.txt";
const string initfile = "/Users/huangyukun/Desktop/test/test/init_data.txt";
const string orbitfile = "/Users/huangyukun/Desktop/test/test/orbit_data.txt";
const string r_point_file = "/Users/huangyukun/Desktop/test/test/r_point_data2.txt";
ofstream outputFile;
int ind = 0;

int in_out(double r, double x3);

double func1(double theta, void *params);
double func2(double theta, void *params);
double func3(double theta, void *params);
double func4(double theta, void *params);
double func5(double theta, void *params);
double func6(double theta, void *params);
double func7(double theta, void *params);

double I_000(double a, double r, double z);
double I_010(double a, double r, double z);
double I_001(double a, double r, double z);
double I_011(double a, double r, double z);

double potential_func(const double &r, const double &x3);
double potential_func2(const double &r, const double &x3);
double hamilton(const state_type5 &x);
double solve_z_dot(double r, double z, double rdot, double hamilton);
double u_p_circle(double eta, double zeta, void *params);
double u_p_int1(double eta, void *params);

double gz_f(const double &r, const double &x3);
double gr_f(const double &r, const double &x3);
double u_rr(const double &r, const double &x3);
double u_zz(const double &r, const double &x3);
double u_rz(const double &r, const double &x3);
void grav_f(const double &r, const double &x3, state_type2 &f, double step_size = 1e-4);

struct f_params
{
    double r;
    double x3;
};
struct f_params2
{
    double r;
    double x3;
    double eta;
};
double r_force(const double &r, const double &lam);

double search_lam0(double lam0_init, double lam_up, double lam_low);
double search_lam1(double lam0_init, double lam_up, double lam_low);
double search_r_eq(double lam1, double lam2);

void run(int i, double r, double r_dot, double hamilton, string dotfile);
void run_correction(int i, double r, double r_dot, double hamilton, string initfile);

class torus_ode
{
    double lam;

  public:
    torus_ode(double gam) : lam(gam) {}
    void operator()(const state_type5 &x, state_type5 &dxdt, double t)
    {
        state_type2 force;
        //        grav_f(x[0],x[1],force);
        force[0] = gr_f(x[0], x[1]);
        force[1] = gz_f(x[0], x[1]);

        if (abs(x[0]) < 1e-16)
        {
            dxdt[0] = x[2];
            dxdt[1] = x[3];
            dxdt[2] = force[0];
            dxdt[3] = force[1];
            dxdt[4] = 0;
        }
        else
        {

            dxdt[0] = x[2];
            dxdt[1] = x[3];
            dxdt[2] = force[0] + lam * lam / pow(x[0], 3);
            dxdt[3] = force[1];
            dxdt[4] = lam / (x[0] * x[0]);
        }
    }
};

class torus_ode21
{
    double lam;

  public:
    torus_ode21(double gam) : lam(gam) {}
    void operator()(const state_type21 &x, state_type21 &dxdt, double t)
    {
        state_type2 force;
        double urr;
        double uzz;
        double urz;
#pragma omp parallel sections num_threads(5)
        {
#pragma omp section
            {
                force[0] = gr_f(x[0], x[1]);
            }
#pragma omp section
            {
                force[1] = gz_f(x[0], x[1]);
            }
#pragma omp section
            {
                urr = u_rr(x[0], x[1]) + 3 * lam * lam / pow(x[0], 4);
            }
#pragma omp section
            {
                uzz = u_zz(x[0], x[1]);
            }
#pragma omp section
            {
                urz = u_rz(x[0], x[1]);
            }
        }

        if (abs(x[0]) < 1e-16)
        {
            dxdt[0] = x[2];
            dxdt[1] = x[3];
            dxdt[2] = force[0];
            dxdt[3] = force[1];
            dxdt[4] = 0;
        }
        else
        {
            dxdt[0] = x[2];
            dxdt[1] = x[3];
            dxdt[2] = force[0] + lam * lam / pow(x[0], 3);
            dxdt[3] = force[1];
            dxdt[4] = lam / (x[0] * x[0]);
        }

        for (int i = 5; i <= 12; i++)
        {
            dxdt[i] = x[i + 8];
        }
        for (int i = 13; i <= 16; i++)
        {
            dxdt[i] = -x[i - 8] * urr - x[i - 4] * urz;
        }

        for (int i = 17; i <= 20; i++)
        {
            dxdt[i] = -x[i - 12] * urz - x[i - 8] * uzz;
        }
    }
};

void write_cout(const state_type5 &x, const double t)
{
    //    double E;
    //    E = hamilton(x);
    //    outputFile << setprecision(12) << t << '\t' << x[0] << '\t' << x[1] <<
    //    '\t' << x[2] << '\t' << x[3] << '\t'<< x[4] <<'\t'<< E << '\t'  << endl;
    //    cout << t << '\t'<< x[0] << '\t' << x[1] << '\t' << E << '\t'  << endl;

    outputFile << ind << '\t';
    outputFile << setprecision(12) << t << '\t';

    for (int i = 0; i < 4; i++)
    {
        outputFile << x[i] << '\t';
    }
    double E;
    state_type5 x_temp;
    for (int i = 0; i < 5; i++)
    {
        x_temp[i] = x[i];
    }
    E = hamilton(x_temp);

    outputFile << E << '\t' << endl;

    cout << ind << '\t' << t << '\t' << x[0] << '\t' << x[1] << '\t' << endl;
}

void read_file_and_plot(state_type6 array, string initfile)
{
    ifstream file;
    string line;
    file.open(initfile);
    while (getline(file, line))
    {
        stringstream ss(line);
        for (int i = 0; i < 5; i++)
        {
            ss >> array[i];
        }

        outputFile.open(orbitfile, ios::app);
        state_type5 x;
        for (int i = 0; i < 5; i++)
        {
            x[i] = array[i + 1];
        }

        runge_kutta_dopri5<state_type5> stepper;
        integrate_const(stepper, torus_ode(lam), x, 0.0, array[0] + 0.01, 0.01, write_cout);
        ind++;
        outputFile.close();
    }
}

int main(int argc, const char *argv[])
{
    const clock_t begin_time = omp_get_wtime();
    const clock_t start = clock();

    // state_type6 data;
    // read_file_and_plot(data,initfile);

    //
    //    cout << I_000(0.25, 0.1, 0) << endl;
    //    cout << I_000(0.25, 0.1, 0.1) << endl;
    //    cout << I_000(0.25, 0.1, -0.1) << endl;

    //    double r_dot_range = -1.5;
    //    double r1[32] = {0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65, 0.70, 0.73, 1.26, 1.30,1.35,1.40,1.45,1.50,1.55,1.60,1.65,
    //        1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25};
    //    double E = -0.48;
    //    int num_i = 2;
    //    int j;
    //    double r_dot;
    //
    //    #pragma omp parallel for private(j,r_dot) num_threads(8)
    //    for (int i=0;i<num_i;i++) {
    //        int id = omp_get_thread_num();
    //        for (j=0; j<31; j++){
    //            r_dot = r_dot_range + j * 0.1;
    //            run(i*100+j, r1[i], r_dot, E, dotfile);
    //            }
    //        cout <<"Threads: " << id << '\t' << i << '\t' << "Done!" << endl;
    //    }
    //

    //    outputFile.open(filename);
    //    state_type5 x;
    //    double zdot;
    //    double r = 0.5;
    //    double r_dot = 0.0;
    //
    //    zdot = -solve_z_dot(r, 0,r_dot, -0.35);
    //    cout << potential_func2(r,0) << endl;
    //
    //    cout << zdot << endl;
    //
    //    x = {r, 0 , r_dot , zdot, 0};
    //
    //    runge_kutta_dopri5< state_type5 > stepper;
    //    integrate_const(stepper, torus_ode{lam} , x , 0.0 , 100.0 , 0.01 , write_cout);

    //    outputFile.open(initfile,ios::app);
    //    state_type21 x;
    //    double zdot;
    //    double r = 1.90;
    //    double r_dot = 0;

    ////    zdot = -solve_z_dot(r, 0,r_dot, -0.48);
    ////    cout << potential_func2(r,0) << endl;
    //
    //    zdot = -0.424850328945;
    //
    //    x = {r,0,r_dot,zdot,0, 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    ////    x = {r+0.0001,0,r_dot,zdot,0, 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
    //
    //    runge_kutta_dopri5< state_type21 > stepper;
    //    integrate_const(stepper, torus_ode21{lam} , x , 0.0 , 30.0 , 0.01 , write_cout);

    //    run_correction(0, 2.12, 0.00, -0.47, initfile);

    //    run_correction(0, 1.9, 0.0, -0.48, initfile);

    //    上面的这些代码是算torus下质点运行的轨道的

    //    for (int i=1; i<=12; i++) {
    //        r0 = r0 + 0.001;
    //        ofstream lamdata(lamfile,ios::app);
    //        lamdata << r0 << '\t';
    //        lamdata << setprecision(12) << search_lam0(1.53,1.50,1.55) << '\t' << search_lam1(2,1.0,3)<< endl;
    //        lamdata.close();
    //    }

    //    double interv = lam_value2 - lam_value1;
    //    double num = 1000;
    //    double step = interv / num /1000;
    //    ofstream r_point_data(r_point_file,ios::app);
    //    for (int i=1; i<=num; i++) {
    //        double lam = lam_value2 + 0.001 * i;
    //        double r_low = 1.61812;
    //        double r_up = 1+r0;
    //        double r = 1.5;

    //        for (int j=1; j<=100; j++) {
    //            if (r_force(r,lam) > 1e-8) {
    //                r_up = r;
    //                r = (r + r_low)/2.;
    //                //cout << r << '\t' << "step forward" << endl;
    //            }else if (r_force(r, lam) < -1e-8){
    //                r_low = r;
    //                r = (r + r_up)/2.;
    //                //cout << r << '\t' << "step back" << endl;
    //            }
    //            else {
    //                cout << setprecision(15)<< r << '\t' << "end~!" <<endl;
    //                break;
    //
    //                }
    //            }
    //        r_point_data << setprecision(15) << lam << '\t' <<  r << '\t';

    //        r_low = 2.5;
    //        r_up = 8;
    //        r = 2.6;
    //
    //        for (int j=1; j<=100; j++) {
    //            if (r_force(r,lam) > 1e-8) {
    //                r_up = r;
    //                r = (r + r_low)/2.;
    //                //cout << r << '\t' << "step forward" << endl;
    //            }else if (r_force(r, lam) < -1e-8){
    //                r_low = r;
    //                r = (r + r_up)/2.;
    //                //cout << r << '\t' << "step back" << endl;
    //            }
    //            else {
    //                cout << setprecision(15)<< r << '\t' << "end~!" <<endl;
    //                break;
    //
    //            }
    //
    //        }
    //        r_point_data << setprecision(15) << lam << '\t' << r << '\t' << endl;
    //
    //    }

    //  上面这些代码是计算角动量变化带来的周期轨道分叉的

    cout << setprecision(12) << endl
         << potential_func(0.7, 0.01) << endl;
    cout << setprecision(12) << endl
         << potential_func2(0.7, 0.01) << endl;

    //    cout <<endl<<"Time:"<< float( omp_get_wtime() - begin_time ) << endl;
    //    cout << endl << "Cpu Time: " << static_cast<double>(clock() - start)/CLOCKS_PER_SEC << endl;
}

void run(int i, double r, double r_dot, double hamilton, string dotfile)
{
    ofstream dotFile;
    dotFile.open(dotfile, ios::app);

    double tmax = 500.0; // 最大积分时间
    double step = 0.01;  // 平常积分采用的步长
    double dt = step;    // 实际积分采用的步长
    double t0 = 0.0;     // 初始时刻
    double current_t;    // 积分的当前时刻
    int cross_count;     // 正向（zdot<0）穿越赤道面的次数
    int max_count = 25;  // 最大穿越次数
    double zdot;
    double r_1 = r;

    double E = hamilton;

    //        outputFile.open(filename);
    dt = step;
    current_t = 0.0;
    cross_count = 0;

    zdot = -solve_z_dot(r_1, 0, r_dot, E);
    if (zdot == 0.0)
        return;
    state_type5 x = {{r_1, 0, r_dot, zdot, 0}};
    state_type5 x_last(x);

    dotFile << i << '\t' << setprecision(12) << current_t << '\t' << x[0] << '\t' << x[2] << '\t' << x[3] << '\t' << endl;

    //    double lam = x[0] * x[0] * lam0; // Lambda = r ^ 2 * lambda_dot0

    runge_kutta_dopri5<state_type5> stepper;
    while ((in_out(x[0], x[1]) == 1) and (current_t <= tmax) and cross_count < max_count)
    {
        if ((x_last[1] * x[1] >= 0.) or x_last[1] < 0.)
        {
            //                write_cout(x , current_t);
            x_last = x;
            stepper.do_step(torus_ode(lam), x, t0, dt);
            current_t = current_t + dt;
        }
        else if (abs(x_last[1]) > 1e-7)
        {
            current_t = current_t - dt;
            dt = dt / 10.;
            x = x_last;
            stepper.do_step(torus_ode(lam), x, t0, dt);
            current_t = current_t + dt;
        }
        else if (abs(x_last[1]) < 1e-7)
        {
            dt = step;
            x_last = x;
            cross_count++;
            dotFile << i << '\t' << setprecision(12) << current_t << '\t' << x[0] << '\t' << x[2] << '\t' << x[3] << '\t' << endl;
            cout << i << '\t' << setprecision(12) << current_t << '\t' << x[0] << '\t' << x[2] << '\t' << x[3] << '\t' << endl;
        }
    }
}

void run_correction(int i, double r, double r_dot, double hamilton, string initfile)
{
    ofstream initFile;
    initFile.open(initfile, ios::app);

    double tmax = 500.0; // 最大积分时间
    double step = 0.01;  // 平常积分采用的步长
    double dt = step;    // 实际积分采用的步长
    double t0 = 0.0;     // 初始时刻
    double current_t;    // 积分的当前时刻
    int cross_count;     // 正向（zdot<0）穿越赤道面的次数
    int max_count = 2;   // 最大穿越次数
    double zdot;
    state_type21 x, x_last, x0;

    double r_1 = r;
    double E = hamilton;

    //    outputFile.open(filename);
    dt = step;
    current_t = 0.0;
    cross_count = 0;

    zdot = -solve_z_dot(r_1, 0, r_dot, E);
    if (zdot == 0.0)
        return;
    x0 = {r, 0, r_dot, zdot, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1}; // initial conditions
    x = x0;
    x_last = x0;

    //    initFile << i << '\t' << setprecision(12) << current_t << '\t' << x[0] << '\t' << x[2] << '\t' << x[3] << '\t' << endl;

    //    double lam = x[0] * x[0] * lam0; // Lambda = r ^ 2 * lambda_dot0

    runge_kutta_dopri5<state_type21> stepper;
    while ((in_out(x[0], x[1]) == 1) and (current_t <= tmax) and cross_count < max_count)
    {
        if ((x_last[1] * x[1] >= 0.) or x_last[1] < 0.)
        {
            //            write_cout(x , current_t);
            x_last = x;
            stepper.do_step(torus_ode21(lam), x, t0, dt);
            current_t = current_t + dt;
            if (in_out(x[0], x[1]) != 1)
                cout << "Crash!" << endl;
        }
        else if (abs(x_last[1]) > 1e-7)
        {
            current_t = current_t - dt;
            dt = dt / 10.;
            x = x_last;
            stepper.do_step(torus_ode21(lam), x, t0, dt);
            current_t = current_t + dt;
        }
        else if (abs(x_last[1]) < 1e-7)
        {
            cross_count++;
            double delta_r = x[0] - r;
            double delta_r_dot = x[2] - r_dot;
            double dist;
            dist = sqrt(pow(delta_r, 2) + pow(delta_r_dot, 2));
            if (cross_count == max_count)
            {
                if (dist > 1e-6)
                {
                    cross_count = 0;
                    double r_double_dot = gr_f(x[0], x[1]) + lam * lam / pow(x[0], 3);

                    double a = x[7] - x[11] * x[2] / x[3];
                    double b = x[8] - x[12] * x[2] / x[3];
                    double c = x[15] - x[11] * r_double_dot / x[3];
                    double d = x[16] - x[12] * r_double_dot / x[3];
                    double corr_r_dot = (b * delta_r_dot - d * delta_r) / (a * d - b * c);
                    double corr_z_dot = (-a * delta_r_dot + c * delta_r) / (a * d - b * c);
                    x0[2] = x0[2] + corr_r_dot;
                    x0[3] = x0[3] + corr_z_dot;

                    dt = step;
                    x = x0;
                    x_last = x;

                    //                    cout << setprecision(12)<< a << '\t' << b << '\t' << c << '\t' << d << '\t' << endl;
                    cout << setprecision(12) << "Error is: " << delta_r << '\t' << delta_r_dot << endl;
                    cout << setprecision(12) << "Correction is: " << corr_r_dot << '\t' << corr_z_dot << '\t' << endl;
                    cout << i << '\t' << setprecision(12) << current_t << '\t' << x[0] << '\t' << x[2] << '\t' << x[3] << '\t' << endl;
                    cout << endl;
                    current_t = 0.0;
                }
                else if (dist <= 1e-7)
                {
                    //                    initFile << i << '\t'<< setprecision(12) << current_t << '\t' << x[0] << '\t' << x[2] << '\t' << x[3] << '\t' << endl;
                    //                    write_cout(x, current_t);
                    cout << "Succeed!" << endl;
                    cout << setprecision(12) << "Error is: " << delta_r << '\t' << delta_r_dot << endl;
                    cout << i << '\t' << setprecision(12) << current_t << '\t' << x[0] << '\t' << x[2] << '\t' << x[3] << '\t' << endl;
                    break;
                }
            }
            else
            {
                dt = step;
                x_last = x;
            }
        }
    }
}

int in_out(double r, double x3)
{
    r = abs(r);
    if ((1 - r) * (1 - r) + x3 * x3 < r0 * r0)
    {
        return 0;
    }
    else if ((1 - r) * (1 - r) + x3 * x3 > r0 * r0)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

double func1(double theta, void *params)
{
    // 用于积分的一个关于 theta 的函数，将 r 和 x3 作为参量输入

    struct f_params *value = (struct f_params *)params;
    double r = (value->r);
    double x3 = (value->x3);
    r = abs(r);

    double R1, a, b, c, nu, k_prime, k1, result, temp;
    R1 = 1 + r0 * cos(theta);
    a = 2 * (pow(r, 2) + pow((x3 - r0 * sin(theta)), 2));
    b = a / 2. - pow(R1, 2) + sqrt(pow((a / 2 - pow(R1, 2)), 2) + 4 * pow(R1, 2) * pow((x3 - r0 * sin(theta)), 2));
    c = a / 2. - pow(R1, 2) - sqrt(pow((a / 2 - pow(R1, 2)), 2) + 4 * pow(R1, 2) * pow((x3 - r0 * sin(theta)), 2));
    nu = (a - b) / (2 * pow(r, 2));
    k_prime = sqrt((pow((r - R1), 2) + pow((x3 - r0 * sin(theta)), 2)) / (pow((r + R1), 2) + pow((x3 - r0 * sin(theta)), 2)));
    k1 = (1 - k_prime) / (1 + k_prime);

    if (1 - nu < 1e-12)
    {
        temp = 1e8;
    }
    else
    {
        temp = ellint_3(k1, nu);
    }

    result = ((c + 2 * (pow(R1, 2) - pow(r, 2))) * ellint_1(k1) + (a - c) * ellint_2(k1) - 2 * pow((x3 - r0 * sin(theta)), 2) * temp) * cos(theta) / sqrt(a - c);
    //    cout << setprecision(12) << result << endl << nu <<endl<< k1 << endl << theta <<endl;
    return result * sqrt(2) / (pi * pi * r0);
}

double func2(double theta, void *params)
{
    // 用于积分的一个关于 theta 的函数，将 r 和 x3 作为参量输入

    struct f_params *value = (struct f_params *)params;
    double r = (value->r);
    double x3 = (value->x3);
    r = abs(r);

    double R1, p, q, k, k_p, z, phi, result;
    int sgn;
    R1 = 1 + r0 * cos(theta);
    z = x3 - r0 * sin(theta);

    p = sqrt(pow((R1 + r), 2) + z * z);
    q = sqrt(pow((R1 - r), 2) + z * z);
    k = 2 / p * sqrt(r * R1);
    k_p = sqrt(1 - k * k);
    phi = asin(abs(z) / q);

    if (R1 - r > 0)
    {
        sgn = 1;
    }
    else if (R1 - r < 0)
    {
        sgn = -1;
    }

    result = (-p * ellint_2(k) - (R1 * R1 - r * r) / p * ellint_1(k) + abs(z) * (pi / 2 + pi / 2 * sgn) - abs(z) * sgn * (ellint_2(k) * ellint_1(k_p, phi) + ellint_1(k) * ellint_2(k_p, phi) - ellint_1(k) * ellint_1(k_p, phi)));

    return -result * cos(theta) / (pi * pi * r0);
}

double func3(double theta, void *params)
{
    // 用于积分的一个关于 theta 的函数，将 r 和 x3 作为参量输入

    struct f_params *value = (struct f_params *)params;
    double r = (value->r);
    double x3 = (value->x3);

    double R1, p, q, k, k_p, z, phi, result;
    int sgn = 0;

    z = r * sin(theta);
    r = sqrt(pow(1 - r * cos(theta), 2) + x3 * x3);
    R1 = r0;

    p = sqrt(pow((R1 + r), 2) + z * z);
    q = sqrt(pow((R1 - r), 2) + z * z);
    k = 2 / p * sqrt(r * R1);
    k_p = sqrt(1 - k * k);
    phi = asin(abs(z) / q);

    if (R1 - r > 0)
    {
        sgn = 1;
    }
    else if (R1 - r < 0)
    {
        sgn = -1;
    }

    result = (-p * ellint_2(k) - (R1 * R1 - r * r) / p * ellint_1(k) + abs(z) * (pi / 2 + pi / 2 * sgn) - abs(z) * sgn * (ellint_2(k) * ellint_1(k_p, phi) + ellint_1(k) * ellint_2(k_p, phi) - ellint_1(k) * ellint_1(k_p, phi)));

    return -result / (pi * pi * r0 * r0);
}

double func4(double theta, void *params)
{
    struct f_params *value = (struct f_params *)params;
    double r = (value->r);
    double x3 = (value->x3);

    double a, gz;

    a = 1 + r0 * cos(theta);
    x3 = x3 - r0 * sin(theta);

    if (x3 > 0.0)
    {
        gz = 1 / (pi * r0) * a * I_000(a, r, x3) * sin(theta);
    }
    else if (x3 < 0.0)
    {
        x3 = abs(x3);
        gz = 1 / (pi * r0) * a * I_000(a, r, x3) * sin(theta);
    }

    return gz;
}

double func5(double theta, void *params)
{
    struct f_params *value = (struct f_params *)params;
    double r = (value->r);
    double x3 = (value->x3);

    double a, gr;
    a = 1 + r0 * cos(theta);
    x3 = x3 - r0 * sin(theta);

    if (x3 > 0.0)
    {
        gr = 1 / (pi * r0) * a * I_010(a, r, x3) * sin(theta);
    }
    else if (x3 < 0.0)
    {
        double x3_u;
        x3_u = abs(x3);
        gr = 1 / (pi * r0) * a * (-I_010(a, r, x3_u) + 2 * I_010(a, r, 0.)) * sin(theta);
    }
    return gr;
}

double func6(double theta, void *params)
{
    struct f_params *value = (struct f_params *)params;
    double r = (value->r);
    double x3 = (value->x3);

    double a, u_zz;
    a = 1 + r0 * cos(theta);
    x3 = x3 - r0 * sin(theta);

    if (x3 > 0.0)
    {
        u_zz = -1 / (pi * r0) * a * I_001(a, r, x3) * sin(theta);
    }
    else if (x3 < 0.0)
    {
        x3 = abs(x3);
        u_zz = 1 / (pi * r0) * a * I_001(a, r, x3) * sin(theta);
    }

    return u_zz;
}

double func7(double theta, void *params)
{
    struct f_params *value = (struct f_params *)params;
    double r = (value->r);
    double x3 = (value->x3);

    double a, u_rz;
    a = 1 + r0 * cos(theta);
    x3 = x3 - r0 * sin(theta);

    if (x3 > 0.0)
    {
        u_rz = -1 / (pi * r0) * a * I_011(a, r, x3) * sin(theta);
    }
    else if (x3 < 0.0)
    {
        double x3_u;
        x3_u = abs(x3);
        u_rz = -1 / (pi * r0) * a * I_011(a, r, x3_u) * sin(theta);
    }
    return u_rz;
}

double I_000(double a, double r, double z)
{
    double k, result;
    k = sqrt((4 * a * r) / ((a + r) * (a + r) + z * z));
    if (r == 0.)
    {
        return 0.0;
    }
    else
    {
        result = k / (pi * sqrt(a * r)) * ellint_1(k);
        return result;
    }
}

double I_010(double a, double r, double z)
{
    double k, result, c, alpha, beta;
    k = sqrt((4 * a * r) / ((a + r) * (a + r) + z * z));
    c = z * z / ((r - a) * (r - a) + z * z);

    alpha = asin(k);
    beta = asin(sqrt(c));

    if (r == 0.)
    {
        return 0.0;
    }

    if (k == 1.)
    {
        return 1 / (2 * r);
    }

    else if (r - a >= 1e-6)
    {
        result = -k * z / (4 * r * sqrt(a * r)) * ellint_1(k) * 2 / pi - 1 / (2 * r) * heuman_lambda(k, beta) + 1 / r;
    }
    else if (r - a < 1e-6 && r - a >= -1e-6)
    {
        result = -k * z / (4 * r * r) * ellint_1(k) * 2 / pi + 1 / (2 * r);
    }
    else if (r - a < -1e-6)
    {
        result = -k * z / (4 * r * sqrt(a * r)) * ellint_1(k) * 2 / pi + 1 / (2 * r) * heuman_lambda(k, beta);
    }

    return result;
}

double I_001(double a, double r, double z)
{
    double k;
    k = sqrt((4 * a * r) / ((a + r) * (a + r) + z * z));
    if (r == 0.0)
    {
        return 0.0;
    }
    else
    {
        return z * k * k * k * ellint_2(k) / (4 * pi * (1 - k * k) * sqrt(pow((a * r), 3)));
    }
}

double I_011(double a, double r, double z)
{
    double k, k_p2;
    k = sqrt((4 * a * r) / ((a + r) * (a + r) + z * z));
    k_p2 = 1 - (4 * a * r) / ((a + r) * (a + r) + z * z);
    if (r == 0.0)
    {
        return 0.0;
    }
    else
    {
        return 1 / pi * (k * k * k * (r * r - a * a - z * z) / (8 * r * k_p2 * sqrt(pow((a * r), 3))) * ellint_2(k) + k / (2 * r * sqrt(r * a)) * ellint_1(k));
    }
}

double gz_f(const double &r, const double &x3)
{
    struct f_params params;
    params.r = abs(r);
    params.x3 = x3;
    return -gauss_legendre(128, func4, &params, 0, 2 * pi);
}

double gr_f(const double &r, const double &x3)
{
    struct f_params params;
    params.r = r;
    params.x3 = x3;
    if (r > 0.)
    {
        return -gauss_legendre(128, func5, &params, 0, 2 * pi);
    }
    else
    {
        params.r = -r;
        return gauss_legendre(128, func5, &params, 0, 2 * pi);
    }
}

double u_rr(const double &r, const double &x3)
{
    double urr;
    double uzz, ur;
    uzz = u_zz(r, x3);
    ur = -gr_f(r, x3);

    if (abs(r) < 1e-12)
    {
        urr = -uzz;
    }

    else
    {

        if (in_out(r, x3) == 1)
        {
            urr = -1 / r * ur - uzz;
        }
        else if (in_out(r, x3) == 0)
        {
            urr = -1 / r * ur - uzz;
        }
    }
    return urr;
}

double u_zz(const double &r, const double &x3)
{
    struct f_params params;
    params.r = abs(r);
    params.x3 = x3;
    return gauss_legendre(128, func6, &params, 0, 2 * pi);
}

double u_rz(const double &r, const double &x3)
{
    struct f_params params;
    double sgn = r > 0.0 ? 1 : -1;
    params.r = abs(r);
    params.x3 = x3;
    return sgn * gauss_legendre(128, func7, &params, 0, 2 * pi);
}

double potential_func(const double &r, const double &x3)
{
    // 计算 torus 的引力势能函数

    double u_p;
    if (abs(r) < 1e-12)
    {
        double lam, k;
        lam = 2 * r0 * sqrt(1 + pow(x3, 2)) / (1 + pow(r0, 2) + pow(x3, 2));
        k = sqrt(2 * lam / (1 + lam));

        cout << lam << endl
             << k << endl;

        u_p = (4.0 / 3.0 / (pi * r0) * sqrt(4 - 2 * pow(k, 2)) / pow(k, 2) * sqrt(1 + pow(r0, 2) / (1 + pow(x3, 2))) *
               (ellint_2(k) - 2 * (1 - pow(k, 2)) / (2 - pow(k, 2)) * ellint_1(k)));
    }
    else
    {

        //        gsl_integration_workspace * w
        //        = gsl_integration_workspace_alloc (1000);

        //        double result, error;
        struct f_params params;
        params.r = r;
        params.x3 = x3;

        //        gsl_function F;
        //        F.function = &func2;
        //        F.params = &params;
        //
        //        gsl_integration_qags (&F, 0, 2*pi, 0, 1e-10, 1000,
        //                              w, &result, &error);
        //        gsl_integration_workspace_free(w);
        //        return result;

        return gauss_legendre(128, func2, &params, 0, 2 * pi);
    }
    return u_p;
}

double hamilton(const state_type5 &x)
{
    if (x[0] != 0.0)
    {
        return 0.5 * (x[2] * x[2] + lam * lam / (x[0] * x[0]) + x[3] * x[3]) - potential_func2(x[0], x[1]);
    }
    else
    {
        return 0.5 * (x[2] * x[2] + x[3] * x[3]) - potential_func2(x[0], x[1]);
    }
}

double solve_z_dot(double r, double z, double rdot, double hamilton)
{
    double result = (hamilton + potential_func2(r, z)) * 2 - rdot * rdot;
    if (result >= 0)
    {
        return sqrt(result);
    }
    else
    {
        return 0.; //means error
    }
}

double u_p_circle(double zeta, void *params)
{
    struct f_params2 *value = (struct f_params2 *)params;
    double r = (value->r);
    double x3 = (value->x3);
    double eta = (value->eta);
    return (2 * ellint_1(sqrt(4 * r * (1 + eta) / (pow((1 + eta + r), 2) + pow((x3 - zeta), 2)))) * sqrt(pow((1 + eta), 2) / (pow((1 + eta + r), 2) + pow((x3 - zeta), 2))));
}

double u_p_int1(double eta, void *params)
{
    struct f_params *value = (struct f_params *)params;
    struct f_params2 params2;
    params2.r = (value->r);
    params2.x3 = (value->x3);
    params2.eta = eta;

    return gauss_legendre(256, u_p_circle, &params2, -sqrt(r0 * r0 - eta * eta), sqrt(r0 * r0 - eta * eta));
}

double potential_func2(const double &r, const double &x3)
{
    // 以线圆环作为基元进行积分 (Old)

    struct f_params params;
    params.r = abs(r);
    params.x3 = x3;
    return gauss_legendre(256, u_p_int1, &params, -r0, r0) / (pi * pi * r0 * r0);
}

void grav_f(const double &r, const double &x3, state_type2 &v, double step_size)
{
    // 根据 torus 的引力势计算任意一点的重力
    if (abs(r) < 1e-12 and abs(x3) < 1e-12)
    {
        v[0] = 0.0;
        v[1] = 0.0;
    }
    else
    {
        double v1, v2, v3, v4, v5, v6, v7, v8;
        omp_set_num_threads(8);
#pragma omp parallel sections
        {
#pragma omp section
            {
                v1 = potential_func(r + 2 * step_size, x3);
            }

#pragma omp section
            {
                v2 = potential_func(r + step_size, x3);
            }

#pragma omp section
            {
                v3 = potential_func(r - step_size, x3);
            }

#pragma omp section
            {
                v4 = potential_func(r - 2 * step_size, x3);
            }

#pragma omp section
            {
                v5 = potential_func(r, x3 + 2 * step_size);
            }

#pragma omp section
            {
                v6 = potential_func(r, x3 + step_size);
            }

#pragma omp section
            {
                v7 = potential_func(r, x3 - step_size);
            }

#pragma omp section
            {
                v8 = potential_func(r, x3 - 2 * step_size);
            }
        }
        v[0] = (-v1 + 8 * v2 - 8 * v3 + v4) / (12 * step_size);
        v[1] = (-v5 + 8 * v6 - 8 * v7 + v8) / (12 * step_size);
    }
}

double r_force(const double &r, const double &lam)
{
    state_type2 force;
    grav_f(r, 0, force);
    return -force[0] - lam * lam / pow(r, 3);
}

double search_lam0(double lam0_init, double lam_up, double lam_low)
{
    double x, number, lam, y, temp, result, step;
    lam = lam0_init;
    number = 10000;
    step = 0.1;
    for (int j = 1; j <= 100; j++)
    {
        result = 1 + r0 + 0.02;
        temp = r_force(1 + r0, lam);
        for (int i = 1; i <= number; i++)
        {
            x = 1 + r0 + i / number * 2;
            y = r_force(x, lam);
            if (y < temp)
            {
                temp = y;
            }
            else
            {
                result = temp;
                break;
            }
            if (i == number)
            {
                result = temp;
            }
        }
        if (r_force(x, lam) > 1e-8)
        {
            lam_up = lam;
            lam = (lam + lam_low) / 2.;
            cout << lam << '\t' << "step forward" << endl;
        }
        else if (r_force(x, lam) < -1e-8)
        {
            lam_low = lam;
            lam = (lam + lam_up) / 2.;
            cout << lam << '\t' << "step back" << endl;
        }
        else
        {
            cout << setprecision(15) << lam << '\t' << "end~!" << endl;
            break;
        }
    }
    return lam;
}

double search_lam1(double lam0_init, double lam_up, double lam_low)
{
    double lam, x;
    lam = lam0_init;
    x = r0 + 1;
    for (int j = 1; j <= 100; j++)
    {

        if (r_force(x, lam) > 1e-8)
        {
            lam_up = lam;
            lam = (lam + lam_low) / 2.;
            cout << lam << '\t' << "step forward" << endl;
        }
        else if (r_force(x, lam) < -1e-8)
        {
            lam_low = lam;
            lam = (lam + lam_up) / 2.;
            cout << lam << '\t' << "step back" << endl;
        }
        else
        {
            cout << setprecision(15) << lam << '\t' << "end~!" << endl;
            break;
        }
    }
    return lam;
}
