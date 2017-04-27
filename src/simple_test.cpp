#include <fftw3.h>
#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

void dump_vector(int n, double *vec) {
    for (int i = 0; i < n; i++)
        cout << setprecision(15) << vec[i] << '\t';
    printf("\n");
}
int main() {
    int n = 6;
    std::vector<double> p = {1, 1, 1, 2, 5,6};
    
    for (auto ele : p)
    {
        ele /= n;
    };
    double *a = &p[0];
    // double a[] = {1, 1, 1, 2,5};
    // double b[] = {0, 0, 0, 0, 5};
    printf("Original vector\n");
    dump_vector(n, a);
    fftw_plan plan = fftw_plan_r2r_1d(n, a, a, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_execute(plan);
    printf("DCT\n");
    dump_vector(n, a);
    fftw_plan plani = fftw_plan_r2r_1d(n, a, a, FFTW_REDFT01, FFTW_ESTIMATE);
    fftw_execute(plani);
    printf("IDCT\n");
    dump_vector(n, a);
    return 0;
}