#include <fftw3.h>
#include <vector>
void dump_vector(int n, double* vec) {
    for(int i = 0; i < n; i++)
        printf("%f ", vec[i]);
    printf("\n");
}
int main()
{
    std::vector<double> p = {1, 1, 1, 2,5};
	double* a = &p[0];
    // double a[] = {1, 1, 1, 2,5};
    double b[] = {0, 0, 0, 0,5};
    printf("Original vector\n");
    dump_vector(5, a);
    fftw_plan plan = fftw_plan_r2r_1d(5, a, a, FFTW_REDFT10, FFTW_ESTIMATE);
    fftw_execute(plan);
    printf("DCT\n");
    dump_vector(5, a);
    fftw_plan plani = fftw_plan_r2r_1d(5, a, a, FFTW_REDFT01, FFTW_ESTIMATE);
    fftw_execute(plani);
    printf("IDCT\n");
    dump_vector(5, a);
    return 0;
}