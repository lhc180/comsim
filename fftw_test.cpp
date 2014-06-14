#include <iostream>
#include <fftw3.h>
#define N 4

using namespace std;

int main(void)
{
    double *in, *in2, *in_m;
    fftw_complex *out, *out2, *out_m;
    fftw_plan p, p2, ip;

    in = (double*) fftw_malloc(sizeof(double) * N);
    in2 = (double*) fftw_malloc(sizeof(double) * N);
    in_m = (double*) fftw_malloc(sizeof(double) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out_m = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//    p = fftw_plan_r2r_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    p2 = fftw_plan_dft_r2c_1d(N, in2, out2, FFTW_ESTIMATE);
//    ip = fftw_plan_dft_c2r_1d(N, out, in, FFTW_ESTIMATE);
    ip = fftw_plan_dft_c2r_1d(N, out_m, in_m, FFTW_ESTIMATE);

    in[0] = 1.0;
    in[1] = 0;
    in[2] = 1.0;
    in[3] = 0;

    in2[0] = 1.0;
    in2[1] = 1.0;
    in2[2] = 0;
    in2[3] = 0;

    fftw_execute(p); /* repeat as needed */
    fftw_execute(p2);
//    for (int i = 0; i < N; ++i)
//        cout << out[i][0] << " + j" << out[i][1] << endl;

    for (int i = 0; i < N/2 + 1; ++i) {
        out_m[i][0] = (out[i][0]*out2[i][0] - out[i][1]*out2[i][1])/N;
        out_m[i][1] = (out[i][0]*out2[i][1] + out[i][1]*out2[i][0])/N;
    }
    fftw_execute(ip);
    for (int i = 0; i < N; ++i)
        cout << in_m[i] << " ";
    cout << endl;

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(in2);
    fftw_free(in_m);
    fftw_free(out);
    fftw_free(out2);
    fftw_free(out_m);
    return 0;
}
