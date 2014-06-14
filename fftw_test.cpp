#include <iostream>
#include <fftw3.h>
#define N 4

using namespace std;

int main(void)
{
    double *in;
    fftw_complex *out;
    fftw_plan p, ip;

    in = (double*) fftw_malloc(sizeof(double) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
//    p = fftw_plan_r2r_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
    ip = fftw_plan_dft_c2r_1d(N, out, in, FFTW_ESTIMATE);

    in[0] = 1.0;
    in[1] = 0;
    in[2] = 0;
    in[3] = 0;

    fftw_execute(p); /* repeat as needed */
    for (int i = 0; i < N; ++i)
        cout << out[i][0] << " + j" << out[i][1] << endl;

    fftw_execute(ip);
    for (int i = 0; i < N; ++i)
        cout << in[i]/N << " ";
    cout << endl;

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    return 0;
}
