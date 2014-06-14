#pragma GCC diagnostic ignored "-fpermissive"
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <fftw3.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "main.h"

#define  TS         0.001
// 采样频率
#define  TDELTA     0.000001 
#define  NTMINT     ((int)(TS / TDELTA))
#define  F1         10000
#define  F2         40000
#define  PI         3.1415926535897932
#define  AMP        5

double freq_one[NTMINT], freq_zero[NTMINT];
double deviation = 50;
int pixel_cnt;

using namespace std;
using namespace cv;
//using namespace comsim;

static void modulate(const uchar intensity, double *txData[]);
static void signal_channel(double *txData[]);
static uchar demodulate(const double *txData[]);
static void filter(const double *in, const double *filt, double *outp);
static double gaussrand(void);

int main(void)
{
//    Mat image = imread("lena.pgm");
//    if (!image.data) {
//        cerr << "Could not open or find the image" << endl;
//        return 1;
//    }
//    namedWindow("Original image", WINDOW_AUTOSIZE);
//    imshow("Original image", image);

    Mat src(1, 1, CV_8UC1, Scalar(155));

    srand(time(NULL));

    ofstream ff1("f1.dat");
    ofstream ff2("f2.dat");
    for (int i = 0; i < NTMINT; ++i) {
        freq_one[i] = AMP*cos(2*PI*F1*i*TDELTA);
        ff1 << i*TS/NTMINT << " " << freq_one[i] << endl;
        freq_zero[i] = AMP*cos(2*PI*F2*i*TDELTA);
        ff2 << i*TS/NTMINT << " " << freq_zero[i] << endl;
    }
    ff1.close();
    ff2.close();

    double *txData[8];
    for (int r = 0; r < src.rows; ++r) {
        for (int c = 0; c < src.cols; ++c) {

            modulate(src.at<uchar>(r, c), txData);
            ofstream ftx("tx.dat");
            for (int i = 0; i < 8; ++i) {
                for (int j = 0; j < NTMINT; ++j)
                    ftx << i*TS+j*TS/NTMINT << " " << txData[i][j] << endl;
            }
            ftx.close();

            signal_channel(txData);
            ofstream ftx_pn("tx_pn.dat");
            for (int i = 0; i < 8; ++i) {
                for (int j = 0; j < NTMINT; ++j)
                    ftx_pn << i*TS+j*TS/NTMINT << " " << txData[i][j] << endl;
            }
            ftx_pn.close();

            uchar rxData = 0x00;
            rxData = demodulate(txData);

            for (int i = 0; i < 8; ++i)
                fftw_free(txData[i]);
            ++pixel_cnt;
            if (pixel_cnt == src.rows * src.cols) {
                cout << (unsigned int)rxData << endl;
                pixel_cnt = 0;
            }
        }
    }

    return 0;
} static void modulate(const uchar intensity, double *txData[])
{
    for (int i = 0; i < 8; ++i) {
        txData[i] = (double *)fftw_malloc(sizeof(double) * NTMINT * 8);
        if (intensity & (1 << i))
            memcpy(txData[i], freq_one, NTMINT*sizeof(double));
        else
            memcpy(txData[i], freq_zero, NTMINT*sizeof(double));
    }
}

static void signal_channel(double *txData[])
{
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < NTMINT; ++j)
            txData[i][j] += gaussrand();
    }
}

static uchar demodulate(const double *txData[])
{
    static double c1[8][NTMINT], c2[8][NTMINT];
    static uchar pixel = 0x00;
    for (int i = 0; i < 8; ++i) {
        // 上路F1解调
        // 带通滤波器
        filter(txData[i], bpf1, c1[i]);
        // 乘法器
        for (int j = 0; j < NTMINT; ++j)
            c1[i][j] *= freq_one[j];
        // 低通滤波器
        filter(c1[i], lpf, c1[i]);

        // 下路F2解调
        filter(txData[i], bpf2, c2[i]);
        for (int j = 0; j < NTMINT; ++j)
            c2[i][j] *= freq_zero[j];
        filter(c2[i], lpf, c2[i]);

        if (c1[i][500] > c2[i][500])
            pixel |= 1 << i;
    }

    ofstream fch_outp1("ch_outp1.dat");
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < NTMINT; ++j)
            fch_outp1 << i*TS+j*TS/NTMINT << " " << c1[i][j] << endl;
    }
    fch_outp1.close();
    ofstream fch_outp2("ch_outp2.dat");
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < NTMINT; ++j)
            fch_outp2 << i*TS+j*TS/NTMINT << " " << c2[i][j] << endl;
    }
    fch_outp2.close();

    return pixel;
}

static void filter(const double *in, const double *filt, double *outp)
{
    static fftw_complex ffx[NTMINT], fff[NTMINT], ffm[NTMINT];
    static fftw_plan px, pf, ip;

    px = fftw_plan_dft_r2c_1d(NTMINT, in, ffx, FFTW_ESTIMATE);
    pf = fftw_plan_dft_r2c_1d(NTMINT, filt, fff, FFTW_ESTIMATE);
    fftw_execute(px);
    fftw_execute(pf);
    for (int j = 0; j < NTMINT; ++j) {
        ffm[j][0] = (ffx[j][0]*fff[j][0] - ffx[j][1]*fff[j][1])/NTMINT;
        ffm[j][1] = (ffx[j][0]*fff[j][1] + ffx[j][1]*fff[j][0])/NTMINT;
    }
    ip = fftw_plan_dft_c2r_1d(NTMINT, ffm, outp, FFTW_ESTIMATE);
    fftw_execute(ip);
}

static double gaussrand(void)
{
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if(phase == 0) {
        do {
            double U1 = (double)rand() / RAND_MAX;
            double U2 = (double)rand() / RAND_MAX;

            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);

        X = V1 * sqrt(-2 * log(S) / S);
    } else {
        X = V2 * sqrt(-2 * log(S) / S);
    }
    phase = 1 - phase;

    return X * deviation;
}
