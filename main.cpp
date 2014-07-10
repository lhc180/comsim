#pragma GCC diagnostic ignored "-fpermissive"
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <fftw3.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include "main.h"

#define  DEBUG
#define  TS         0.001
// 采样频率
#define  TDELTA     0.000001 
#define  NTMINT     ((int)(TS / TDELTA))
#define  F1         10000
#define  F2         40000
#define  PI         3.1415926535897932
#define  FILE_MODE  ofstream::app

double freq_one[NTMINT], freq_zero[NTMINT];
int amp = 5, deviation = 0;
bool start_up = true;
double npixel, nerror;
uchar txPixel;

using namespace std;
using namespace cv;

static void modulate(const uchar intensity, double *txData[]);
static void signal_channel(double *txData[]);
static uchar demodulate(const double *txData[]);
static void filter(const double *in, const double *filt, double *outp);
static double gaussrand(void);
void reset_err_pro(int, void*);

int main(void)
{
//#if 0
//    Mat src = imread("chessboard.png", CV_LOAD_IMAGE_GRAYSCALE);
    Mat src = imread("rand.png", CV_LOAD_IMAGE_GRAYSCALE);
    Mat dst = src.clone();
    if (!src.data) {
        cerr << "Could not open or find the image" << endl;
        return 1;
    }
    namedWindow("Original image", WINDOW_NORMAL);
    createTrackbar("Amplitude", "Original image", &amp, 5, reset_err_pro);
    createTrackbar("Deviation", "Original image", &deviation, 10, reset_err_pro);
    namedWindow("Image received", WINDOW_NORMAL);
    imshow("Original image", src);
//    imshow("image Received", dst);
//    cout << "src: " << endl << src << endl << endl;
//#endif

//    Mat src(3, 3, CV_8UC1, Scalar(155));
//    Mat dst = src.clone();

    double *txData[8];
    int pixel_cnt = 0;
    for (int i = 0; i < 8; ++i)
        txData[i] = (double *)fftw_malloc(sizeof(double) * NTMINT);
    for (int i = 0; i < NTMINT; ++i) {
        freq_one[i] = amp*cos(2*PI*F1*i*TDELTA);
        freq_zero[i] = amp*cos(2*PI*F2*i*TDELTA);
    }
#ifdef DEBUG
    ofstream ff1("data/f1.dat");
    ofstream ff2("data/f2.dat");
    for (int i = 0; i < NTMINT; ++i) {
        ff1 << i*TS/NTMINT << " " << freq_one[i] << endl;
        ff2 << i*TS/NTMINT << " " << freq_zero[i] << endl;
    }
    ff1.close();
    ff2.close();
#endif

    srand(time(NULL));
//    while (waitKey(100) < 0) {
    while (amp && waitKey(1) < 0) {
    for (int r = 0; r < src.rows; ++r) {
        for (int c = 0; c < src.cols; ++c) {

            txPixel = src.at<uchar>(r, c);
            modulate(src.at<uchar>(r, c), txData);
#ifdef DEBUG
            ofstream ftx("data/tx.dat", FILE_MODE);
            for (int i = 0; i < 8; ++i) {
                for (int j = 0; j < NTMINT; ++j)
                    ftx << i*TS+j*TS/NTMINT << " " << txData[i][j] << endl;
            }
            ftx.close();
#endif

            signal_channel(txData);
#ifdef DEBUG
            ofstream ftx_pn("data/tx_pn.dat");
            for (int i = 0; i < 8; ++i) {
                for (int j = 0; j < NTMINT; ++j)
                    ftx_pn << i*TS+j*TS/NTMINT << " " << txData[i][j] << endl;
            }
            ftx_pn.close();
#endif

            uchar rxData;
            rxData = demodulate(txData);
//            cout << (unsigned int)rxData << " ";
            dst.at<uchar>(r, c) = rxData;
//            cout << (unsigned int)rxData << " " << (unsigned int)dst.at<uchar>(r, c) << endl;;
//            dst.at<uchar>(r, c) = demodulate(txData);

            ++pixel_cnt;
            if (pixel_cnt == src.rows * src.cols) {
                pixel_cnt = 0;
                imshow("Image received", dst);
                cout << nerror / npixel << endl;
//                cout << endl << endl;
//                cout << "src:" << endl << src << endl << endl;
//                cout << "dst:" << endl << dst << endl << endl;
//                waitKey(0);
            }
#ifdef DEBUG
            if (pixel_cnt == 1)
            getchar();
#endif
        }
    }
    }
    for (int i = 0; i < 8; ++i)
        fftw_free(txData[i]);
    return 0;
}

static void modulate(const uchar intensity, double *txData[])
{
    for (int i = 0; i < 8; ++i) {
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
    static double c1[8][NTMINT], c2[8][NTMINT], inp[NTMINT], tmp[NTMINT];
//    static double c1[NTMINT], c2[NTMINT], inp[NTMINT], tmp[NTMINT];
    uchar pixel = 0x00;
    for (int i = 0; i < 8; ++i) {
        // 上路F1解调
        // 带通滤波器
        memcpy(inp, txData[i], NTMINT*sizeof(double));
        filter(inp, bpf1, c1[i]);
        // 乘法器
        for (int j = 0; j < NTMINT; ++j)
            c1[i][j] *= freq_one[j];
        // 低通滤波器
        memcpy(inp, c1[i], NTMINT*sizeof(double));
        filter(inp, lpf, tmp);
        memcpy(c1[i], tmp, NTMINT*sizeof(double));

        // 下路F2解调
        memcpy(inp, txData[i], NTMINT*sizeof(double));
        filter(inp, bpf2, c2[i]);
        for (int j = 0; j < NTMINT; ++j)
            c2[i][j] *= freq_zero[j];
        memcpy(inp, c2[i], NTMINT*sizeof(double));
        filter(inp, lpf, tmp);
        memcpy(c2[i], tmp, NTMINT*sizeof(double));

//        if (c1[i][500] > c2[i][500])
        if (c1[i][NTMINT/2] > c2[i][NTMINT/2])
            pixel |= 1 << i;
        if (pixel&(1<<i) ^ (txPixel&(1<<i)))
            nerror = nerror + 1;
        npixel = npixel + 1;
    }

#ifdef DEBUG
    ofstream fch_outp1("data/ch_outp1.dat");
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < NTMINT; ++j)
            fch_outp1 << i*TS+j*TS/NTMINT << " " << c1[i][j] << endl;
    }
    fch_outp1.close();
    ofstream fch_outp2("data/ch_outp2.dat");
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < NTMINT; ++j)
            fch_outp2 << i*TS+j*TS/NTMINT << " " << c2[i][j] << endl;
    }
    fch_outp2.close();
#endif

    return pixel;
}

static void filter(const double *in, const double *filt, double *outp)
{
    static fftw_complex ffx[NTMINT], fff[NTMINT], ffm[NTMINT];
    static fftw_plan px, pf1, pf2, pfl;
    fftw_plan ip;

    if (start_up) {
        start_up = false;
        px = fftw_plan_dft_r2c_1d(NTMINT, in, ffx, FFTW_ESTIMATE);
        pf1 = fftw_plan_dft_r2c_1d(NTMINT, bpf1, fff, FFTW_ESTIMATE);
        pf2 = fftw_plan_dft_r2c_1d(NTMINT, bpf2, fff, FFTW_ESTIMATE);
        pfl = fftw_plan_dft_r2c_1d(NTMINT, lpf, fff, FFTW_ESTIMATE);
    }
    fftw_execute(px);
    if (filt == bpf1)
        fftw_execute(pf1);
    else if (filt == bpf2)
        fftw_execute(pf2);
    else
        fftw_execute(pfl);
    for (int j = 0; j < NTMINT; ++j) {
        ffm[j][0] = (ffx[j][0]*fff[j][0] - ffx[j][1]*fff[j][1])/NTMINT;
        ffm[j][1] = (ffx[j][0]*fff[j][1] + ffx[j][1]*fff[j][0])/NTMINT;
    }
    ip = fftw_plan_dft_c2r_1d(NTMINT, ffm, outp, FFTW_ESTIMATE);
    fftw_execute(ip);

    fftw_destroy_plan(ip);
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

void reset_err_pro(int, void *)
{
    npixel = 0;
    nerror = 0;
}
