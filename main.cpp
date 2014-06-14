#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <fftw3.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#define  TS         0.001
// 采样频率
#define  TDELTA     0.000001 
#define  NTMINT     ((int)(TS / TDELTA))
#define  F1         10000
#define  F2         40000
#define  PI         3.1415926535897932

double freq_one[NTMINT], freq_zero[NTMINT];
double deviation = 1;

using namespace std;
using namespace cv;
//using namespace comsim;

static void modulate(const uchar intensity, double *txData[]);
static void signal_channel(double *txData[]);
static uchar demodulate(const double *txData[]);
static double gaussrand(void);

int main(void)
{
//    Mat image = imread("lena.pgm");
//    if (!image.data) {
//        cerr << "Could not open or find the image" << endl;
//        return 1;
//    }
//    namedWindow( "Original image", WINDOW_AUTOSIZE);
//    imshow( "Original image", image);

    Mat src(1, 1, CV_8UC1, Scalar(170));

    srand(time(NULL));

    ofstream ff1("f1.dat");
    ofstream ff2("f2.dat");
    for (int i = 0; i < NTMINT; ++i) {
        freq_one[i] = 5*cos(2*PI*F1*i*TDELTA);
        ff1 << i*TS/NTMINT << " " << freq_one[i] << endl;
        freq_zero[i] = 5*cos(2*PI*F2*i*TDELTA);
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
//
//            uchar rxData;
//            rxData = demodulate(txData);

            for (int i = 0; i < 8; ++i)
                fftw_free(txData[i]);
        }
    }

    return 0;
}

static void modulate(const uchar intensity, double *txData[])
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
