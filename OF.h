#include "TRandom.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TDecompChol.h"
#include "TMatrixDSymEigen.h"
#include <sys/time.h>
#include <time.h>
#include <vector>
#include <cmath>
#include "TComplex.h"
#include "TVirtualFFT.h"
#include <fstream>
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TMatrixD.h"
#include "TMinuit.h"
#include "TLegend.h"

using namespace std;

void ComplexToRealIFFT(const vector<TComplex>& inComp, vector<double>& outRe);

void RealToComplexFFT(const vector<double>& pulsevector, vector<TComplex>& outComp);

void RealToComplexFFT(const vector<double>& pulsevector, vector<TComplex>& outComp);

void generate_template(vector<double> &template_vec_9152);

void generate_whitenoise_time(double *noise9152, double *noise131072);

void getNoise_time(double *noise9152, double *noise131072, double *noisePSD);

void PSD2VarianceVec(double *PSD, double *var);

void getHarmonic(double *noise_9152, double *noise_131072, int mode);

void makeGMatrix();

void makeCovMatrix();

void compressTime(double *time131072, double *time9152);

void findAmplitude_9152_full_covariance(double* template_time, double *signal_time, TMatrixD* covInv, TMatrixD* cov, double *a, double *chi2);

void findAmplitude_2template_9152_full_covariance(double* template_time_slow, double* template_time_fast, double *signal_time, TMatrixD* covInv, double *a, double *b, double *chi2);


