#ifndef MEEG_HH
#define MEEG_HH
#include <TCanvas.h>

#define SAMPLE_INTERVAL 24.0
void doStats(int n, int nmin, int nmax, int *y, int &count, double &center, double &spread);
void doStats(int n, int nmin, int nmax, short int *y, int &count, double &center, double &spread);
void doStats_mean(int n, int nmin, int nmax, int *y, int &count, double &center, double &spread);
void doStats_mean(int n, int nmin, int nmax, short int *y, int &count, double &center, double &spread);
void plotResults(const char *title, const char *name, const char *filename, int n, double *x, double *y, TCanvas *canvas);
void plotResults2(const char *title, const char *name, const char *name2, const char *filename, int n, double *x, double *y, double *y2, TCanvas *canvas);

void myText(Double_t x,Double_t y,const char *text, Double_t tsize,Color_t color);
void myText(Double_t x,Double_t y,TString text, Double_t tsize,Color_t color);


#endif
