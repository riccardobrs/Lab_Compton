#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TMath.h>

#include "myLib.h"
using namespace std;


double rand_range(double Min, double max) {
  return Min + (max-Min)*rand()/(1. * RAND_MAX);
}

double f_Cauchy(double x) {
  return 1/(M_PI*(1+x*x));
}

//Try&Catch
double rand_TAC(double xMin, double xMax,double yMin,double yMax, double pdf(double(x))) {
  // pdf è una generica funzione densità di probabilità
  double x = 0. , y = 0.;
  do {
  x = rand_range(xMin,xMax);
  y = rand_range(yMin,yMax);
  } while( y > pdf(x) );
  return x;
}

//Distribuzione Esponenziale 
double f_Exp (double x, double parameter) {
  return exp(-x/parameter)/parameter;
}

//Generazione di numeri casuali distribuiti secondo distribuzione Esponenziale
//Funzione Cumulativa Inversa
double rand_FCI_Exp(double parameter) {
  double y;
  do {
    y = rand_range(0,1);
  } while( y == 1.0 ); //devo campionare nell'intervallo [0,1).
  
  return -parameter*log(1-y);
}

double f_gauss(double x, double mean, double sigma) {
  return exp(-pow(x-mean,2)/(2.*pow(sigma,2)))/(sqrt(M_PI*2.) * sigma);
}

//Try&Catch Gaussiana
double rand_TAC_gauss(double xMin, double xMax,double yMin,double yMax, double mean, double sigma ) {
  // la pdf è una la Gaussiana
  double x = 0. , y = 0.;
  do {
  x = rand_range(xMin,xMax);
  y = rand_range(yMin,yMax);
  } while( y > f_gauss(x,mean,sigma) );
  return x;
}


//Hit or Miss Gaussiana
bool HitMiss_gauss(double xMin, double xMax,double yMin,double yMax, double mean, double sigma ) {
  double x = 0. , y = 0.;
  x = rand_range(xMin,xMax);
  y = rand_range(yMin,yMax);
  if ( y < f_gauss(x,mean,sigma) ) return true;
  else return false;
}

//Crude Monte Carlo Gauss
double Crude_MC_Gauss (double min, double max,double mean, double sigma) {
  double random = rand_range(min,max);
  double y = f_gauss(random,mean,sigma);
  return y;
}

//Funzione che determina minimo e massimo di un vector
void Max_Min(std::vector<double> vec,double& min,double&  max) {
  for (int i = 0; i < vec.size() ; i++) {
    if ( i == 0 ) {
      min = vec[i];
      max = vec[i];
    }
    if ( vec[i] < min ) {
      min = vec[i];
    } else if ( vec[i] > max ) {
      max = vec[i];
    }
  }
}

void Media_DevStd(std::vector<double> x,double& media,double& devStd) {
  double mean = 0.;
  double SumSq = 0;
  int N = x.size();
  for ( int i = 0; i < N ; i++ ) {
    mean += x[i];
    SumSq += pow(x[i],2);
  }
  double bias = N / ( N - 1.);
  media = mean / N; 
  devStd = sqrt ( bias * ( SumSq / N - pow(media,2) ) );
}

void Covarianza(std::vector<double> x,std::vector<double> y,double& covarianza) {
  double mean_x = 0.;
  double mean_y = 0.;
  double Sum_xy = 0;
  int N = x.size();
  for ( int i = 0; i < N ; i++ ) {
    mean_x += x[i];
    mean_y += y[i];
    Sum_xy += x[i] * y[i];
  }
  mean_x /= N;
  mean_y /= N;  
  covarianza =  Sum_xy / N - mean_x * mean_y ;
}

double Gauss(double* x, double* par) {
  return par[0] * exp( - pow( (x[0] - par[1] ) ,2) / ( 2 * pow(par[2],2) ) );
}

double Binormale(double* x,double* par) {
  double A = 1. / ( 2 * M_PI * par[2] * par[3] * sqrt( 1 - pow(par[4],2) ) );
  double G = ( 1. / ( 1 - pow(par[4],2) ) * ( pow (( x[0] - par[0] )/par[2],2) - 2 * par[4] * ( (x[0] - par[0] ) / par[2] ) * ( (x[1] - par[1] ) / par[3] ) +  pow (( x[1] - par[1] )/par[3],2) ) );

  return A * exp( - 0.5 * G );
}

void Media_Pesata(std::vector <double> x, std::vector <double> err_x, double& media,double& errore) {
  double result = 0.;
  double weight = 0.;
  double Sum_weight = 0.;
  for ( int i = 0; i < x.size() ; i++) {
    weight = 1. / pow ( err_x[i] ,2);
    result += weight * x[i];
    Sum_weight += weight; 
  }
  media = result / Sum_weight;
  errore = sqrt ( 1. / Sum_weight );
}




