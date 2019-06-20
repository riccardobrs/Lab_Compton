#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <TMultiGraph.h>
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

double gaus_pol2(double* x, double* par) {
  double g = par[0] * exp( - pow( (x[0] - par[1] ) ,2) / ( 2 * pow(par[2],2) ) );
  double p = par[3]+x[0]*par[4]+x[0]*x[0]*par[5];
  return g + p;
}

double gaus_pol1(double* x, double* par) {
  double g = par[0] * exp( - pow( (x[0] - par[1] ) ,2) / ( 2 * pow(par[2],2) ) );
  double r = par[3]+x[0]*par[4];
  return g + r;
}

double gaussian (double* x, double* par) {
  double g = par[0] * exp( - pow( (x[0] - par[1] ) ,2) / ( 2 * pow(par[2],2) ) );
  return g;
}

double gaus2 (double* x, double* par) {
  double g1 = par[0] * exp( - pow( (x[0] - par[1] ) ,2) / ( 2 * pow(par[2],2) ) );
  double g2 = par[3] * exp( - pow( (x[0] - par[4] ) ,2) / ( 2 * pow(par[5],2) ) );
  return g1+g2;
}

double crystal(double* x, double* par) {
    double g = par[0] * exp( - pow( (x[0] - par[1] ) ,2) / ( 2 * pow(par[2],2) ) );
    double e = par[3] * exp(par[4]*x[0]) + par[5];
    return g+e;
}

vector <double> ris (string filename, double fileset []) { //la funzione restituisce |Ris Picco1|Err Ris1|Ris Picco2|Err Ris2|
     
    ifstream in (filename.c_str());
    if (in.good() == false) {
            cout << "Errore di apertura file" << endl;
            exit(EXIT_FAILURE);
    }

    double x, a, b, c, d, e, min, max;
    string line;
    int N = 0, NBin = 1;
    vector <double> vx, va, vb, vc, vd, ve;
    
    while (getline(in,line)) {
        if ( N > 24 ) {
            in >> x >> a >> b >> c >> d >> e;
            if (in.eof() == true) break;
            vx.push_back(x);
            va.push_back(a);
            vb.push_back(b);
            vc.push_back(c);
            vd.push_back(d);
            ve.push_back(e);
            NBin++;
        }
        N++;
    }
    in.close();
    Max_Min(vx, min, max);
    
    TH1D * histo = new TH1D("nome", "nome", NBin-7, min, max);
    
    for(int i=0; i<vx.size(); i++) {
        for(int j =0; j<va[i];j++) histo->Fill(vx[i]);
        for(int j =0; j<vb[i];j++) histo->Fill(vx[i]+1);
        for(int j =0; j<vc[i];j++) histo->Fill(vx[i]+2);
        for(int j =0; j<vd[i];j++) histo->Fill(vx[i]+3);
        for(int j =0; j<ve[i];j++) histo->Fill(vx[i]+4);
    }

    TF1 * f1 = new TF1 ("511_Gaus+pol2", gaus_pol2, fileset[1], fileset[2], 6);
    f1 -> SetParameter(1,fileset[0]);
    f1 -> SetParameter(2, 100);
    f1 -> SetParameter(3, 1620);
    f1 -> SetParameter(4, -0.95);
    f1 -> SetParameter(5, 1.59e-04);
    f1 -> SetParName (0, "Amp_{1}" );
    f1 -> SetParName (1, "#mu_{1}" );
    f1 -> SetParName (2, "#sigma_{1}" );
    f1 -> SetParName (3, "a_{0}" );
    f1 -> SetParName (4, "a_{1}" );
    f1 -> SetParName (5, "a_{2}" );
    f1 -> SetLineColor(kRed);
    
    TF1 * f2 = new TF1 ("1274_Gaus+pol2", gaus_pol2, fileset[4], fileset[5], 6);
    f2 -> SetParameter(0,278);
    f2 -> SetParameter(1,fileset[3]);
    f2 -> SetParameter(2, 268);
    f2 -> SetParName (0, "Amp_{2}" );
    f2 -> SetParName (1, "#mu_{2}" );
    f2 -> SetParName (2, "#sigma_{2}" );
    f2 -> SetParName (3, "b_{0}" );
    f2 -> SetParName (4, "b_{1}" );
    f2 -> SetParName (5, "b_{2}" );
    f2 -> SetLineColor(6);
    
    histo->Fit("511_Gaus+pol2", "R");
    histo->Fit("1274_Gaus+pol2", "R");
    /*TFitResultPtr r1 = histo->Fit("511_Gaus+pol2", "R S");
    TMatrixDSym covariance_matrix_2_1 = r1 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_2_1 = r1 -> GetCorrelationMatrix();
    TFitResultPtr r2 = histo->Fit("1274_Gaus+pol2", "R S +");
    TMatrixDSym covariance_matrix_2_2 = r2 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_2_2 = r2 -> GetCorrelationMatrix();*/
    
    vector <double> resolution;
    double A = f1->GetParameter(1);
    double B = f1->GetParameter(2);
    double C = f2->GetParameter(1);
    double D = f2->GetParameter(2);
    double err_A = f1->GetParError(1);
    double err_B = f1->GetParError(2);
    double err_C = f2->GetParError(1);
    double err_D = f2->GetParError(2);
    resolution.push_back(2.35*B/A);
    resolution.push_back(2.35*sqrt(pow(err_B/A,2)+pow((B*err_A)/(A*A),2)));
    resolution.push_back(2.35*D/C);
    resolution.push_back(2.35*sqrt(pow(err_D/C,2)+pow((D*err_C)/(C*C),2)));
    
    return resolution;
}




