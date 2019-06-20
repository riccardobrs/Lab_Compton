#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFitResult.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TMath.h>
#include <TH1D.h>
#include <TFile.h>
#include <TMatrixDSym.h>
#include <TPaveStats.h>
#include "myLib.h"

using namespace std;

TCanvas * create_canva (int num) {
    TCanvas * c = new TCanvas(Form("c%d", num), "", 100, 200, 700, 500);
    return c;
}

double comp (double* x, double* par) {
    //double c = 299792458;
    return 511./( 1+par[0]*(1-cos(x[0]*M_PI/180)));
}


int main (int argc, char ** argv) {
    
    gStyle->SetOptFit(1112);
    gStyle->SetOptStat(11); //print only name of histogram and number of entries
 
    TApplication* myApp = new TApplication ("myApp", NULL, NULL);
    TCanvas* canva = create_canva(1);
    TGraphErrors * g = new TGraphErrors();
    TGraph * g1 = new TGraph();
    g->SetPoint(0,30,465.2090);
    g->SetPoint(1,45,400.8025);
    g->SetPoint(2,60,343.6553);
    g->SetPoint(3,75,314.6184);
    g->SetPoint(4,130,194.9182);
    g->SetPointError(0, 6.,25.5012 );
    g->SetPointError(1, 6.,21.9693);
    g->SetPointError(2, 6.,18.9059);
    g->SetPointError(3, 6.,17.2723);
    g->SetPointError(4, 6.,10.7344);
    g1->SetPoint(0,180,0);
    g1->SetPoint(1,0,0);
    g1->GetXaxis()->SetTitle("#theta");
    g1->GetYaxis()->SetTitle("E' (keV)");
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.6);
    g->SetMarkerColor(kRed);

    TF1 * compton = new TF1("compton", comp, 0, 180, 1);
    compton->SetLineColor(kBlue);
    canva->cd();
    g1->Draw("AP");
    g->Draw("P SAME");
    TFitResultPtr r1 = g->Fit("compton", "R S", "sames");
    TMatrixDSym covariance_matrix_1 = r1 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_1 = r1 -> GetCorrelationMatrix();
    g1->GetYaxis()->SetRangeUser(0,180);
    g1->GetYaxis()->SetRangeUser(0,550);
    
    myApp -> Run();
    
    return 0;

}