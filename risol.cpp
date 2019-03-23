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
#include "myLib.h"

using namespace std;

int main () {
    
    TGraphErrors * g1 = new TGraphErrors;
    TGraphErrors * g2 = new TGraphErrors;
    
    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(0.4);
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(0.4);
    
    string line;
    char nuovo = 'y';
    double v, min1, max1, min2, max2;
    int i=0;
    while(true) {
        cout << "Inserire nuovi dati? (y/n)" << endl;
        cin >> nuovo;
        if(nuovo == 'n') break;
        cout << "Inserire file *.txt" << endl;
        cin >> line;
        cout << "Inserire tensione di bias (V)" << endl;
        cin >> v;
        cout << "Inserire range fit 1° picco" << endl;
        cin >> min1;
        cin >> max1;
        cout << "Inserire range fit 2° picco" << endl;
        cin >> min2;
        cin >> max2;
        double range [] = {min1, max1, min2, max2};
        g1->SetPoint(i,v,ris(line, range)[0]);
        g1->SetPointError(i,0.2,ris(line, range)[1]);
        g2->SetPoint(i,v,ris(line, range)[2]);
        g2->SetPointError(i,0.2,ris(line, range)[3]);
        i++;
        cout << "Risoluzione 1 = " << ris(line, range)[0] << " +- " << ris(line, range)[1] << endl;
        cout << "Risoluzione 2 = " << ris(line, range)[2] << " +- " << ris(line, range)[3] << endl;
    }

    TApplication* myApp = new TApplication ("myApp", NULL, NULL);
    TCanvas* canva = new TCanvas("c2","canva",100,200,700,500);
    canva->Divide(2,1);
    canva->cd(1);
    g1->Draw("AP");
    canva->cd(2);
    g2->Draw("AP");
    
    myApp->Run();          
    
    return 0;
}