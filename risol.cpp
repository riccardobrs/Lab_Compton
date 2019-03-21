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
    
    TApplication* myApp = new TApplication ("myApp", NULL, NULL);
    TCanvas* canva = new TCanvas("c2","canva",100,200,700,500);  
    TGraph * g1 = new TGraph;
    TGraph * g2 = new TGraph;
    
    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(0.4);
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(0.4);
    
    string line;
    vector <double> x,y;
    char nuovo = 'y';
    int i=0;
    while(true) {
        cout << "Inserire nuovi dati? (y/n)" << endl;
        cin >> nuovo;
        if(nuovo == 'n') break;
        cin >> line;
        x.push_back(ris(line)[0]);
        y.push_back(ris(line)[1]);
        g1->SetPoint(i, 1000, x[i]);
        g2->SetPoint(i, 1000, y[i]);
        i++;
        cout << "Risoluzione 1 = " << x[0] << "\t\t" << "Risoluzione 2 = " << y[0] << endl;
    }
    
    canva->Divide(2,1);
    canva->cd(1);
    g1->Draw("AP");
    canva->cd(2);
    g2->Draw("AP");
    
    myApp->Run();          
    
    return 0;
}