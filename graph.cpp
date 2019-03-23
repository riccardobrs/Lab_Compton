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

int main (int argc, char ** argv) {
    
    gStyle->SetOptFit(1112);
 
    TApplication* myApp = new TApplication ("myApp", NULL, NULL);
    TCanvas* canva = new TCanvas("c2","2nd column",100,200,700,500);    
    
    string fileInput = argv[1];
    ifstream in (fileInput.c_str());
    if (in.good() == false) {
            cout << "Errore di apertura file" << endl;
            return 1;
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
            NBin+=5;
        }
        N++;
    }
    in.close();
    Max_Min(vx, min, max);
    cout << "NBIN = " << NBin << endl;
    
    TH1D * histo = new TH1D("Histogram", "", NBin-7, min, max);
    
    for(int i=0; i<vx.size(); i++) {
        for(int j =0; j<va[i];j++) histo->Fill(vx[i]);
        for(int j =0; j<vb[i];j++) histo->Fill(vx[i]+1);
        for(int j =0; j<vc[i];j++) histo->Fill(vx[i]+2);
        for(int j =0; j<vd[i];j++) histo->Fill(vx[i]+3);
        for(int j =0; j<ve[i];j++) histo->Fill(vx[i]+4);
    }

    canva->cd();
    histo->Draw();

    myApp -> Run();  
    
    return 0;

}