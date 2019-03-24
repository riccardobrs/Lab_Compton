/*
    ./risol fileset.txt
    
    Il programma apre una TApplication con 2 TGraphErrors che mostrano
    l'andamento della risoluzione (calcolata per i 2 picchi) in dipendenza dalla V_bias.
*/

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

int main (int argc, char ** argv) { //inserire come argv[1] il file dei set per i fit

    string fileSet = argv[1];
    ifstream in (fileSet.c_str());
    if (in.good() == false) {
            cout << "Errore di apertura file" << endl;
            return 1;
    }
    
    TGraphErrors * g1 = new TGraphErrors;
    TGraphErrors * g2 = new TGraphErrors;
    
    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(0.4);
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(0.4);
    
    string line, nometxt;
    char nuovo = 'y';
    double v, mu1, min1, max1, mu2, min2, max2;
    double range [6];
    int i=0;
           
    while(true) {
        
        cout << "Inserire nuovi dati? (y/n)" << endl;
        cin >> nuovo;
        if(nuovo == 'n') break;
        cout << "Inserire file *.txt" << endl;
        cin >> line;
        cout << "Inserire tensione di bias (V)" << endl;
        cin >> v;
        
        while(true) {
            in >> nometxt >> mu1 >> min1 >> max1 >> mu2 >> min2 >> max2;
            if(in.eof() == true) break;
            if (nometxt.compare(fileSet) == 0) break;
        }
        in.close();
        
        range[0] = mu1;
        range[1] = min1;
        range[2] = max1;
        range[3] = mu2;
        range[4] = min2;
        range[5] = max2;
        
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