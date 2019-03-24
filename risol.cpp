/*
    ./risol resolutions.txt
    
    Il programma apre una TApplication con 2 TGraphErrors che mostrano
    l'andamento della risoluzione (calcolata per i 2 picchi) in dipendenza dalla V_bias.
    
    NB è necessario come argv[1] "resolutions.txt", che aggiornato con i nuovi parametri stimati
    ad ogni esecuzione di "fit.cpp"
    
    "resolutions.txt" è così strutturato:
    
    |V_bias|V_bias_err|Risoluz_picco1|Risoluz_picco1_err|Risoluz_picco2|Risoluz_picco2_err|
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

    TGraphErrors * g1 = new TGraphErrors;
    TGraphErrors * g2 = new TGraphErrors;
    
    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(0.7);
    g1->SetMarkerColor(kBlue);
    g1->GetXaxis()->SetTitle("V_{bias} (V)");
    g1->GetYaxis()->SetTitle("Resolution");
    g1->GetYaxis()->SetTitleOffset(1.57);
    
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(0.7);
    g2->SetMarkerColor(kRed);
    g2->GetXaxis()->SetTitle("V_{bias} (V)");
    g2->GetYaxis()->SetTitle("Resolution");
    g2->GetYaxis()->SetTitleOffset(1.57);
    
    double v_bias, v_bias_err, res1, res1_err, res2, res2_err;
    int i = 0;
    
    string inputFile = argv[1];
    ifstream in (inputFile.c_str());
    if (in.good() == false) {
            cout << "Errore di apertura file" << endl;
            return 1;
    }
                  
    while(true) {
        in >> v_bias >> v_bias_err >> res1 >> res1_err >> res2 >> res2_err;
        if(in.eof() == true) break;
        g1->SetPoint(i, v_bias, res1);
        g1->SetPointError(i, v_bias_err, res1_err);
        g2->SetPoint(i, v_bias, res2);
        g2->SetPointError(i, v_bias_err, res2_err);
        i++;
    }
    in.close();

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