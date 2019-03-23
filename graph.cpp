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
//    ofstream outfile ("fitResult_"+fileInput);
    
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
    
    TH1D * histo = new TH1D("nome", "nome", NBin-7, min, max);
    
    for(int i=0; i<vx.size(); i++) {
        for(int j =0; j<va[i];j++) histo->Fill(vx[i]);
        for(int j =0; j<vb[i];j++) histo->Fill(vx[i]+1);
        for(int j =0; j<vc[i];j++) histo->Fill(vx[i]+2);
        for(int j =0; j<vd[i];j++) histo->Fill(vx[i]+3);
        for(int j =0; j<ve[i];j++) histo->Fill(vx[i]+4);
    }
/*    
    TF1 * f1 = new TF1 ("511_Gaus+pol2", gaus_pol2, 1800, 3128, 6);
    f1 -> SetParameter(1,2500);
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
    
    TF1 * f2 = new TF1 ("1274_Gaus+pol2", gaus_pol2, 5200, 6700, 6);
    f2 -> SetParameter(0,278);
    f2 -> SetParameter(1,5943);
    f2 -> SetParameter(2, 268);
    f2 -> SetParName (0, "Amp_{2}" );
    f2 -> SetParName (1, "#mu_{2}" );
    f2 -> SetParName (2, "#sigma_{2}" );
    f2 -> SetParName (3, "b_{0}" );
    f2 -> SetParName (4, "b_{1}" );
    f2 -> SetParName (5, "b_{2}" );
    f2 -> SetLineColor(6);
*/      
    canva->cd();
    histo->Draw();
/*
    TFitResultPtr r1 = histo->Fit("511_Gaus+pol2", "R S");
    TMatrixDSym covariance_matrix_1 = r1 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_1 = r1 -> GetCorrelationMatrix();
    TFitResultPtr r2 = histo->Fit("1274_Gaus+pol2", "R S +");
    TMatrixDSym covariance_matrix_2 = r2 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_2 = r2 -> GetCorrelationMatrix();

    outfile << "\nMatrice di covarianza 1° PICCO" << endl;
    for (int i = 0; i < 6; i++) {
        outfile << setw(15) << f1->GetParName(i);
    }
    outfile << endl;
    
    for (int i=0; i<6; i++) {
        outfile << f1->GetParName(i);
        for (int j=0; j<6; j++) {
        double sigma_ij = covariance_matrix_1(i,j);
        outfile << setw(15)<< sigma_ij;
        }
    outfile << endl;
    }
	
    outfile << "\nMatrice di correlazione 1° PICCO" << endl;
    for (int i = 0; i < 6; i++) {
        outfile << setw(15) << f1->GetParName(i);
    }
    outfile << endl;
       
    for (int i=0; i<6; i++) {
        outfile << f1->GetParName(i);
        for (int j=0; j<6; j++) {
	    double ro_ij = correlation_matrix_1(i,j);
	    outfile << setw(15) << ro_ij;
        }
        outfile << endl;
    }
    
    outfile << "\nMatrice di covarianza 2° PICCO" << endl;
    for (int i = 0; i < 6; i++) {
        outfile << setw(15) << f2->GetParName(i);
    }
    outfile << endl;
    
    for (int i=0; i<6; i++) {
        outfile << f2->GetParName(i);
        for (int j=0; j<6; j++) {
        double sigma2_ij = covariance_matrix_2(i,j);
        outfile << setw(15)<< sigma2_ij;
        }
    outfile << endl;
    }
	
    outfile << "\nMatrice di correlazione 2° PICCO" << endl;
    for (int i = 0; i < 6; i++) {
        outfile << setw(15) << f2->GetParName(i);
    }
    outfile << endl;
       
    for (int i=0; i<6; i++) {
        outfile << f2->GetParName(i);
        for (int j=0; j<6; j++) {
	    double ro2_ij = correlation_matrix_2(i,j);
	    outfile << setw(15) << ro2_ij;
        }
        outfile << endl;
    }

    outfile << endl << "*********** PRIMO PICCO **************" << endl;
    outfile << endl << "VALORI RESTITUITI DA ROOT" << endl;
    outfile << "Gradi di libertà = " << f1->GetNDF() << endl;
    outfile << "Chi2 = " << f1->GetChisquare() << endl;
    outfile << "Reduced Chi2 = " << f1->GetChisquare() / f1->GetNDF() << endl;
    outfile << "p-value = " << f1->GetProb() << endl;	

    outfile << "** STIMA PARAMETRI**" << endl;
    outfile << "Amp    //ofstream outfile ("fitResult_"+fileInput);_{1} = " << f1 -> GetParameter(0) << " +- " << f1 -> GetParError(0) << endl;
    outfile << "#mu_{1} = " << f1 -> GetParameter(1) << " +- " << f1 -> GetParError(1) << endl;
    outfile << "#sigma_{1} = " << f1 -> GetParameter(2) << " +- " << f1 -> GetParError(2) << endl;
    outfile << "a_{0} = " << f1 -> GetParameter(3) << " +- " << f1 -> GetParError(3) << endl;
    outfile << "a_{1} = " << f1 -> GetParameter(4) << " +- " << f1 -> GetParError(4) << endl;
    outfile << "a_{2} = " << f1 -> GetParameter(5) << " +- " << f1 -> GetParError(5) << endl;
    
    outfile << endl << "*********** SECONDO PICCO **************" << endl;
    outfile << endl << "VALORI RESTITUITI DA ROOT" << endl;
    outfile << "Gradi di libertà = " << f2->GetNDF() << endl;
    outfile << "Chi2 = " << f2->GetChisquare() << endl;
    outfile << "Reduced Chi2 = " << f2->GetChisquare() / f2->GetNDF() << endl;
    outfile << "p-value = " << f2->GetProb() << endl;	

    outfile << "** STIMA PARAMETRI**" << endl;
    outfile << "Amp_{2} = " << f2 -> GetParameter(0) << " +- " << f2 -> GetParError(0) << endl;
    outfile << "#mu_{2} = " << f2 -> GetParameter(1) << " +- " << f2 -> GetParError(1) << endl;
    outfile << "#sigma_{2} = " << f2 -> GetParameter(2) << " +- " << f2 -> GetParError(2) << endl;
    outfile << "b_{0} = " << f2 -> GetParameter(3) << " +- " << f2 -> GetParError(3) << endl;
    outfile << "b_{1} = " << f2 -> GetParameter(4) << " +- " << f2 -> GetParError(4) << endl;
    outfile << "b_{2} = " << f2 -> GetParameter(5) << " +- " << f2 -> GetParError(5) << endl;
*/
    myApp -> Run();  
    
    return 0;

}