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
    TCanvas* c2 = new TCanvas("c2","2nd column",100,200,700,500);
    TCanvas* c3 = new TCanvas("c3","3rd column",100,200,700,500);
    TCanvas* c4 = new TCanvas("c4","4th column",100,200,700,500);
    TCanvas* c5 = new TCanvas("c5","5th column",100,200,700,500);
    TCanvas* c6 = new TCanvas("c6","6th column",100,200,700,500);
    
    TMultiGraph * MG2 = new TMultiGraph();
    TMultiGraph * MG3 = new TMultiGraph();
    TMultiGraph * MG4 = new TMultiGraph();
    TMultiGraph * MG5 = new TMultiGraph();
    TMultiGraph * MG6 = new TMultiGraph();
    
    string fileInput = argv[1];
    ifstream in (fileInput.c_str());
    if (in.good() == false) {
            cout << "Errore di apertura file" << endl;
            return 1;
    }
    ofstream outfile ("fitResult_"+fileInput);
    
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
    
    TGraph * gr2x = new TGraph;
    TGraph * gr2y = new TGraph;
    TGraph * gr2z = new TGraph;
    TGraph * gr3x = new TGraph;
    TGraph * gr3y = new TGraph;
    TGraph * gr3z = new TGraph;
    TGraph * gr4x = new TGraph;
    TGraph * gr4y = new TGraph;
    TGraph * gr4z = new TGraph;
    TGraph * gr5x = new TGraph;
    TGraph * gr5y = new TGraph;
    TGraph * gr5z = new TGraph;
    TGraph * gr6x = new TGraph;
    TGraph * gr6y = new TGraph;
    TGraph * gr6z = new TGraph;
    
    for(int i=0; i<(5500/5); i++) {
        gr2x->SetPoint(i, vx[i], va[i]);
        gr3x->SetPoint(i, vx[i], vb[i]);
        gr4x->SetPoint(i, vx[i], vc[i]);
        gr5x->SetPoint(i, vx[i], vd[i]);
        gr6x->SetPoint(i, vx[i], ve[i]);
    }
    for(int i=(5500/5); i<(6500/5); i++) {
        gr2y->SetPoint(i-(5500/5), vx[i], va[i]);
        gr3y->SetPoint(i, vx[i], vb[i]);
        gr4y->SetPoint(i, vx[i], vc[i]);
        gr5y->SetPoint(i, vx[i], vd[i]);
        gr6y->SetPoint(i, vx[i], ve[i]);
    }
    for(int i=(6500/5); i<vx.size(); i++) {
        gr2z->SetPoint(i-(6500/5), vx[i], va[i]);
        gr3z->SetPoint(i, vx[i], vb[i]);
        gr4z->SetPoint(i, vx[i], vc[i]);
        gr5z->SetPoint(i, vx[i], vd[i]);
        gr6z->SetPoint(i, vx[i], ve[i]);
    }
    
    gr2x->SetMarkerStyle(20);
    gr2x->SetMarkerSize(0.34);
    gr2x->SetMarkerColor(4);
    gr2y->SetMarkerStyle(20);
    gr2y->SetMarkerSize(0.34);
    gr2y->SetMarkerColor(4);
    gr2z->SetMarkerStyle(20);
    gr2z->SetMarkerSize(0.34);
    gr2z->SetMarkerColor(4);
    
    gr3x->SetMarkerStyle(20);
    gr3x->SetMarkerSize(0.34);
    gr3x->SetMarkerColor(4);
    gr3y->SetMarkerStyle(20);
    gr3y->SetMarkerSize(0.34);
    gr3y->SetMarkerColor(4);
    gr3z->SetMarkerStyle(20);
    gr3z->SetMarkerSize(0.34);
    gr3z->SetMarkerColor(4);
    
    gr4x->SetMarkerStyle(20);
    gr4x->SetMarkerSize(0.34);
    gr4x->SetMarkerColor(4);
    gr4y->SetMarkerStyle(20);
    gr4y->SetMarkerSize(0.34);
    gr4y->SetMarkerColor(4);
    gr4z->SetMarkerStyle(20);
    gr4z->SetMarkerSize(0.34);
    gr4z->SetMarkerColor(4);
    
    gr5x->SetMarkerStyle(20);
    gr5x->SetMarkerSize(0.34);
    gr5x->SetMarkerColor(4);
    gr5y->SetMarkerStyle(20);
    gr5y->SetMarkerSize(0.34);
    gr5y->SetMarkerColor(4);
    gr5z->SetMarkerStyle(20);
    gr5z->SetMarkerSize(0.34);
    gr5z->SetMarkerColor(4);
    
    gr6x->SetMarkerStyle(20);
    gr6x->SetMarkerSize(0.34);
    gr6x->SetMarkerColor(4);
    gr6y->SetMarkerStyle(20);
    gr6y->SetMarkerSize(0.34);
    gr6y->SetMarkerColor(4);
    gr6z->SetMarkerStyle(20);
    gr6z->SetMarkerSize(0.34);
    gr6z->SetMarkerColor(4);
    
    TF1 * f1 = new TF1 ("Gauss1", Gauss, 2240, 2700, 3);
    f1 -> SetParameter(1,2500);
    f1 -> SetParameter(2, 100);
    f1 -> SetParName (0, "Amp_{1}" );
    f1 -> SetParName (1, "#mu_{1}" );
    f1 -> SetParName (2, "#sigma_{1}" );
    f1 -> SetLineColor(kRed);
    
    TF1 * f2 = new TF1 ("Gauss2", Gauss, 5500, 6500, 3);
    f2 -> SetParameter(0,278);
    f2 -> SetParameter(1,5943);
    f2 -> SetParameter(2, 268);
    f2 -> SetParName (0, "Amp_{2}" );
    f2 -> SetParName (1, "#mu_{2}" );
    f2 -> SetParName (2, "#sigma_{2}" );
    f2 -> SetLineColor(6);
    
    MG2->Add(gr2x);
    MG2->Add(gr2y);
    MG2->Add(gr2z);
    MG3->Add(gr3x);
    MG3->Add(gr3y);
    MG3->Add(gr3z);
    MG4->Add(gr4x);
    MG4->Add(gr4y);
    MG4->Add(gr4z);
    MG5->Add(gr5x);
    MG5->Add(gr5y);
    MG5->Add(gr5z);
    MG6->Add(gr6x);
    MG6->Add(gr6y);
    MG6->Add(gr6z);
    
    c2->cd();
    MG2->Draw("AP");
    TFitResultPtr r2_1 = gr2x->Fit("Gauss1", "S");
    TMatrixDSym covariance_matrix_2_1 = r2_1 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_2_1 = r2_1 -> GetCorrelationMatrix();
    TFitResultPtr r2_2 = gr2y->Fit("Gauss2", "S");
    TMatrixDSym covariance_matrix_2_2 = r2_2 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_2_2 = r2_2 -> GetCorrelationMatrix();

    c3->cd();
    MG3->Draw("AP");
    TFitResultPtr r3_1 = gr3x->Fit("Gauss1", "S");
    TMatrixDSym covariance_matrix_3_1 = r3_1 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_3_1 = r3_1 -> GetCorrelationMatrix();
    TFitResultPtr r3_2 = gr3y->Fit("Gauss2", "S");
    TMatrixDSym covariance_matrix_3_2 = r3_2 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_3_2 = r3_2 -> GetCorrelationMatrix();

    c4->cd();
    MG4->Draw("AP");
    TFitResultPtr r4_1 = gr4x->Fit("Gauss1", "S");
    TMatrixDSym covariance_matrix_4_1 = r4_1 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_4_1 = r4_1 -> GetCorrelationMatrix();
    TFitResultPtr r4_2 = gr4y->Fit("Gauss2", "S");
    TMatrixDSym covariance_matrix_4_2 = r4_2 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_4_2 = r4_2 -> GetCorrelationMatrix();
 
    c5->cd();
    MG5->Draw("AP");
    TFitResultPtr r5_1 = gr5x->Fit("Gauss1", "S");
    TMatrixDSym covariance_matrix_5_1 = r5_1 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_5_1 = r5_1 -> GetCorrelationMatrix();
    TFitResultPtr r5_2 = gr5y->Fit("Gauss2", "S");
    TMatrixDSym covariance_matrix_5_2 = r5_2 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_5_2 = r5_2 -> GetCorrelationMatrix();
 
    c6->cd();
    MG6->Draw("AP");
    TFitResultPtr r6_1 = gr6x->Fit("Gauss1", "S");
    TMatrixDSym covariance_matrix_6_1 = r6_1 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_6_1 = r6_1 -> GetCorrelationMatrix();
    TFitResultPtr r6_2 = gr6y->Fit("Gauss2", "S");
    TMatrixDSym covariance_matrix_6_2 = r6_2 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_6_2 = r6_2 -> GetCorrelationMatrix();
    /*
    outfile << "\nMatrice di covarianza " << endl;
    for (int i = 0; i < 3; i++) {
        outfile << setw(15) << f1->GetParName(i);
    }
    outfile << endl;
    
    for (int i=0; i<3; i++) {
        outfile << f1->GetParName(i);
        for (int j=0; j<3; j++) {
        double sigma_ij = covariance_matrix_2_1(i,j);
        outfile << setw(15)<< sigma_ij;
        }
    outfile << endl;
    }
	
    outfile << "\nMatrice di correlazione " << endl;
    for (int i = 0; i < 3; i++) {
        outfile << setw(15) << f1->GetParName(i);
    }
    outfile << endl;
       
    for (int i=0; i<3; i++) {
        outfile << f1->GetParName(i);
        for (int j=0; j<3; j++) {
	    double ro_ij = correlation_matrix_2_1(i,j);
	    outfile << setw(15) << ro_ij;
        }
        outfile << endl;
    }

    outfile << endl << "VALORI RESTITUITI DA ROOT" << endl;
    outfile << "Gradi di libertà = " << f1->GetNDF() << endl;
    outfile << "Chi2 = " << f1->GetChisquare() << endl;
    outfile << "Reduced Chi2 = " << f1->GetChisquare() / f1->GetNDF() << endl;
    outfile << "p-value = " << f1->GetProb() << endl;	

    outfile << "** STIMA PARAMETRI**" << endl;
    outfile << "Amp_{1} = " << f1 -> GetParameter(0) << " +- " << f1 -> GetParError(0) << endl;
    outfile << "#mu_{1} = " << f1 -> GetParameter(1) << " +- " << f1 -> GetParError(1) << endl;
    outfile << "#sigma_{1} = " << f1 -> GetParameter(2) << " +- " << f1 -> GetParError(2) << endl;
    */
    myApp -> Run();  
    
    return 0;

}