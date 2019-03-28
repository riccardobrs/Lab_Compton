/*
    ./fit <nome_del_file_dati>.Txt
    
    Il programma necessita inoltre la presenza di un file "fitset2.txt" che contiene i set per fit
    
    Organizzazione di fitset2.txt
    
    | <nome_del_file_dati>.Txt | mu1 | min1 | max 1 | mu2 | min2 | max2 |
    
    Spiegazione variabili
    
    mu1 = media primo picco
    min1 = min del range di fit primo picco
    max2 = max del range di fit primo picco
    mu2, min2, max2 analoghi per il secondo picco
    
    Il programma salva di default gli istogrammi fittati in formato png. Tuttavia viene aperta
    anche una TApplication per eventuali modifiche ad hoc.
    
    Inoltre aggiorna il contenuto di un file "resolutions.txt" da passare come argv[1] di "risol.cpp"
    NB: il file "resolutions.txt" DEVE essere già esistente nella cartella, perché oltre ad essere aperto
    in modalità scrittura, viene aperto anche in modalità lettura. Alla prima esecuzione di fit.cpp dovrà
    quindi essere creato il file di testo vuoto.
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

int main (int argc, char ** argv) {
    
    //gStyle->SetOptFit(1112);
 
    TApplication* myApp = new TApplication ("myApp", NULL, NULL);
    TCanvas* canva = new TCanvas("canva","canva",100,200,700,500);    
    
    string fileInput = argv[1];
    ifstream in (fileInput.c_str());
    if (in.good() == false) {
            cout << "Errore di apertura file" << endl;
            return 1;
    }
//    ofstream outfile ("fitResult_"+fileInput);
    
    double x, a, b, c, d, e, min, max;
    string line, txtdata;
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
    
    double mu1, mu2, min1, min2, max1, max2;
    
    string fileSet = "fitset2.txt";
    ifstream fs(fileSet.c_str());
    if (fs.good() == false) {
            cout << "Errore di apertura file" << endl;
            return 1;
    }
    while(true) {
        fs >> txtdata >> mu1 >> min1 >> max1 >> mu2 >> min2 >> max2;
        if(fs.eof() == true) break;
        if(txtdata.compare(fileInput) == 0) break;
    }
    fs.close();
    
    TF1 * f1 = new TF1 ("511_Gaus+pol2", gaus_pol2, min1, max1, 6);
    f1 -> SetParameter(1,mu1);
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
    
    TF1 * f2 = new TF1 ("1274_Gaus+pol2", gaus_pol2, min2, max2, 6);
    //f2 -> SetParameter(0,278);
    f2 -> SetParameter(1,mu2);
    f2 -> SetParameter(2, 268);
    f2 -> SetParName (0, "Amp_{2}" );
    f2 -> SetParName (1, "#mu_{2}" );
    f2 -> SetParName (2, "#sigma_{2}" );
    f2 -> SetParName (3, "b_{0}" );
    f2 -> SetParName (4, "b_{1}" );
    f2 -> SetParName (5, "b_{2}" );
    f2 -> SetLineColor(6);
      
    canva->cd();
    histo->Draw();
    string file_in;
    const char * argv1_name;
    double v_bias;
    double v_bias_err = 0.2; // errore fissato "a mano"
    file_in = fileInput.replace(16, 4, ".png");
    if(fileInput[0]==0) {
        string fileInput2 = fileInput;
        fileInput.clear();
        for(int j=0; j<3; j++)
            fileInput[j]=fileInput2[j+1];
        v_bias = stoi(fileInput);
    }
    else
        v_bias = stoi(fileInput.replace(4, 15, ""));
    argv1_name = file_in.c_str();
    
    TFitResultPtr r1 = histo->Fit("511_Gaus+pol2", "R S");
    TMatrixDSym covariance_matrix_1 = r1 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_1 = r1 -> GetCorrelationMatrix();
    TFitResultPtr r2 = histo->Fit("1274_Gaus+pol2", "R S +");
    TMatrixDSym covariance_matrix_2 = r2 -> GetCovarianceMatrix();
    TMatrixDSym correlation_matrix_2 = r2 -> GetCorrelationMatrix();
    
    double res1 = (2.35*f1->GetParameter(2)) / (f1->GetParameter(1));
    double res2 = (2.35*f2->GetParameter(2)) / (f2->GetParameter(1));
    double res1_err = 2.35 * sqrt(pow(f1->GetParError(2)/f1->GetParameter(1),2)+pow(f1->GetParError(1)*f1->GetParameter(2)/pow(f1->GetParameter(1),2),2));
    double res2_err = 2.35 * sqrt(pow(f2->GetParError(2)/f2->GetParameter(1),2)+pow(f2->GetParError(1)*f2->GetParameter(2)/pow(f2->GetParameter(1),2),2));
    
/*
    Aggiungo ad un file "resolutions.txt" delle righe così costruite:
    
    |V_bias|V_bias_err|Risoluz_picco1|Risoluz_picco1_err|Risoluz_picco2|Risoluz_picco2_err|
*/

    string resolutions = "resolutions.txt";
    bool notwrite = false;
    double o1, o2, o3, o4, o5, o6;
    
    fstream out_in (resolutions.c_str(), ios::in); //apro il file in modalità input (lettura)
    if (out_in.good() == false) {
            cout << "File 'resolutions.txt' inesistente o corrotto" << endl;
            return 1;
    }
    while(true) {
        out_in >> o1 >> o2 >> o3 >> o4 >> o5 >> o6;
        if(out_in.eof() == true) break;
        if(o1==v_bias) notwrite = true; //se è già presente non va scritto nuovamente
    }
    out_in.close();
    fstream out_app (resolutions.c_str(), ios::app); //apro il file in modalità append (scrittura: aggiungo righe alle preesistenti)
    if(notwrite == false)
        out_app << v_bias << "\t" << v_bias_err << "\t" << res1 << "\t" << res1_err << "\t" << res2 << "\t" << res2_err << endl;
    out_app.close();
/*
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
    canva->Print(argv1_name, "png");
    myApp -> Run();
    
    return 0;

}