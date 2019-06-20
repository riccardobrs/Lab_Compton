/*
    ./fit fitset2.Txt
    
    Organizzazione di fitset2.txt
    
    | <nome_del_file_dati>.Txt | mu1 | min1 | max 1 | reb | min2 | max2 |
    
    Spiegazione variabili
    
    mu1 = media primo picco
    min1 = min del range di fit primo picco
    max2 = max del range di fit primo picco
    reb, min2, max2 analoghi per il secondo picco
    
    Il programma salva di default gli istogrammi fittati in formato png. Tuttavia viene aperta
    anche una TApplication per eventuali modifiche ad hoc.
    
    Inoltre aggiorna il contenuto di un file "resolutions2.txt" da passare come argv[1] di "risol.cpp"
    NB: il file "resolutions2.txt" DEVE essere già esistente nella cartella, perché oltre ad essere aperto
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
#include <TFile.h>
#include <TMatrixDSym.h>
#include <TPaveStats.h>
#include "myLib.h"

using namespace std;

TCanvas * create_canva (int num) {
    TCanvas * c = new TCanvas(Form("c%d", num), "", 100, 200, 700, 500);
    return c;
}

TH1D * create_histo (int num, int nbin, double min, double max) {
    TH1D * h = new TH1D ("Picco ambientale", "", nbin, min, max);
    return h;
}

vector <TPaveStats *> create_pave (TH1D * h1, TH1D * h2, int num) { //Creo un box statistico per ogni picco
    TPaveStats * p1 = (TPaveStats *)h1 -> FindObject("stats");
    p1->SetName(Form("pave1_%d",num));
    p1->SetY1NDC(0.5); //per impostare la posizione del lato inferiore del box statistico p1
    p1->SetY2NDC(0.9); //per impostare la posizione del lato superiore del box statistico p1
    p1->SetTextColor(kRed);
    TPaveStats * p2 = (TPaveStats *)h2 -> FindObject("stats");
    p2->SetName(Form("pave2_%d",num));
    p2->SetY1NDC(0.1);
    p2->SetY2NDC(0.5);
    p2->SetTextColor(6);
    vector <TPaveStats *> p;
    p.push_back(p1);
    p.push_back(p2);
    return p;
}

int main (int argc, char ** argv) {
    
    gStyle->SetOptFit(1112);
    gStyle->SetOptStat(11); //print only name of histogram and number of entries
 
    TApplication* myApp = new TApplication ("myApp", NULL, NULL);
    TCanvas* canva;
    TH1D * histo;
    TH1D * histo2; //utilizzato successivamente per creare doppie box statistiche
    
    string fileInput = argv[1];
    ifstream infile (fileInput.c_str());
    if (infile.good() == false) {
            cout << "Errore di apertura file" << endl;
            return 1;
    }
    
    ofstream outfile ("gauss-fondo_"+fileInput);

    //Dichiarazione variabili
    string datitxt, line, txtdata, macro;
    double n1, n2, n3, n4, n5;
    double x, a, b, c, d, e, min, max;
    vector <double> vx, va, vb, vc, vd, ve;
    int N, NBin, numero = 1;
    double mu1, reb, min1, max1, livetime;
    string fileSet, file_in, inp_file, fitfun;
    const char * argv1_name;
    const char *  macro_name;
    double v_bias;
    double v_bias_err = 0.2; // errore fissato "a mano"
    double cov1, cov2, res1, res1_err, cov;
    bool notwrite;
    double o1, o2, o3, o4, o5, o6;
    string resolutions = "resolutions2.txt";
    double ming, maxg;
    double Amp, Sig, NEv, Rate;
    double Amp_err, Sig_err, NEv_err, Rate_err;
    
    double A, K, effin_2, theta, nn;
    double A_err, effin_2_err, nn_err;
    double K_err = 0.00001;
    double theta_err = 6.;
    double effin_1 = 0.196;//0.0905353;
    double effin_1_err = 0.001;//0.000132296;
    double eff_2E = 0.39;//0.180115;
    double eff_2E_err = 0.01;//0.000265647;
    double NN = 11340*(6.02214076e+23)*82/207.2;
    double ang_1 = 0.5*(1-cos(atan(0.0254/(2*0.18))));
    double ang_2 = 0.5*(1-cos(atan(0.0254/(0.26))));
    double ang_1_err = sin(atan(0.0254/(2*0.18)))*0.5*0.005*0.18*0.18*(0.0254/2)/(0.18*0.18*(0.18*0.18+pow((0.0254/2),2)));
    double ang_2_err = sin(atan(0.0254/(0.26)))*0.5*0.005*0.26*0.26*0.0254/(0.26*0.26*(0.26*0.26+pow((0.0254),2)));
    double urto, urto_err;
    
    
    
    while(true) { //lettura di fitset.txt
        infile >> datitxt >> n1 >> n2 >> n3 >> n4 >> n5;
        if (infile.eof() == true) break;
        ifstream in (datitxt.c_str());
    
        N = 0;
        NBin = 1;
        
        
        while (getline(in,line)) { //lettura del file dati
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
        
        histo = create_histo(numero, NBin-7, min, max);
        
        for(int i=0; i<vx.size(); i++) {
            for(int j =0; j<va[i];j++) histo->Fill(vx[i]);
            for(int j =0; j<vb[i];j++) histo->Fill(vx[i]+1);
            for(int j =0; j<vc[i];j++) histo->Fill(vx[i]+2);
            for(int j =0; j<vd[i];j++) histo->Fill(vx[i]+3);
            for(int j =0; j<ve[i];j++) histo->Fill(vx[i]+4);
        }
        
        fileSet = argv[1];
        ifstream fs(fileSet.c_str());
        if (fs.good() == false) {
                cout << "Errore di apertura file" << endl;
                return 1;
        }
        while(true) {
            fs >> txtdata >> mu1 >> min1 >> max1 >> reb >> livetime;
            if(fs.eof() == true) break;
            if(txtdata.compare(datitxt) == 0) break;
        }
        fs.close();
    
        fitfun = "511_Gaus+pol2";
        TF1 * f1 = new TF1 (fitfun.c_str(), gaus_pol2, min1, max1, 6);
        f1 -> SetParameter(1,mu1);
        f1 -> SetParameter(2, 150);
        //f1 -> SetParameter(3, 1620);
        //f1 -> SetParameter(4, -0.95);
        //f1 -> SetParameter(5, 1.59e-04);
        f1 -> SetParName (0, "Amp_{1}" );
        f1 -> SetParName (1, "#mu_{1}" );
        f1 -> SetParName (2, "#sigma_{1}" );
        f1 -> SetParName (3, "a_{0}" );
        f1 -> SetParName (4, "a_{1}" );
        f1 -> SetParName (5, "a_{2}" );
        f1 -> SetLineColor(kRed);
        /*if(txtdata=="030_sorg039_1.Txt") {
            cout << "****************************************" << endl;
            f1->SetParameter(0,250);
            f1->SetParameter(3,50/(2000*2000));
            f1->SetParameter(4,-8000*50/(2000*2000));
            f1->SetParameter(5,4000*4000*50/(2000*2000));
        }*/
        


        /*fitfun = "crystalball";
        TF1 * f1 = new TF1 (fitfun.c_str(), crystal, min1, max1, 6);
        f1 -> SetParameter(1,mu1);
        //f1 -> SetParameter(2, 100);
        //f1 -> SetParameter(3, 1620);
        //f1 -> SetParameter(4, -0.95);
        //f1 -> SetParameter(5, 1.59e-04);
        f1 -> SetParName (0, "Amp_{1}" );
        f1 -> SetParName (1, "#mu_{1}" );
        f1 -> SetParName (2, "#sigma_{1}" );
        f1 -> SetParName (3, "a_{0}" );
        f1 -> SetParName (4, "a_{1}" );
        f1 -> SetParName (5, "a_{2}" );
        f1 -> SetLineColor(kRed);*/
                    
        canva = create_canva(numero);
        canva->cd();
        histo->GetXaxis()->SetTitle("Channel");
        histo->GetYaxis()->SetTitle("Events");
        //histo->SetFillColor(kBlue);
        //histo->SetFillStyle(3001);

        inp_file = datitxt;
        file_in = datitxt.replace(13, 4, ".png"); //il primo numero è la posizione dell'ultimo "." ---> modificare se necessario
        macro = datitxt.replace(13, 4, ".root");
        TFile * rootf = new TFile(macro.c_str(), "RECREATE");
        
        histo->Rebin(reb);
        histo->Draw();
        /*
        if(datitxt[0]==0) {
            string fileInput2 = datitxt;
            datitxt.clear();
            for(int j=0; j<3; j++)
                datitxt[j]=fileInput2[j+1];
            v_bias = stoi(datitxt);
        }
        else {
            v_bias = stoi(datitxt.replace(4, 15, ""));
        }
        */
        v_bias = 000;
        argv1_name = file_in.c_str();
                
        TFitResultPtr r1 = histo->Fit(fitfun.c_str(), "R S", "sames");
        TMatrixDSym covariance_matrix_1 = r1 -> GetCovarianceMatrix();
        TMatrixDSym correlation_matrix_1 = r1 -> GetCorrelationMatrix();
        
        canva->Update();
        
        
        cov1 = covariance_matrix_1(1,2);
        res1 = (2.35*f1->GetParameter(2)) / (f1->GetParameter(1));
    
        res1_err = 2.35 * sqrt(pow(f1->GetParError(2)/f1->GetParameter(1),2)+pow(f1->GetParError(1)*f1->GetParameter(2)/pow(f1->GetParameter(1),2),2)+(2/f1->GetParameter(1))*(f1->GetParameter(2)/pow(f1->GetParameter(1),2))*cov1);
      
        
    /*
        Aggiungo ad un file "resolutions2.txt" delle righe così costruite:
        
        |V_bias|V_bias_err|Risoluz_picco1|Risoluz_picco1_err|Risoluz_picco2|Risoluz_picco2_err|
    */

        notwrite = false;
        
        fstream out_in (resolutions.c_str(), ios::in); //apro il file in modalità input (lettura)
        if (out_in.good() == false) {
                cout << "File 'resolutions2.txt' inesistente o corrotto" << endl;
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
        //out_app << v_bias << "\t" << v_bias_err << "\t" << res1 << "\t" << res1_err << "\t" << res2 << "\t" << res2_err << endl;
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
        
        ming = f1->GetParameter(1)-5*f1->GetParameter(2);
        maxg = f1->GetParameter(1)+5*f1->GetParameter(2);
        TF1 * fb = new TF1("fb", gaussian, ming, maxg, 3);
        Amp = f1->GetParameter(0);
        Sig = f1->GetParameter(2);
        Amp_err = f1->GetParError(0);
        Sig_err = f1->GetParError(2);
        fb -> SetParameter(0,Amp);
        fb -> SetParameter(1, f1->GetParameter(1));
        fb -> SetParameter(2, Sig);
        fb->SetLineColor(kGreen+1);
        //fb->Draw("same");
        cov = covariance_matrix_1(0,2);
        NEv = sqrt(2*M_PI) * Amp * Sig;
        NEv_err = sqrt(2*M_PI*(pow(Amp*Sig_err,2)+pow(Sig*Amp_err,2)+2*Amp*Sig*cov));
        cout << inp_file << endl;
        Rate = NEv/livetime;
        Rate_err = NEv_err/livetime;
        
        //Rate = 0.001;
        
        if(inp_file == "000_sorgBoh_1.Txt") {
            theta = 0;
            A = 66702;
            A_err = 189;
            effin_2 = eff_2E * 1.;
            effin_2_err = eff_2E_err * 1.;
            K = 0.00165;
            nn = 2*A*effin_1*ang_1;
            nn_err = 2*sqrt( pow(effin_1_err*ang_1,2)+pow(effin_1*ang_1_err,2));
        }
        else if(inp_file == "030_sorg039_1.Txt") {
            theta = 30;
            A = 66702;
            A_err = 189;
            effin_2 = eff_2E * 1.1;
            effin_2_err = eff_2E_err * 1.1;
            K = 0.00159;
            nn = 2*A*effin_1*ang_1;
            nn_err = 2*sqrt( pow(effin_1_err*ang_1,2)+pow(effin_1*ang_1_err,2));
        }
        else if(inp_file == "045_arancio_1.Txt") {
            theta = 45;
            A = 428025;
            A_err = 85605;
            effin_2 = eff_2E * 1.16;
            effin_2_err = eff_2E_err * 1.16;
            K = 0.00149;
            nn = 2*A*effin_1*ang_1;
            nn_err = 2*sqrt( pow(effin_1_err*ang_1,2)+pow(effin_1*ang_1_err,2));
        }
        else if(inp_file == "060_arancio_2.Txt") {
            theta = 60;
            A = 428025;
            A_err = 85605;
            effin_2 = eff_2E * 1.27;
            effin_2_err = eff_2E_err * 1.27;
            K = 0.00131;
            nn = 2*A*effin_1*ang_1;
            nn_err = 2*sqrt( pow(effin_1_err*ang_1,2)+pow(effin_1*ang_1_err,2));
        }
        else if(inp_file == "075_arancio_2.Txt") {
            theta = 75;
            A = 428025;
            A_err = 85605;
            effin_2 = eff_2E * 1.38;
            effin_2_err = eff_2E_err * 1.38;
            K = 0.00111;
            nn = 2*A*effin_1*ang_1;
            nn_err = 2*sqrt( pow(effin_1_err*ang_1,2)+pow(effin_1*ang_1_err,2));
        }
        else if(inp_file == "130_arancio_1.Txt") {
            theta = 130;
            A = 428025;
            A_err = 85605;
            effin_2 = eff_2E * 1.58;
            effin_2_err = eff_2E_err * 1.58;
            K = 0.00075;
            nn = 2*A*effin_1*ang_1;
            nn_err = 2*sqrt( pow(effin_1_err*ang_1,2)+pow(effin_1*ang_1_err,2));
        }
        
        //Rate= 0.001;
        //if (inp_file == "030_sorg039_1.Txt") {
            //Rate = 0.10132617;
            //Rate_err = 0.0042282;
        //}
        urto = Rate/(nn*NN*ang_2*4*M_PI*effin_2*K);
        urto_err = sqrt(pow(Rate_err/(nn*NN*ang_2*4*M_PI*effin_2*K),2)+pow(Rate*nn_err/(nn*nn*NN*ang_2*4*M_PI*effin_2*K),2)+pow(Rate*ang_2_err/(nn*NN*ang_2*ang_2*4*M_PI*effin_2*K),2)+pow(Rate*effin_2_err/(nn*NN*ang_2*4*M_PI*effin_2*effin_2*K),2)+pow(Rate*K_err/(nn*NN*ang_2*4*M_PI*effin_2*K*K),2));

        cout << "FOTONI RIVELATI = " << NEv << " +- " << NEv_err << endl;
        cout << "RATE = " << Rate << " +- " << Rate_err << endl;
        outfile << "*****************************************" << endl; 
        outfile << "FOTONI RIVELATI = " << NEv << " +- " << NEv_err << endl;
        outfile << "RATE = " << Rate << " +- " << Rate_err << endl;
        outfile << "PARAMETRI gauss" << endl;
        outfile << inp_file << endl << "Amp = " << Amp << " +- " << Amp_err << endl;
        outfile << "Mean = " << f1->GetParameter(1) << " +- " << f1->GetParError(1) << endl;
        outfile << "Sigma = " << Sig << " +- " << Sig_err << endl;
        outfile << "Sez urto = " << urto << " +- " << urto_err << endl;
        outfile << "*****************************************" << endl << endl;
        
        canva->Print(argv1_name, "png");
        macro_name = (macro).c_str();
        //canva->SaveAs(macro_name);
        histo->Write();
        rootf->Close();

               
        //svuotamento dei vector
        vx.clear();
        va.clear();
        vb.clear();
        vc.clear();
        vd.clear();
        ve.clear();
        numero++;
        
    }
    infile.close();
    myApp -> Run();
    
    return 0;

}