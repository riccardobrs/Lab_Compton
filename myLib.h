#ifndef myLib_h
#define myLib_h

using namespace std;

double rand_range(double min, double max);
double f_Cauchy(double x);
double rand_TAC(double xMin, double xMax,double yMin,double yMax, double pdf(double x));
double f_Exp (double x, double parameter);
double rand_FCI_Exp(double parameter);
double f_gauss(double x, double mean, double sigma);
double rand_TAC_gauss(double xMin, double xMax,double yMin,double yMax,double mean, double sigma);
bool HitMiss_gauss(double xMin, double xMax,double yMin,double yMax, double mean, double sigma);
double Crude_MC_Gauss (double min, double max,double mean, double sigma);

void Max_Min(std::vector <double> vec,double&  min,double&  max);
void Media_DevStd(std::vector <double> x,double& media,double& devStd);
void Covarianza(std::vector <double> x,std::vector<double> y,double& covarianza);
double Gauss(double* x, double* par);
double Binormale(double* x,double* par);
void Media_Pesata(std::vector <double> x,std::vector <double> err_x, double& media,double& errore);
double gaus_pol2(double* x, double* par);
vector <double> ris (string filename, double range []); //la funzione restituisce |Ris Picco1|Err Ris1|Ris Picco2|Err Ris2|

#endif

