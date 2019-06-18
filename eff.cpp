#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>

using namespace std;

int main (int argc, char ** argv) {
    
    double A = 6.67E+04;
    double A_err = 100;
    double BR = 0.905;
    double r_coinc = 5.20829;
    double r_coinc_err = 0.00666156;
    double r_1 = 28.9165;
    double r_1_err = 0.145727;
    double r_2 = 57.5277;
    double r_2_err = 0.201624;
    double eff_1 = r_coinc/r_2;
    double eff_2 = r_coinc/r_1;
    double eff_1_err = sqrt(pow((r_coinc_err/r_2),2)+pow((r_coinc*r_2_err*r_2_err/(r_2*r_2)),2));
    double eff_2_err = sqrt(pow((r_coinc_err/r_1),2)+pow((r_coinc*r_1_err*r_1_err/(r_1*r_1)),2));
    double num = r_1*r_2;
    double num_err = sqrt( pow((r_1*r_2_err),2)+pow((r_2*r_1_err),2));
    double den = 2*A*BR*r_coinc;
    double den_err = 2*BR*sqrt(pow((A*r_coinc_err),2)+pow((r_coinc*A_err),2));
    double ang_solido_1 = num/den;
    double ang_solido_1_err = sqrt(pow((num_err/den),2)+pow((num*den_err/(den*den)),2));
    
    double eff_ass_1 = eff_1*ang_solido_1;
    double eff_ass_1_err = sqrt(pow((eff_1*ang_solido_1_err),2)+pow((eff_1_err*ang_solido_1),2));
    double eff_ass_2 = eff_2*ang_solido_1;
    double eff_ass_2_err = sqrt(pow((eff_2*ang_solido_1_err),2)+pow((eff_2_err*ang_solido_1),2));
    
    cout << "(Angolo solido 1pol / 4pi) = " << ang_solido_1 << " +- " << ang_solido_1_err << endl;
    cout << "Efficienza assoluta riv 1 18cm = " << eff_ass_1 << " +- " << eff_ass_1_err << endl;
    cout << "Efficienza assoluta riv 2 36cm = " << eff_ass_2 << " +- " << eff_ass_2_err << endl;
    
    return 0;

}