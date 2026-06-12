#include <iostream>
#include <cmath>
#include <vector> 
#include "interpolation.h"

double lininterpolate(double x, double x1, double y1, double x2, double y2){return y1 + (x - x1) * (y2 - y1) / (x2 - x1);}

//IONIZATION COEFFICIENTS

double GaNionizationConsts(double front, double top, double Field){ return front * std::exp(-top/Field);}

double frontGaN(double a, double b, double Temp){return a * (1 + b * (Temp - 298));}

double topGaN(double c, double d, double Temp){return -1 * c * (1 + d * (Temp - 298));}

double AlGaNConsts(double a, double b, double Field){ return a * std::exp(-b/Field);}

double alpha(double alFrac, double Field, double Temp){

    //GaN electron constsants

    double a_n =  2.69e7; // cm^-1
    double b_n = 2.00e-3; // K^-1
    double c_n = 2.27e7; // V * cm^-1
    double d_n = 5.00e-4; // K^-1

    double frontGaN_n = topGaN(a_n, b_n, Temp);
    double topGaN_n = frontGaN(c_n, d_n, Temp);

    //0.65 AlGaN electron constants

    double frontAlpha = 7.82e6; // cm^-1
    double topAlpha = 3.7e7; // V * cm^-1

    return lininterpolate(alFrac, 0, GaNionizationConsts(frontGaN_n, topGaN_n, Field), 0.65, AlGaNConsts(frontAlpha, topAlpha, Field));
}

double beta(double alFrac, double Field, double Temp){

    //GaN holes constsants
    double a_p =  4.32e6; // cm^-1
    double b_p = 2.00e-3; // K^-1
    double c_p = 1.31e7; // V * cm^-1
    double d_p = 9.00e-4; // K^-1
    
    double frontGaN_p = frontGaN(a_p, b_p, Temp);

    double topGaN_p = topGaN(c_p, d_p, Temp);

    //0.65 AlGaN hole constants

    double frontBeta = 5.65e4; // cm^-1
    double topBeta = 7.04e6; // V * cm^-1

    return lininterpolate(alFrac, 0, GaNionizationConsts(frontGaN_p, topGaN_p, Field), 0.65, AlGaNConsts(frontBeta, topBeta, Field));
}

//ABSORPTIONS

double gamma(double AlFrac, double Energy){
    double mole_fracs[5] = {0, 0.27, 0.34, 0.38, 1};

    double sampled_energies[4][2] = {{4.1, 4.7}, {4.5, 5.1}, {4.6, 5.1}, {4.7, 5.3}};

    double gamma1;
    double gamma2;

    double sample_abs1 = 14e4; // per cm
    double sample_abs2 = 19e4; // per cm

    for(int i = 0; i < 5; i++){
        //This section of the code interpolates absorption values between mole fractions
        if ((mole_fracs[i] <= AlFrac <= mole_fracs[i+1]) && (AlFrac < 0.38)){
            gamma1 = lininterpolate(Energy, sampled_energies[i][0], sample_abs1, sampled_energies[i][1], sample_abs2);
            gamma2 = lininterpolate(Energy, sampled_energies[i+1][0], sample_abs1, sampled_energies[i+1][1], sample_abs2);
            return lininterpolate(AlFrac, mole_fracs[i], gamma1, mole_fracs[i+1], gamma2);
        }
        else if ((0.38 <= AlFrac) && (i == 3)){
        //This section of the code extrapolates to AlN since absorptions past its bandgap are unknown from the paper
            gamma1 = lininterpolate(Energy, sampled_energies[i-1][0], sample_abs1, sampled_energies[i-1][1], sample_abs2);
            gamma2 = lininterpolate(Energy, sampled_energies[i][0], sample_abs1, sampled_energies[i][1], sample_abs2);
            return lininterpolate(AlFrac, mole_fracs[i-1], gamma1, mole_fracs[i], gamma2); 
        }
        return 0;
    }
    return 0;
}