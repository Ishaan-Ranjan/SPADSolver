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

    double frontGaN_n = frontGaN(a_n, b_n, Temp);
    double topGaN_n = topGaN(c_n, d_n, Temp);

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
    //We are going to a log based interpolation, meaning we interpolate to find an exponent, and use ten as the base
    const double mole_fracs[6] = {0, 0.11, 0.20, 0.38, 0.5, 0.86};

    //Graphs look like they can be separated into 3 lines, so we will need 4 sampled energies/abs coeff for each mole_frac data point
    const double sampled_energies[6][4] = {{2.25,3.125,3.25,3.5},{2.625,3.25,3.5,3.7},{2.25,3,3.5,3.875},{2.5,3.5,3.75,4.125},{2.25,3.5,4,4.625},{2.875,4,4.75,5.5}};
    const double sampled_abs_asy[6] = {4.85,4.85,4.8,4.8,4.8,4.9}; //The absorbtion coeff level off above the last energy sample
    const double sampled_abs_exps[6][4] = {{1.85,2.4,2.7,4.85},{2.4,2.84,3.17,4.85},{2,2.78,2.95,4.8},{2.3,2.95,3.3,4.8},{1.95,2.48,2.9,4.8},{2.48,2.85,3.3,4.9}};

    double gamma_exp1;
    double gamma_exp2;

    double molLow;
    double molHigh;

    double gamma_final = 0;

    for(int i = 0; i < 5; i++){ //Loops through mole fraction intervals until the correct one is found
        if (mole_fracs[i] <= AlFrac && AlFrac <= mole_fracs[i+1]){
            molLow = mole_fracs[i];
            molHigh = mole_fracs[i+1];
        } else if (i == 4){ //If the mole fraction is outside the possible intervals, extrapolate using the last two points
            molLow = mole_fracs[4];
            molHigh = mole_fracs[5];
        } else{
            continue;
        }
        //Now that we have the correct mole fraction interval, we need to find the correct energy and absorption coeff exponents
        //Also, we can now interpolate and get the lower and upper bound absorption coefficients

        //Doesn't account for above energy intervals yet
        for (int j = 0; j<3; j++){//For the lower bound
            if (Energy <= sampled_energies[i][0]){//If energy is lower than all intervals, extrapolate
                gamma_exp1 = lininterpolate(Energy, sampled_energies[i][j], sampled_abs_exps[i][j], sampled_energies[i][j+1], sampled_abs_exps[i][j+1]);
                break;
            }

            if (sampled_energies[i][j] <= Energy && Energy <= sampled_energies[i][j+1]){ 
                gamma_exp1 = lininterpolate(Energy, sampled_energies[i][j], sampled_abs_exps[i][j], sampled_energies[i][j+1], sampled_abs_exps[i][j+1]);
                break;
            } 
            gamma_exp1 = sampled_abs_asy[i]; //If the energy is above some threshold, set it to the asymptotic value
        } 

        for (int k = 0; k<3; k++){//For the upper bound
            if (Energy <= sampled_energies[i+1][0]){//If energy is lower than all intervals, extrapolate
                gamma_exp2 = lininterpolate(Energy, sampled_energies[i+1][k], sampled_abs_exps[i+1][k], sampled_energies[i+1][k+1], sampled_abs_exps[i+1][k+1]);
                break;
            }

            if (sampled_energies[i+1][k] <= Energy && Energy <= sampled_energies[i+1][k+1]){  
                gamma_exp2 = lininterpolate(Energy, sampled_energies[i+1][k], sampled_abs_exps[i+1][k], sampled_energies[i+1][k+1], sampled_abs_exps[i+1][k+1]);
                break;
            }
            gamma_exp2 = sampled_abs_asy[i+1]; //If energy is above a threshold
        }

        gamma_final = lininterpolate(AlFrac, molLow, gamma_exp1, molHigh, gamma_exp2); //Calculates the final exponent
        gamma_final = std::pow(10, gamma_final); //Calculates gamma_final
        break;
    }
    return gamma_final;
}