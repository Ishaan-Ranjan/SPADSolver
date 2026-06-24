#pragma once

#include <cmath>

namespace reference {

inline double lininterpolate(double x, double x1, double y1, double x2, double y2) {
    return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

inline double gan_ionization(double front, double top, double field) {
    return front * std::exp(-top / field);
}

inline double algan_ionization(double a, double b, double field) {
    return a * std::exp(-b / field);
}

inline double alpha_gan_endpoint(double field, double temp) {
    double a_n = 2.69e7;
    double b_n = 2.00e-3;
    double c_n = 2.27e7;
    double d_n = 5.00e-4;

    double front = a_n * (1.0 + b_n * (temp - 298.0));
    double top = -1.0 * c_n * (1.0 + d_n * (temp - 298.0));
    return gan_ionization(front, top, field);
}

inline double alpha_algan_endpoint(double field) {
    return algan_ionization(7.82e6, 3.7e7, field);
}

inline double alpha(double al_frac, double field, double temp) {
    return lininterpolate(al_frac, 0.0, alpha_gan_endpoint(field, temp), 0.65,
                          alpha_algan_endpoint(field));
}

inline double beta_gan_endpoint(double field, double temp) {
    double a_p = 4.32e6;
    double b_p = 2.00e-3;
    double c_p = 1.31e7;
    double d_p = 9.00e-4;

    double front = a_p * (1.0 + b_p * (temp - 298.0));
    double top = -1.0 * c_p * (1.0 + d_p * (temp - 298.0));
    return gan_ionization(front, top, field);
}

inline double beta_algan_endpoint(double field) {
    return algan_ionization(5.65e4, 7.04e6, field);
}

inline double beta(double al_frac, double field, double temp) {
    return lininterpolate(al_frac, 0.0, beta_gan_endpoint(field, temp), 0.65,
                          beta_algan_endpoint(field));
}

inline double gamma_exp_at_index(int index, double energy) {
    const double sampled_energies[6][4] = {
        {2.25, 3.125, 3.25, 3.5},   {2.625, 3.25, 3.5, 3.7}, {2.25, 3, 3.5, 3.875},
        {2.5, 3.5, 3.75, 4.125},   {2.25, 3.5, 4, 4.625},    {2.875, 4, 4.75, 5.5}};
    const double sampled_abs_asy[6] = {4.85, 4.85, 4.8, 4.8, 4.8, 4.9};
    const double sampled_abs_exps[6][4] = {
        {1.85, 2.4, 2.7, 4.85}, {2.4, 2.84, 3.17, 4.85}, {2, 2.78, 2.95, 4.8},
        {2.3, 2.95, 3.3, 4.8},  {1.95, 2.48, 2.9, 4.8},  {2.48, 2.85, 3.3, 4.9}};

    double gamma_exp = 0.0;
    for (int j = 0; j < 3; ++j) {
        if (energy <= sampled_energies[index][0]) {
            gamma_exp = lininterpolate(energy, sampled_energies[index][j],
                                       sampled_abs_exps[index][j],
                                       sampled_energies[index][j + 1],
                                       sampled_abs_exps[index][j + 1]);
            break;
        }
        if (sampled_energies[index][j] <= energy &&
            energy <= sampled_energies[index][j + 1]) {
            gamma_exp = lininterpolate(energy, sampled_energies[index][j],
                                       sampled_abs_exps[index][j],
                                       sampled_energies[index][j + 1],
                                       sampled_abs_exps[index][j + 1]);
            break;
        }
        gamma_exp = sampled_abs_asy[index];
    }
    return gamma_exp;
}

inline double gamma(double al_frac, double energy) {
    const double mole_fracs[6] = {0, 0.11, 0.20, 0.38, 0.5, 0.86};

    double gamma_final = 0.0;
    for (int i = 0; i < 5; ++i) {
        double mol_low = 0.0;
        double mol_high = 0.0;

        if (mole_fracs[i] <= al_frac && al_frac <= mole_fracs[i + 1]) {
            mol_low = mole_fracs[i];
            mol_high = mole_fracs[i + 1];
        } else if (i == 4) {
            mol_low = mole_fracs[4];
            mol_high = mole_fracs[5];
        } else {
            continue;
        }

        const double gamma_exp1 = gamma_exp_at_index(i, energy);
        const double gamma_exp2 = gamma_exp_at_index(i + 1, energy);
        gamma_final =
            lininterpolate(al_frac, mol_low, gamma_exp1, mol_high, gamma_exp2);
        gamma_final = std::pow(10.0, gamma_final);
        break;
    }
    return gamma_final;
}

}  // namespace reference
