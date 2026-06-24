#include <cmath>

#include "interpolation.h"
#include "tests/interpolation_reference.h"
#include "tests/test_helpers.h"

namespace {

constexpr double k_field = 3.0e6;  // V/cm
constexpr double k_temp = 298.0;   // K
constexpr double k_tol = 1e-3;
constexpr double k_gamma_tol = 1.0;

void test_alpha_endpoints() {
    CHECK_NEAR(alpha(0.0, k_field, k_temp),
               reference::alpha_gan_endpoint(k_field, k_temp), k_tol);
    CHECK_NEAR(alpha(0.65, k_field, k_temp),
               reference::alpha_algan_endpoint(k_field), k_tol);
}

void test_beta_endpoints() {
    CHECK_NEAR(beta(0.0, k_field, k_temp),
               reference::beta_gan_endpoint(k_field, k_temp), k_tol);
    CHECK_NEAR(beta(0.65, k_field, k_temp),
               reference::beta_algan_endpoint(k_field), k_tol);
}

void test_alpha_midpoint_interpolation() {
    const double al_frac = 0.325;
    const double expected = reference::alpha(al_frac, k_field, k_temp);
    CHECK_NEAR(alpha(al_frac, k_field, k_temp), expected, k_tol);
}

void test_beta_midpoint_interpolation() {
    const double al_frac = 0.325;
    const double expected = reference::beta(al_frac, k_field, k_temp);
    CHECK_NEAR(beta(al_frac, k_field, k_temp), expected, k_tol);
}

void test_alpha_decreases_with_al_fraction() {
    const double low = alpha(0.1, k_field, k_temp);
    const double high = alpha(0.5, k_field, k_temp);
    CHECK_TRUE(low > high);
}

void test_beta_is_positive_for_typical_inputs() {
    const double value = beta(0.3, k_field, k_temp);
    CHECK_TRUE(value > 0.0);
    CHECK_TRUE(std::isfinite(value));
}

void test_alpha_and_beta_differ() {
    const double a = alpha(0.3, k_field, k_temp);
    const double b = beta(0.3, k_field, k_temp);
    CHECK_TRUE(std::abs(a - b) > 1.0);
}

void test_temperature_dependence() {
    const double cold = beta(0.2, k_field, 250.0);
    const double hot = beta(0.2, k_field, 350.0);
    CHECK_TRUE(std::abs(cold - hot) > 1.0);
}

void test_gamma_at_sampled_energy() {
    const double al_frac = 0.0;
    const double energy = 3.125;
    const double expected = reference::gamma(al_frac, energy);
    CHECK_NEAR(gamma(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_energy_interpolation_at_gan() {
    const double al_frac = 0.0;
    const double energy = 3.0;
    const double expected = reference::gamma(al_frac, energy);
    CHECK_NEAR(gamma(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_mole_fraction_interpolation() {
    const double al_frac = 0.055;
    const double energy = 3.0;
    const double expected = reference::gamma(al_frac, energy);
    CHECK_NEAR(gamma(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_asymptotic_high_energy() {
    const double al_frac = 0.0;
    const double energy = 4.0;
    const double expected = std::pow(10.0, 4.85);
    CHECK_NEAR(gamma(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_extrapolation_below_lowest_energy() {
    const double al_frac = 0.0;
    const double energy = 2.0;
    const double expected = reference::gamma(al_frac, energy);
    CHECK_NEAR(gamma(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_extrapolation_above_highest_mole_frac() {
    const double al_frac = 0.9;
    const double energy = 4.5;
    const double expected = reference::gamma(al_frac, energy);
    CHECK_NEAR(gamma(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_mid_interval() {
    const double al_frac = 0.29;
    const double energy = 3.6;
    const double expected = reference::gamma(al_frac, energy);
    CHECK_NEAR(gamma(al_frac, energy), expected, k_gamma_tol);
}

void test_gamma_is_positive() {
    const double value = gamma(0.38, 3.75);
    CHECK_TRUE(value > 0.0);
    CHECK_TRUE(std::isfinite(value));
}

}  // namespace

int main() {
    test_alpha_endpoints();
    test_beta_endpoints();
    test_alpha_midpoint_interpolation();
    test_beta_midpoint_interpolation();
    test_alpha_decreases_with_al_fraction();
    test_beta_is_positive_for_typical_inputs();
    test_alpha_and_beta_differ();
    test_temperature_dependence();
    test_gamma_at_sampled_energy();
    test_gamma_energy_interpolation_at_gan();
    test_gamma_mole_fraction_interpolation();
    test_gamma_asymptotic_high_energy();
    test_gamma_extrapolation_below_lowest_energy();
    test_gamma_extrapolation_above_highest_mole_frac();
    test_gamma_mid_interval();
    test_gamma_is_positive();

    std::cout << test_count() << " checks, " << test_failures() << " failures\n";
    return test_failures() == 0 ? 0 : 1;
}
