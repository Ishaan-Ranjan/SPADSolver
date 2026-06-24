#include <cmath>
#include <vector>

#include "GaussianElimination.h"
#include "tests/test_helpers.h"

namespace {

constexpr double k_tol = 1e-9;

std::vector<double> make_augmented_matrix(
    const std::vector<std::vector<double>>& rows, int bins) {
    std::vector<double> matrix(bins * (bins + 1));
    for (int i = 0; i < bins; ++i) {
        for (int j = 0; j <= bins; ++j) {
            matrix[i * (bins + 1) + j] = rows[i][j];
        }
    }
    return matrix;
}

std::vector<double> solve_system(std::vector<double> matrix, int bins) {
    rowEchelon(matrix, bins);
    return solVector(matrix, bins);
}

bool solution_satisfies_system(const std::vector<std::vector<double>>& coeffs,
                               const std::vector<double>& rhs,
                               const std::vector<double>& sol, double tol) {
    const int bins = static_cast<int>(rhs.size());
    for (int i = 0; i < bins; ++i) {
        double residual = -rhs[i];
        for (int j = 0; j < bins; ++j) {
            residual += coeffs[i][j] * sol[j];
        }
        if (std::abs(residual) > tol) {
            return false;
        }
    }
    return true;
}

void check_solution(const std::vector<std::vector<double>>& coeffs,
                    const std::vector<double>& rhs,
                    const std::vector<double>& expected) {
    const int bins = static_cast<int>(rhs.size());
    std::vector<std::vector<double>> rows(bins);
    for (int i = 0; i < bins; ++i) {
        rows[i] = coeffs[i];
        rows[i].push_back(rhs[i]);
    }

    const std::vector<double> sol = solve_system(make_augmented_matrix(rows, bins), bins);

    for (int i = 0; i < bins; ++i) {
        CHECK_NEAR(sol[i], expected[i], k_tol);
    }
    CHECK_TRUE(solution_satisfies_system(coeffs, rhs, sol, k_tol));
}

void test_1x1_system() {
    check_solution({{5.0}}, {10.0}, {2.0});
}

void test_2x2_system() {
    // 2x + y = 5
    // x + 3y = 9
    check_solution({{2.0, 1.0}, {1.0, 3.0}}, {5.0, 9.0}, {1.2, 2.6});
}

void test_2x2_diagonal_system() {
    check_solution({{4.0, 0.0}, {0.0, 5.0}}, {8.0, 15.0}, {2.0, 3.0});
}

void test_3x3_banded_system() {
    // 4x + z = 9
    // 2x + 2y = 6
    // 3y + z = 7
    check_solution({{4.0, 0.0, 1.0}, {2.0, 2.0, 0.0}, {0.0, 3.0, 1.0}},
                   {9.0, 6.0, 7.0},
                   {11.0 / 7.0, 10.0 / 7.0, 19.0 / 7.0});
}

void test_sol_vector_extracts_last_column() {
    const int bins = 2;
    std::vector<double> matrix = {1.0, 0.0, 3.0, 0.0, 1.0, 5.0};

    const std::vector<double> sol = solVector(matrix, bins);

    CHECK_NEAR(sol[0], 3.0, k_tol);
    CHECK_NEAR(sol[1], 5.0, k_tol);
}

void test_row_echelon_reduced_form_2x2() {
    const int bins = 2;
    std::vector<double> matrix = {2.0, 1.0, 5.0, 1.0, 3.0, 9.0};

    rowEchelon(matrix, bins);

    CHECK_NEAR(matrix[0 * 3 + 0], 1.0, k_tol);
    CHECK_NEAR(matrix[0 * 3 + 1], 0.0, k_tol);
    CHECK_NEAR(matrix[0 * 3 + 2], 1.2, k_tol);
    CHECK_NEAR(matrix[1 * 3 + 0], 0.0, k_tol);
    CHECK_NEAR(matrix[1 * 3 + 1], 1.0, k_tol);
    CHECK_NEAR(matrix[1 * 3 + 2], 2.6, k_tol);
}

void test_row_echelon_reduced_form_3x3() {
    const int bins = 3;
    std::vector<double> matrix = {4.0, 0.0, 1.0, 9.0, 2.0, 2.0, 0.0, 6.0,
                                  0.0, 3.0, 1.0, 7.0};
    const std::vector<double> expected = {11.0 / 7.0, 10.0 / 7.0, 19.0 / 7.0};

    rowEchelon(matrix, bins);

    for (int i = 0; i < bins; ++i) {
        CHECK_NEAR(matrix[i * 4 + i], 1.0, k_tol);
        CHECK_NEAR(matrix[i * 4 + 3], expected[i], k_tol);
    }
}

}  // namespace

int main() {
    test_1x1_system();
    test_2x2_system();
    test_2x2_diagonal_system();
    test_3x3_banded_system();
    test_sol_vector_extracts_last_column();
    test_row_echelon_reduced_form_2x2();
    test_row_echelon_reduced_form_3x3();

    std::cout << test_count() << " checks, " << test_failures() << " failures\n";
    return test_failures() == 0 ? 0 : 1;
}
