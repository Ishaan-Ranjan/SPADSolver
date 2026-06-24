#pragma once

#include <cmath>
#include <iostream>
#include <string>

inline int& test_failures() {
    static int failures = 0;
    return failures;
}

inline int& test_count() {
    static int count = 0;
    return count;
}

inline void check_true(bool condition, const std::string& message) {
    ++test_count();
    if (!condition) {
        ++test_failures();
        std::cerr << "FAIL: " << message << '\n';
    }
}

inline void check_near(double actual, double expected, double tolerance,
                       const std::string& message) {
    ++test_count();
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        ++test_failures();
        std::cerr << "FAIL: " << message << " (got " << actual << ", expected "
                  << expected << ", tolerance " << tolerance << ")\n";
    }
}

#define CHECK_TRUE(cond) check_true((cond), #cond)
#define CHECK_NEAR(actual, expected, tol) \
    check_near((actual), (expected), (tol), #actual " == " #expected)
