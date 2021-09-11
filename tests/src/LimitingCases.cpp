#include <cmath>

#include "doctest.h"
#include "JustMath/Beta.hpp"

TEST_SUITE("LimitingCases") {

TEST_CASE("LimitingCases.x=0") {
    /* double */
    double a_d = 10.0;
    double b_d = 10.0;
    double res_d = 0;
    CHECK_EQ(JustMath::incbeta(0.0, a_d, b_d), res_d);

    /* float */
    float a_f = 10.0;
    float b_f = 10.0;
    float res_f = 0;
    CHECK_EQ(JustMath::incbeta(0.0f, a_f, b_f), res_f);

    /* int */
    CHECK_EQ(JustMath::incbeta(0, 2, 5), 0);
}

TEST_CASE("LimitingCases.x=1") {
    /* double */
    double a_d = 0.5;
    double b_d = 0.7;
    double res_d = JustMath::beta(a_d, b_d);
    CHECK_EQ(JustMath::incbeta(1.0, a_d, b_d), res_d);

    /* float */
    float a_f = 0.5;
    float b_f = 0.7;
    float res_f = JustMath::beta(a_f, b_f);;
    CHECK_EQ(JustMath::incbeta(1.0f, a_f, b_f), res_f);

    /* float-double-double */
    CHECK_EQ(JustMath::incbeta(1.0f, a_d, b_d), res_d);
}

TEST_CASE("LimitingCases.a=1") {
    int n_x = 1000;
    int n_b = 10000;
    double db = 0.01;
    double dx = 1.0 / n_x;

    /* fine mesh along x */
    for (int ix = 1; ix < n_x; ix++) {
        double x = ix * dx;
        for (int ib = 1; ib < n_b; ib++) {
            double b = ib * db;
            double exact_res = (1 - std::pow(1 - x, b)) / b;
            CHECK(JustMath::incbeta(x, 1, b) == doctest::Approx(exact_res));
        }
    }

    /* large values along b */
    n_x /= 20;
    n_b *= 20;
    db *= 1000;
    dx = 1.0 / n_x;
    for (int ix = 1; ix < n_x; ix++) {
        double x = ix * dx;
        for (int ib = 1; ib < n_b; ib++) {
            double b = ib * db;
            double exact_res = (1 - std::pow(1 - x, b)) / b;
            CHECK(JustMath::incbeta(x, 1, b) == doctest::Approx(exact_res));
        }
    }
}

TEST_CASE("LimitingCases.b=1") {
    int n_x = 1000;
    int n_a = 10000;
    double da = 0.01;
    double dx = 1.0 / n_x;

    /* fine mesh along x */
    for (int ix = 1; ix < n_x; ix++) {
        double x = ix * dx;
        for (int ia = 1; ia < n_a; ia++) {
            double a = ia * da;
            double exact_res = std::pow(x, a) / a;
            CHECK(JustMath::incbeta(x, a, 1) == doctest::Approx(exact_res));
        }
    }

    /* large values along a */
    n_x /= 20;
    n_a *= 20;
    da *= 1000;
    dx = 1.0 / n_x;
    for (int ix = 1; ix < n_x; ix++) {
        double x = ix * dx;
        for (int ia = 1; ia < n_a; ia++) {
            double a = ia * da;
            double exact_res = std::pow(x, a) / a;
            CHECK(JustMath::incbeta(x, a, 1) == doctest::Approx(exact_res));
        }
    }
}

}