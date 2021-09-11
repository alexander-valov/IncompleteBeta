#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace JustMath {

/************************************************************************
 * @brief Utility structures to promote floating point types.
 ************************************************************************/
namespace detail {
    template<class T, bool = std::is_integral<T>::value>
    struct promote_fp { typedef double type; };

    template<class T>
    struct promote_fp<T, false> {};

    template<>
    struct promote_fp<long double> { typedef long double type; };

    template<>
    struct promote_fp<double> { typedef double type; };

    template<>
    struct promote_fp<float> { typedef float type; };

    template<
        class T1, class T2, 
        class TP1 = typename promote_fp<T1>::type,
        class TP2 = typename promote_fp<T2>::type
    >
    struct promote_fp_2 { typedef std::remove_reference_t<decltype(TP1() + TP2())> type; };

    template<
        class T1, class T2, class T3,
        class TP1 = typename promote_fp<T1>::type,
        class TP2 = typename promote_fp<T2>::type,
        class TP3 = typename promote_fp<T3>::type
    >
    struct promote_fp_3 { typedef std::remove_reference_t<decltype(TP1() + TP2() + TP3())> type; };

    /********************************************************************
     * @brief Implementation of beta function.
     * @param[in] a Argument a > 0
     * @param[in] b Argument b > 0
     * @return beta(a, b)
     ******************************************************************** */
    template<class T>
    T beta_impl(const T& a, const T& b) {
        return std::exp(std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
    }

    /********************************************************************
     * @brief Implementation of incomplete beta function.
     * @param[in, out] x Argument 0 <= x <= 1
     * @param[in] a Argument a > 0
     * @param[in] b Argument b > 0
     * @param[in] TOL Continued fraction approximation tolerance
     * @param[in] MAX_ITER Continued fraction approximation max iterations
     * @return Approximation of B(x, a, b)
     ******************************************************************** */
    template<class T>
    T incbeta_impl(
        const T& x,
        const T& a,
        const T& b,
        const T& TOL = T(1e-8),
        unsigned MAX_ITER = 200
    ) {
        if (x < 0 || x > 1) {
            throw std::domain_error("incbeta(x, a, b): The argument 'x' must be inside [0, 1] interval");
        }
        if (a <= 0) {
            throw std::domain_error("incbeta(x, a, b): The argument 'a' must be greater than zero");
        }
        if (b <= 0) {
            throw std::domain_error("incbeta(x, a, b): The argument 'b' must be greater than zero");
        }

        // -----------------------------------------------------------
        // Limiting cases
        // -----------------------------------------------------------
        if (x == T(0)) { return T(0); }
        if (x == T(1)) { return beta_impl(a, b); }

        // -----------------------------------------------------------
        // The continued fraction converges rapidly for
        //     x < (a + 1) / (a + b + 2)
        // In other cases is used symmetry relation
        //     B(x, a, b) = beta(a, b) - B(1 - x, b, a)
        // -----------------------------------------------------------
        if (x - (a + 1) / (a + b + 2) > std::numeric_limits<T>::epsilon()) {
            return (beta_impl(a, b) - incbeta_impl(1 - x, b, a, TOL, MAX_ITER));
        }

        // -----------------------------------------------------------
        // Calculate the front part of incomplete beta function
        // -----------------------------------------------------------
        T front = std::exp(a * std::log(x) + b * std::log(1 - x)) / a;

        // -----------------------------------------------------------
        // Calculate continued fraction part using Lentz's algorithm
        // -----------------------------------------------------------
        T D = 0;
        T C = 1;
        T F = 1;
        T CD = 0;
        T THRESHOLD = T(1e-30);
        for (unsigned i = 0; i < MAX_ITER; i++) {
            // -----------------------------------------------------------
            // Calculate numerator of continued fraction
            // -----------------------------------------------------------
            T numerator;
            unsigned m = i / 2;
            if (i == 0) {
                numerator = 1;
            } else if (i % 2 == 0) {
                numerator = (m * (b - m) * x) / ((a + 2 * m - 1) * (a + 2 * m));
            } else {
                numerator = -((a + m) * (a + b + m) * x) / ((a + 2 * m) * (a + 2 * m + 1));
            }

            // -----------------------------------------------------------
            // Iteration of Lentz's algorithm
            // -----------------------------------------------------------
            D = 1 + numerator * D;
            if (std::abs(D) < THRESHOLD) { D = THRESHOLD; }
            D = 1 / D;

            C = 1 + numerator / C;
            if (std::abs(C) < THRESHOLD) { C = THRESHOLD; }

            CD = C * D;
            F *= CD;

            // -----------------------------------------------------------
            // Check for convergence
            // -----------------------------------------------------------
            if (std::abs(1 - CD) < TOL) {
                return front * (F - 1);
            }
        }
        throw std::logic_error("incbeta(x, a, b): no convergence for given TOL and MAX_ITER");
    }
} // end of detail namespace

/********************************************************************
 * @brief Beta function.
 * 
 * @see https://en.wikipedia.org/wiki/Beta_function
 * 
 * @param[in] a Argument a > 0
 * @param[in] b Argument b > 0
 * @return beta(a, b)
 ******************************************************************** */
template<class Ta, class Tb>
typename detail::promote_fp_2<Ta, Tb>::type beta(
    const Ta& a,
    const Tb& b
) {
    typedef typename detail::promote_fp_2<Ta, Tb>::type type;
    return detail::beta_impl<type>(a, b);
}

/********************************************************************
 * @brief Incomplete beta function.
 * 
 * Incomplete beta function symmetry:
 *     B(x, a, b) = beta(a, b) - B(1 - x, b, a)
 * 
 * @see https://codeplea.com/incomplete-beta-function-c
 * @see https://dlmf.nist.gov/8.17#SS5.p1
 * 
 * @param[in, out] x Argument 0 <= x <= 1
 * @param[in] a Argument a > 0
 * @param[in] b Argument b > 0
 * @param[in] TOL Continued fraction approximation tolerance
 * @param[in] MAX_ITER Continued fraction approximation max iterations
 * @return Approximation of B(x, a, b)
 ******************************************************************** */
template<class Tx, class Ta, class Tb>
typename detail::promote_fp_3<Tx, Ta, Tb>::type incbeta(
    const Tx& x,
    const Ta& a,
    const Tb& b,
    const typename detail::promote_fp_3<Tx, Ta, Tb>::type& TOL = 
        typename detail::promote_fp_3<Tx, Ta, Tb>::type(1e-8),
    unsigned MAX_ITER = 200
) {
    typedef typename detail::promote_fp_3<Tx, Ta, Tb>::type type;
    return detail::incbeta_impl<type>(x, a, b, TOL, MAX_ITER);
}

/********************************************************************
 * @brief Regularized incomplete beta function:
 *     I(x, a, b) = B(x, a, b) / B(a, b)
 * @param[in, out] x Argument 0 <= x <= 1
 * @param[in] a Argument a > 0
 * @param[in] b Argument b > 0
 * @param[in] TOL Continued fraction approximation tolerance
 * @param[in] MAX_ITER Continued fraction approximation max iterations
 * @return Approximation of I(x, a, b)
 ******************************************************************** */
template<class Tx, class Ta, class Tb>
typename detail::promote_fp_3<Tx, Ta, Tb>::type incbeta_reg(
    const Tx& x,
    const Ta& a,
    const Tb& b,
    const typename detail::promote_fp_3<Tx, Ta, Tb>::type& TOL = 
        typename detail::promote_fp_3<Tx, Ta, Tb>::type(1e-8),
    unsigned MAX_ITER = 200
) {
    return incbeta(x, a, b, TOL, MAX_ITER) / beta(a, b);
}

} // end of JustMath namespace