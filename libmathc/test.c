#include "test.h"
#include "combinatorics.h"
#include "quadratics.h"
#include "integral.h"
#include "trigonometry.h"
#include "exponential.h"
#include "numerical.h"
#include "polynomial.h"
#include "complex.h"

int tests_passed = 0;
int tests_failed = 0;

float test_f_quad (float x) {return x * x - 4;}
float test_df_quad (float x) {return 2 * x;}
float test_f_cubic (float x) {return x * x * x - 27;}
float test_df_cubic (float x) {return 3 * x * x;}
float test_f_parabola (float x) {return x * x;}
float test_f_line (float x) {return 2 * x - 6;}

int main() {
    printf("=== Combinatorics ===\n");
    CHECK_INT(factorial(0),1);
    CHECK_INT(factorial(1),1);
    CHECK_INT(factorial(5),120);
    CHECK_INT(factorial(-1),0);
    CHECK_INT(npr(5,2),20);
    CHECK_INT(npr(5,0),1);
    CHECK_INT(npr(5,6),0);
    CHECK_INT(ncr(5,2),10);
    CHECK_INT(ncr(5,0),1);
    CHECK_INT(ncr(5,5),1);

    printf("\n=== Quadratics ===\n");
    CHECK_NEAR(equation(1,2,1,3),16);
    CHECK_NEAR(formula_pos(1,3,2),-1);
    CHECK_NEAR(formula_neg(1,3,2),-2);

    printf("\n=== Integral ===\n");
    CHECK_NEAR_LOOSE(factorial_float(5),120,5);
    CHECK_NEAR(integral_traditional(0,1,2,0,0.001),0.333333);
    CHECK_NEAR(gamma_new(1),1);
    CHECK_NEAR(gamma_new(2),1);
    CHECK_NEAR(gamma_new(3),2);
    CHECK_NEAR(gamma_new(0.5),1.772454);

    printf("\n=== Trigonometry ===\n");
    CHECK_NEAR(sine(0),0);
    CHECK_NEAR(sine(math_pi/2),1);
    CHECK_NEAR(sine(math_pi),0);
    CHECK_NEAR(cosine(0),1);
    CHECK_NEAR(cosine(math_pi/2),0);
    CHECK_NEAR(cosine(math_pi),-1);
    CHECK_NEAR(tangent(0),0);
    CHECK_NEAR(tangent(math_pi/4),1);
    CHECK_NEAR(arcsine(0),0);
    CHECK_NEAR(arcsine(1),math_pi/2);
    CHECK_NEAR(arccosine(0),math_pi/2);
    CHECK_NEAR(arccosine(1),0);
    CHECK_NEAR(arctangent(0),0);
    CHECK_NEAR(arctangent(1),math_pi/4);
    CHECK_NEAR(arctangent(100),1.560797);
    CHECK_NEAR(arccotangent(1),math_pi/4);
    CHECK_NEAR(arccotangent(-1),3*math_pi/4);
    CHECK_NEAR(arccotangent(0),math_pi/2);

    printf("\n=== Exponential ===\n");
    CHECK_NEAR(exponential(0),1);
    CHECK_NEAR(exponential(1),2.718282);
    CHECK_NEAR(logarithm(1),0);
    CHECK_NEAR(logarithm(math_e),1);
    CHECK_NEAR(power(2,3),8);
    CHECK_NEAR(power(9,0.5),3);
    CHECK_NEAR(logarithm_base(8,2),3);
    CHECK_NEAR(hyperbolic_sine(0),0);
    CHECK_NEAR(hyperbolic_cosine(0),1);
    CHECK_NEAR(hyperbolic_tangent(0),0);
    CHECK_NEAR(inverse_hyperbolic_sine(0),0);
    CHECK_NEAR(inverse_hyperbolic_cosine(1),0);
    CHECK_NEAR(inverse_hyperbolic_tangent(0),0);

    printf("\n=== Numerical ===\n");
    CHECK_NEAR(newton_raphson(test_f_quad,test_df_quad,3,0.0001,100),2);
    CHECK_NEAR(newton_raphson(test_f_cubic,test_df_cubic,5,0.0001,100),3);
    CHECK_NEAR(bisection(test_f_quad,1,3,0.0001,100),2);
    CHECK_NEAR(bisection(test_f_line,0,5,0.0001,100),3);
    CHECK_NEAR(derivative(test_f_parabola,3,0.001),6);
    CHECK_NEAR(second_derivative(test_f_parabola,3,0.1),2);
    CHECK_NEAR(integral_simpson(test_f_parabola,0,1,100),0.333333);

    printf("\n=== Polynomial ===\n");
    float coeffs_quad[] = {-4, 0, 1};
    float coeffs_cubic[] = {-27, 0, 0, 1};
    float deriv_quad[2];
    float deriv_cubic[3];
    CHECK_NEAR(polynomial_eval(coeffs_quad,2,2),0);
    CHECK_NEAR(polynomial_eval(coeffs_quad,2,3),5);
    CHECK_NEAR(polynomial_eval(coeffs_cubic,3,3),0);
    polynomial_derivative(coeffs_quad,2,deriv_quad);
    CHECK_NEAR(deriv_quad[0],0);
    CHECK_NEAR(deriv_quad[1],2);
    polynomial_derivative(coeffs_cubic,3,deriv_cubic);
    CHECK_NEAR(deriv_cubic[0],0);
    CHECK_NEAR(deriv_cubic[1],0);
    CHECK_NEAR(deriv_cubic[2],3);
    CHECK_NEAR(polynomial_newton(coeffs_quad,2,3,0.0001,100),2);
    CHECK_NEAR(polynomial_newton(coeffs_cubic,3,5,0.0001,100),3);

    printf("\n=== Complex ===\n");
    cplx a = {1, 2};
    cplx b = {3, 4};
    cplx r = cplx_add(a, b);
    CHECK_NEAR(r.real, 4);
    CHECK_NEAR(r.imag, 6);
    r = cplx_sub(a, b);
    CHECK_NEAR(r.real, -2);
    CHECK_NEAR(r.imag, -2);
    r = cplx_mul(a, b);
    CHECK_NEAR(r.real, -5);
    CHECK_NEAR(r.imag, 10);
    r = cplx_div((cplx){5, 10}, (cplx){1, 2});
    CHECK_NEAR(r.real, 5);
    CHECK_NEAR(r.imag, 0);
    CHECK_NEAR(cplx_abs((cplx){3, 4}), 5);
    CHECK_NEAR(cplx_arg((cplx){1, 1}), math_pi/4);
    CHECK_NEAR(cplx_arg((cplx){-1, 0}), math_pi);
    r = cplx_conjugate((cplx){1, 2});
    CHECK_NEAR(r.real, 1);
    CHECK_NEAR(r.imag, -2);
    r = cplx_exponential((cplx){0, math_pi});
    CHECK_NEAR(r.real, -1);
    CHECK_NEAR(r.imag, 0);
    r = cplx_logarithm((cplx){math_e, 0});
    CHECK_NEAR(r.real, 1);
    CHECK_NEAR(r.imag, 0);
    r = cplx_power((cplx){math_e, 0}, (cplx){2, 0});
    CHECK_NEAR(r.real, math_e * math_e);
    CHECK_NEAR(r.imag, 0);

    TEST_SUMMARY();
}
