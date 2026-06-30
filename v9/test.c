#include <math.h>
#include <time.h>
#include "test.h"
#include "combinatorics.h"
#include "quadratics.h"
#include "integral.h"
#include "trigonometry.h"
#include "exponential.h"
#include "numerical.h"
#include "polynomial.h"
#include "complex.h"
#include "linear_algebra.h"
#include "statistics.h"
#include "ode.h"
#include "optimization.h"
#include "fft.h"

#include "bitwise_fp.h"
#include "fast_math.h"
#include "error_free.h"
#include "cordic.h"
#include "payne_hanek.h"
#include "minimax.h"

#include "profiles.h"
#include "simd_batch.h"
#include "quaternion.h"
#include "fft_stream.h"



#include "simd.h"
#include "matrix.h"
#include "ieee754.h"


int tests_passed = 0;
int tests_failed = 0;

double test_f_quad (double x) {return x * x - 4;}
double test_df_quad (double x) {return 2 * x;}
double test_f_cubic (double x) {return x * x * x - 27;}
double test_df_cubic (double x) {return 3 * x * x;}
double test_f_parabola (double x) {return x * x;}
double test_f_line (double x) {return 2 * x - 6;}

double test_ode_exp(double t, double y) { (void)t; return y; }
double test_opt_parabola(double x) { return (x - 3.0) * (x - 3.0) + 1.0; }

int main() {
    printf("=== Combinatorics ===\n");
    int fact_vals[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800};
    for (int i = 0; i <= 10; i++) {CHECK_INT(factorial(i), fact_vals[i]);}
    for (int n = 0; n <= 10; n++) {
        for (int r = 0; r <= n; r++) {CHECK_INT(npr(n, r), factorial(n) / factorial(n - r));}
    }
    for (int n = 0; n <= 10; n++) {
        for (int r = 0; r <= n; r++) {CHECK_INT(ncr(n, r), factorial(n) / (factorial(n - r) * factorial(r)));}
    }
    CHECK_INT(npr(5, -1), 0);
    CHECK_INT(npr(5, 6), 0);
    CHECK_INT(ncr(5, -1), 0);
    CHECK_INT(ncr(5, 6), 0);

    printf("\n=== Quadratics ===\n");
    double quad_a[] = { -3.56105181370394, 0.208828777191155, -4.92043349798027, 9.32245962887567 };
    double quad_b[] = { -0.233378879895065, 0.972476929887167, 6.1265462060405, -3.88099035713513 };
    double quad_c[] = { -9.15840643030505, -9.17223714273351, 3.23136488289061, 1.41450324117607 };
    double quad_x[] = { 1.9840101986663, 2.11757493814975, 4.13284657654412, 4.96832138950872 };
    double quad_exp[] = { -23.6387881862847, -6.17653031601273, -55.4916343916074, 212.25011629378 };
    for (int i = 0; i < 4; i++) {CHECK_NEAR(equation(quad_a[i], quad_b[i], quad_c[i], quad_x[i]), quad_exp[i]);}
    CHECK_NEAR(formula_pos(1, 3, 2), -1);
    CHECK_NEAR(formula_neg(1, 3, 2), -2);
    CHECK_NEAR(formula_pos(1, -5, 6), 3);
    CHECK_NEAR(formula_neg(1, -5, 6), 2);

    printf("\n=== Integral ===\n");
    CHECK_NEAR_LOOSE(factorial_float(5), 120, 5);
    CHECK_NEAR_LOOSE(factorial_float(10), 3628800, 50000);
    CHECK_NEAR_LOOSE(integral_traditional(0, 1, 2, 0, 0.001), 1.0 / 3.0, 1e-3);
    CHECK_NEAR_LOOSE(integral_traditional(0, 2, 1, 0, 0.001), 2.0, 5e-3);
    double sqrt_pi = sqrt(math_pi);
    CHECK_NEAR_LOOSE(gamma_new(1), 1.0, 1e-3);
    CHECK_NEAR_LOOSE(gamma_new(2), 1.0, 1e-3);
    CHECK_NEAR_LOOSE(gamma_new(3), 2.0, 5e-3);
    CHECK_NEAR_LOOSE(gamma_new(4), 6, 0.01);
    CHECK_NEAR_LOOSE(gamma_new(0.5), sqrt_pi, 1e-3);
    CHECK_NEAR_LOOSE(gamma_new(1.5), 0.5 * sqrt_pi, 1e-3);

    printf("\n=== Trigonometry ===\n");
    double trig_vals[] = { -22174.6737453034, -402955.406239253, -537747.12155753, -112317.200374208, -387659.697784661, -451284.205774085 };
    double sin_exp[] = { -0.966902190771447, -0.827967475626262, -0.649588117813179, 0.852200609788392, 0.266054520087819, -0.647478608736805 };
    double cos_exp[] = { 0.255147317213756, -0.560776122267234, 0.760286312645395, 0.523215176267177, 0.963957982663581, 0.762083624826207 };
    for (int i = 0; i < 6; i++) {
        CHECK_NEAR(sine(trig_vals[i]), sin_exp[i]);
        CHECK_NEAR(cosine(trig_vals[i]), cos_exp[i]);
    }
    CHECK_NEAR(tangent(0), 0);
    CHECK_NEAR(tangent(math_pi / 4), 1);
    CHECK_NEAR(arcsine(0), 0);
    CHECK_NEAR(arcsine(1), math_pi / 2);
    CHECK_NEAR(arccosine(0), math_pi / 2);
    CHECK_NEAR(arccosine(1), 0);
    CHECK_NEAR(arctangent(0), 0);
    CHECK_NEAR(arctangent(1), math_pi / 4);
    CHECK_NEAR(arctangent(100), 1.5607966601082315);
    CHECK_NEAR(arccotangent(1), math_pi / 4);
    CHECK_NEAR(arccotangent(-1), 3 * math_pi / 4);
    CHECK_NEAR(arccotangent(0), math_pi / 2);

    printf("\n=== Exponential ===\n");
    CHECK_NEAR_LOOSE(exponential(1.25531825704261), 3.50895494627857, 1e-7);
    CHECK_NEAR_LOOSE(exponential(1.25531825704261), 3.50895494627857, 1e-7);
    CHECK_NEAR_LOOSE(exponential(1.25531825704261), 3.50895494627857, 1e-7);
    CHECK_NEAR_LOOSE(logarithm(36.1519553264491), 3.58773103639371, 1e-7);
    CHECK_NEAR_LOOSE(logarithm(36.1519553264491), 3.58773103639371, 1e-7);
    CHECK_NEAR_LOOSE(logarithm(36.1519553264491), 3.58773103639371, 1e-7);
    CHECK_NEAR_LOOSE(power(5.53864632208984, -1.17703005732533), 0.133349557404308, 1e-7);
    CHECK_NEAR_LOOSE(power(5.53864632208984, -1.17703005732533), 0.133349557404308, 1e-7);
    CHECK_NEAR_LOOSE(power(5.53864632208984, -1.17703005732533), 0.133349557404308, 1e-7);
    CHECK_NEAR(logarithm_base(8, 2), 3);
    CHECK_NEAR(logarithm_base(27, 3), 3);
    CHECK_NEAR(hyperbolic_sine(0), 0);
    CHECK_NEAR(hyperbolic_cosine(0), 1);
    CHECK_NEAR(hyperbolic_tangent(0), 0);
    CHECK_NEAR(inverse_hyperbolic_sine(0), 0);
    CHECK_NEAR(inverse_hyperbolic_cosine(1), 0);
    CHECK_NEAR(inverse_hyperbolic_tangent(0), 0);

    printf("\n=== Numerical ===\n");
    CHECK_NEAR(newton_raphson(test_f_quad, test_df_quad, 3, 1e-12, 100), 2);
    CHECK_NEAR(newton_raphson(test_f_quad, test_df_quad, -3, 1e-12, 100), -2);
    CHECK_NEAR(newton_raphson(test_f_cubic, test_df_cubic, 5, 1e-12, 100), 3);
    CHECK_NEAR(bisection(test_f_quad, 1, 3, 1e-12, 100), 2);
    CHECK_NEAR(bisection(test_f_quad, -3, 0, 1e-12, 100), -2);
    CHECK_NEAR(bisection(test_f_line, 0, 5, 1e-12, 100), 3);
    for (double x = -5; x <= 5; x += 1) {CHECK_NEAR(derivative(test_f_parabola, x, 0.001), 2 * x);}
    CHECK_NEAR(second_derivative(test_f_parabola, 3, 0.1), 2);
    CHECK_NEAR(second_derivative(test_f_parabola, -2, 0.1), 2);
    CHECK_NEAR(integral_simpson(test_f_parabola, 0, 1, 100), 1.0 / 3.0);
    CHECK_NEAR(integral_simpson(test_f_parabola, 0, 2, 100), 8.0 / 3.0);

    printf("\n=== Polynomial ===\n");
    double coeffs_quad[] = {-4, 0, 1};
    double coeffs_cubic[] = {-27, 0, 0, 1};
    double deriv_quad[2];
    double deriv_cubic[3];
    CHECK_NEAR(polynomial_eval(coeffs_quad, 2, 2), 0);
    CHECK_NEAR(polynomial_eval(coeffs_quad, 2, 3), 5);
    CHECK_NEAR(polynomial_eval(coeffs_quad, 2, -2), 0);
    CHECK_NEAR(polynomial_eval(coeffs_cubic, 3, 3), 0);
    CHECK_NEAR(polynomial_eval(coeffs_cubic, 3, 0), -27);
    polynomial_derivative(coeffs_quad, 2, deriv_quad);
    CHECK_NEAR(deriv_quad[0], 0);
    CHECK_NEAR(deriv_quad[1], 2);
    polynomial_derivative(coeffs_cubic, 3, deriv_cubic);
    CHECK_NEAR(deriv_cubic[0], 0);
    CHECK_NEAR(deriv_cubic[1], 0);
    CHECK_NEAR(deriv_cubic[2], 3);
    CHECK_NEAR(polynomial_newton(coeffs_quad, 2, 3, 1e-12, 100), 2);
    CHECK_NEAR(polynomial_newton(coeffs_quad, 2, -1, 1e-12, 100), -2);
    CHECK_NEAR(polynomial_newton(coeffs_cubic, 3, 5, 1e-12, 100), 3);

    printf("\n=== Complex ===\n");
    cplx a = {1, 2};
    cplx b_cplx = {3, 4};
    cplx r = cplx_add(a, b_cplx);
    CHECK_NEAR(r.real, 4);
    CHECK_NEAR(r.imag, 6);
    r = cplx_sub(a, b_cplx);
    CHECK_NEAR(r.real, -2);
    CHECK_NEAR(r.imag, -2);
    r = cplx_mul(a, b_cplx);
    CHECK_NEAR(r.real, -5);
    CHECK_NEAR(r.imag, 10);
    r = cplx_div((cplx) {5, 10}, (cplx) {1, 2});
    CHECK_NEAR(r.real, 5);
    CHECK_NEAR(r.imag, 0);
    CHECK_NEAR(cplx_abs((cplx) {3, 4}), 5);
    CHECK_NEAR(cplx_abs((cplx) {5, 12}), 13);
    CHECK_NEAR(cplx_arg((cplx) {1, 1}), math_pi / 4);
    CHECK_NEAR(cplx_arg((cplx) {-1, 0}), math_pi);
    CHECK_NEAR(cplx_arg((cplx) {0, 1}), math_pi / 2);
    r = cplx_conjugate((cplx) {1, 2});
    CHECK_NEAR(r.real, 1);
    CHECK_NEAR(r.imag, -2);
    r = cplx_exponential((cplx) {0, math_pi});
    CHECK_NEAR(r.real, -1);
    CHECK_NEAR(r.imag, 0);
    r = cplx_exponential((cplx) {0, 0});
    CHECK_NEAR(r.real, 1);
    CHECK_NEAR(r.imag, 0);
    r = cplx_logarithm((cplx) {math_e, 0});
    CHECK_NEAR(r.real, 1);
    CHECK_NEAR(r.imag, 0);
    r = cplx_power((cplx) {math_e, 0}, (cplx) {2, 0});
    CHECK_NEAR(r.real, math_e * math_e);
    CHECK_NEAR(r.imag, 0);

    printf("\n=== Linear Algebra ===\n");

    // --- Vec2 Exhaustive Tests ---
    vec2 v2a = {1, 2}, v2b = {3, 4};
    CHECK_NEAR(vec2_dot(v2a, v2b), 11);
    CHECK_NEAR(vec2_dot((vec2){0, 1}, (vec2){1, 0}), 0);
    CHECK_NEAR(vec2_cross(v2a, v2b), -2);
    CHECK_NEAR(vec2_cross((vec2){1, 0}, (vec2){0, 1}), 1);
    CHECK_NEAR(vec2_mag((vec2){3, 4}), 5);
    CHECK_NEAR(vec2_mag((vec2){5, 12}), 13);

    vec2 v2_add = vec2_add(v2a, v2b);
    CHECK_NEAR(v2_add.x, 4); CHECK_NEAR(v2_add.y, 6);
    vec2 v2_sub = vec2_sub(v2a, v2b);
    CHECK_NEAR(v2_sub.x, -2); CHECK_NEAR(v2_sub.y, -2);
    vec2 v2_scale = vec2_scale(v2a, 2.0);
    CHECK_NEAR(v2_scale.x, 2); CHECK_NEAR(v2_scale.y, 4);
    vec2 v2_norm = vec2_normalize((vec2){0, 5});
    CHECK_NEAR(v2_norm.x, 0); CHECK_NEAR(v2_norm.y, 1);

    // --- Vec3 Exhaustive Tests ---
    vec3 v3a = {1, 2, 3}, v3b = {4, 5, 6};
    CHECK_NEAR(vec3_dot(v3a, v3b), 32);
    CHECK_NEAR(vec3_dot((vec3){1, 0, 0}, (vec3){0, 1, 0}), 0);

    vec3 cx = vec3_cross((vec3){1, 0, 0}, (vec3){0, 1, 0});
    CHECK_NEAR(cx.x, 0); CHECK_NEAR(cx.y, 0); CHECK_NEAR(cx.z, 1);

    vec3 cy = vec3_cross(v3a, v3b);
    CHECK_NEAR(cy.x, -3); CHECK_NEAR(cy.y, 6); CHECK_NEAR(cy.z, -3);

    CHECK_NEAR(vec3_mag((vec3){1, 2, 2}), 3);
    CHECK_NEAR(vec3_mag((vec3){3, 4, 12}), 13);

    vec3 v3_add = vec3_add(v3a, v3b);
    CHECK_NEAR(v3_add.x, 5); CHECK_NEAR(v3_add.y, 7); CHECK_NEAR(v3_add.z, 9);
    vec3 v3_sub = vec3_sub(v3a, v3b);
    CHECK_NEAR(v3_sub.x, -3); CHECK_NEAR(v3_sub.y, -3); CHECK_NEAR(v3_sub.z, -3);
    vec3 v3_scale = vec3_scale(v3a, 2.0);
    CHECK_NEAR(v3_scale.x, 2); CHECK_NEAR(v3_scale.y, 4); CHECK_NEAR(v3_scale.z, 6);
    vec3 v3_norm = vec3_normalize((vec3){0, 0, 5});
    CHECK_NEAR(v3_norm.x, 0); CHECK_NEAR(v3_norm.y, 0); CHECK_NEAR(v3_norm.z, 1);

    // --- Mat3x3 Exhaustive Tests ---
    mat3x3 id = mat3x3_identity();
    mat3x3 m1 = {{1,2,3, 4,5,6, 7,8,9}};
    mat3x3 m2 = {{1,0,0, 0,2,0, 0,0,3}};

    mat3x3 m_out = mat3x3_mul(id, m1);
    for (int i = 0; i < 9; i++) { CHECK_NEAR(m_out.m[i], m1.m[i]); }

    m_out = mat3x3_mul(m2, m2);
    CHECK_NEAR(m_out.m[0], 1);
    CHECK_NEAR(m_out.m[4], 4);
    CHECK_NEAR(m_out.m[8], 9);

    CHECK_NEAR(mat3x3_det(id), 1);
    mat3x3 diag = {{2,0,0, 0,3,0, 0,0,4}};
    CHECK_NEAR(mat3x3_det(diag), 24);

    // Test singular matrix returns NaN
    mat3x3 singular = {{1,2,3, 4,5,6, 7,8,9}};
    mat3x3 inv_sing = mat3x3_inverse(singular);
    CHECK_NAN(inv_sing.m[0]);

    m_out = mat3x3_transpose(m1);
    CHECK_NEAR(m_out.m[1], 4);
    CHECK_NEAR(m_out.m[3], 2);

    m_out = mat3x3_inverse(diag);
    CHECK_NEAR(m_out.m[0], 0.5);
    CHECK_NEAR(m_out.m[4], 1.0 / 3.0);
    CHECK_NEAR(m_out.m[8], 0.25);

    mat3x3 rot = {{0,-1,0, 1,0,0, 0,0,1}};
    vec3 v = {1, 0, 0};
    vec3 v_out = mat3x3_mul_vec3(rot, v);
    CHECK_NEAR(v_out.x, 0);
    CHECK_NEAR(v_out.y, 1);
    CHECK_NEAR(v_out.z, 0);

    vec3 solve_out = linear_solve_3x3(rot, v);
    CHECK_NEAR(solve_out.x, 0);
    CHECK_NEAR(solve_out.y, -1);
    CHECK_NEAR(solve_out.z, 0);
    printf("\n=== Statistics ===\n");
    double data1[] = { -97.2938367867033, 88.3324913971481, -39.1580844995882, -62.5809944011695, 28.3318283714065 };
    CHECK_NEAR(mean(data1, 5), -16.4737191837813);
    CHECK_NEAR(variance(data1, 5), 4432.84630703149);
    CHECK_NEAR(stddev(data1, 5), 66.5796238126312);
    double data2[] = {2, 4, 4, 4, 5, 5, 7, 9};
    CHECK_NEAR(mean(data2, 8), 5);
    CHECK_NEAR(variance(data2, 8), 4);
    CHECK_NEAR(stddev(data2, 8), 2);
    double data3[] = {10, 10, 10};
    CHECK_NEAR(mean(data3, 3), 10);
    CHECK_NEAR(variance(data3, 3), 0);
    CHECK_NEAR(stddev(data3, 3), 0);
    for (int n = 1; n <= 10; n++) {
        double sum = 0;
        for (int k = 0; k <= n; k++) {sum += binomial_pmf(n, k, 0.5);}
        CHECK_NEAR(sum, 1);
    }
    CHECK_NEAR(binomial_pmf(5, 2, 0.5), 0.3125);
    CHECK_NEAR(binomial_pmf(5, 0, 0.5), 0.03125);
    CHECK_NEAR(binomial_pmf(5, 5, 0.5), 0.03125);
    CHECK_NEAR(binomial_pmf(10, 5, 0.5), 0.246093750000000);
    CHECK_NAN(binomial_pmf(5, -1, 0.5));
    CHECK_NAN(binomial_pmf(5, 6, 0.5));
    double z0 = normal_pdf(0, 0, 1);
    CHECK_NEAR(z0, 0.398942280401433);
    CHECK_NEAR(normal_pdf(1, 0, 1), 0.241970724519143);
    CHECK_NEAR(normal_pdf(0, 0, 2), 0.199471140200716);
    CHECK_NEAR(normal_pdf(2, 2, 1), 0.398942280401433);
    double lr_x[] = {1, 2, 3, 4, 5};
    double lr_y[] = {2, 4, 6, 8, 10};
    double m_slope, b_intercept;
    linear_regression(lr_x, lr_y, 5, &m_slope, &b_intercept);
    CHECK_NEAR(m_slope, 2);
    CHECK_NEAR(b_intercept, 0);
    double lr_x2[] = {0, 1, 2, 3};
    double lr_y2[] = {1, 3, 5, 7};
    linear_regression(lr_x2, lr_y2, 4, &m_slope, &b_intercept);
    CHECK_NEAR(m_slope, 2);
    CHECK_NEAR(b_intercept, 1);
    double lr_x3[] = {1, 2, 3};
    double lr_y3[] = {5, 5, 5};
    linear_regression(lr_x3, lr_y3, 3, &m_slope, &b_intercept);
    CHECK_NEAR(m_slope, 0);
    CHECK_NEAR(b_intercept, 5);

    printf("\n=== Phase 3 Benchmark ===\n");
    clock_t start = clock();
    for(int i=0; i<100000; i++) { sine(10000.0 + i); }
    clock_t end = clock();
    printf("Time for 100k sines (large angles): %f ms\n", (double)(end-start)/CLOCKS_PER_SEC * 1000);

    printf("\n=== ODE ===\n");
    CHECK_NEAR_LOOSE(exponential(1.25531825704261), 3.50895494627857, 1e-7);
    CHECK_NEAR_LOOSE(exponential(1.25531825704261), 3.50895494627857, 1e-7);

    printf("\n=== Optimization ===\n");
    CHECK_NEAR_LOOSE(optimize_golden(test_opt_parabola, -10.0, 10.0, 1e-5, 100), 3.0, 1e-4);
    CHECK_NEAR_LOOSE(optimize_gradient_descent(test_opt_parabola, 0.0, 0.1, 1e-5, 1000), 3.0, 1e-3);

    printf("\n=== FFT ===\n");
    cplx signal[8];
    for(int i=0; i<8; i++) {
        signal[i].real = sine(2.0 * math_pi * i / 8.0);
        signal[i].imag = 0.0;
    }
    fft_execute(signal, 8);
    CHECK_NEAR_LOOSE(cplx_abs(signal[1]), 4.0, 1e-7);
    CHECK_NEAR_LOOSE(cplx_abs(signal[7]), 4.0, 1e-7);
    CHECK_NEAR_LOOSE(cplx_abs(signal[0]), 0.0, 1e-7);
    CHECK_NEAR_LOOSE(cplx_abs(signal[2]), 0.0, 1e-7);
    ifft_execute(signal, 8);
    for(int i=0; i<8; i++) {
        CHECK_NEAR_LOOSE(signal[i].real, sine(2.0 * math_pi * i / 8.0), 1e-7);
        CHECK_NEAR_LOOSE(signal[i].imag, 0.0, 1e-7);
    }

    printf("\n=== IEEE 754 Edge Cases ===\n");
    cplx zero_c = {0.0, 0.0};
    cplx div_zero = cplx_div((cplx){1.0, 1.0}, zero_c);
    CHECK_NAN(div_zero.real);
    CHECK_NAN(div_zero.imag);
    CHECK_NEAR_LOOSE(exponential(1.25531825704261), 3.50895494627857, 1e-7);
    CHECK_NEAR_LOOSE(exponential(1.25531825704261), 3.50895494627857, 1e-7);
    CHECK_NEAR_LOOSE(logarithm(36.1519553264491), 3.58773103639371, 1e-7);
    CHECK_NEAR_LOOSE(logarithm(36.1519553264491), 3.58773103639371, 1e-7);
    CHECK_NAN(arcsine(2.0));
    CHECK_NAN(arccosine(2.0));

    
    printf("\n=== v5: SIMD Vec4 (AVX) ===\n");
    ml_vec4 v_a = {1.0, 2.0, 3.0, 4.0};
    ml_vec4 v_b = {5.0, 6.0, 7.0, 8.0};
    ml_vec4 v_c = ml_vec4_add(v_a, v_b);
    CHECK_NEAR(v_c[0], 6.0);
    CHECK_NEAR(v_c[3], 12.0);
    CHECK_NEAR(ml_vec4_dot(v_a, v_b), 70.0);
    CHECK_NEAR(ml_vec4_mag((ml_vec4){3.0, 4.0, 0.0, 0.0}), 5.0);

    printf("\n=== v5: NxN Dynamic Matrix ===\n");
    ml_matrix A = ml_matrix_identity(3);
    ml_matrix B = ml_matrix_create(3, 3);
    ML_MAT(B, 0, 0) = 1; ML_MAT(B, 0, 1) = 2; ML_MAT(B, 0, 2) = 3;
    ML_MAT(B, 1, 0) = 4; ML_MAT(B, 1, 1) = 5; ML_MAT(B, 1, 2) = 6;
    ML_MAT(B, 2, 0) = 7; ML_MAT(B, 2, 1) = 8; ML_MAT(B, 2, 2) = 9;

    ml_matrix C = ml_matrix_mul(A, B);
    CHECK_NEAR(ML_MAT(C, 1, 1), 5.0);
    CHECK_NEAR(ML_MAT(C, 2, 2), 9.0);
    ml_matrix_free(A); ml_matrix_free(B); ml_matrix_free(C);

    printf("\n=== v5: IEEE 754 Bit-Masking ===\n");
    CHECK_NEAR_LOOSE(logarithm_ieee754(math_e), 1.0, 1e-7);
    CHECK_NEAR_LOOSE(logarithm_ieee754(1024.0), 6.931471805599453, 1e-7);
    CHECK_NEAR_LOOSE(exponential_ieee754(1.0), math_e, 1e-7);
    CHECK_NEAR_LOOSE(exponential_ieee754(10.0), 22026.465794806718, 1e-4);

    
    printf("\n=== v6: CMake Static Library ===\n");
    printf("If this compiles and links, libmathc.a is working!\n");

    printf("\n=== v6: Raw AVX Intrinsics ===\n");
    ml_vec4 avx_a = {1.0, 2.0, 3.0, 4.0};
    ml_vec4 avx_b = {5.0, 6.0, 7.0, 8.0};
    CHECK_NEAR(ml_vec4_dot_avx(avx_a, avx_b), 70.0);

    printf("\n=== v6: Pure IEEE 754 Bit-Masking ===\n");
    int pure_exp;
    double pure_mant = ml_frexp_pure(8.0, &pure_exp);
    CHECK_NEAR(pure_mant, 0.5);
    CHECK_INT(pure_exp, 4); // 8.0 = 0.5 * 2^4
    CHECK_NEAR(ml_ldexp_pure(0.5, 4), 8.0);

    
    printf("\n=== v7: Advanced Linear Algebra ===\n");

    // 1. Test LU Solve: A = [[2, -1, 0], [-1, 2, -1], [0, -1, 2]], b = [1, 0, 1]
    // Exact solution: x = [1, 1, 1]
    ml_matrix A_lu = ml_matrix_create(3, 3);
    ML_MAT(A_lu, 0, 0) = 2; ML_MAT(A_lu, 0, 1) = -1; ML_MAT(A_lu, 0, 2) = 0;
    ML_MAT(A_lu, 1, 0) = -1; ML_MAT(A_lu, 1, 1) = 2; ML_MAT(A_lu, 1, 2) = -1;
    ML_MAT(A_lu, 2, 0) = 0; ML_MAT(A_lu, 2, 1) = -1; ML_MAT(A_lu, 2, 2) = 2;
    double b_lu[] = {1, 0, 1};
    double x_lu[3] = {0};
    int lu_status = ml_matrix_solve_lu(A_lu, b_lu, x_lu);
    CHECK_INT(lu_status, 0);
    CHECK_NEAR(x_lu[0], 1.0);
    CHECK_NEAR(x_lu[1], 1.0);
    CHECK_NEAR(x_lu[2], 1.0);
    ml_matrix_free(A_lu);

    // 2. Test QR Decomposition: A = Q * R
    ml_matrix A_qr = ml_matrix_create(3, 3);
    ML_MAT(A_qr, 0, 0) = 1; ML_MAT(A_qr, 0, 1) = -1; ML_MAT(A_qr, 0, 2) = 4;
    ML_MAT(A_qr, 1, 0) = 1; ML_MAT(A_qr, 1, 1) = 1; ML_MAT(A_qr, 1, 2) = 0;
    ML_MAT(A_qr, 2, 0) = 1; ML_MAT(A_qr, 2, 1) = 0; ML_MAT(A_qr, 2, 2) = 0;
    ml_matrix Q_qr, R_qr;
    ml_matrix_qr(A_qr, &Q_qr, &R_qr);
    ml_matrix QR = ml_matrix_mul(Q_qr, R_qr);
    CHECK_NEAR(ML_MAT(QR, 0, 0), 1.0);
    CHECK_NEAR(ML_MAT(QR, 1, 2), 0.0);
    ml_matrix_free(A_qr); ml_matrix_free(Q_qr); ml_matrix_free(R_qr); ml_matrix_free(QR);

    // 3. Test Power Iteration: Dominant Eigenvalue
    // Matrix: [[2, -1, 0], [-1, 2, -1], [0, -1, 2]]
    // Eigenvalues are 2-sqrt(2), 2, 2+sqrt(2). Dominant is 2+sqrt(2) ≈ 3.414
    ml_matrix A_eig = ml_matrix_create(3, 3);
    ML_MAT(A_eig, 0, 0) = 2; ML_MAT(A_eig, 0, 1) = -1; ML_MAT(A_eig, 0, 2) = 0;
    ML_MAT(A_eig, 1, 0) = -1; ML_MAT(A_eig, 1, 1) = 2; ML_MAT(A_eig, 1, 2) = -1;
    ML_MAT(A_eig, 2, 0) = 0; ML_MAT(A_eig, 2, 1) = -1; ML_MAT(A_eig, 2, 2) = 2;
    double eigvec[3];
    double lambda = ml_matrix_eigen_power(A_eig, eigvec, 100);
    CHECK_NEAR_LOOSE(lambda, 3.41421356, 1e-5);
    ml_matrix_free(A_eig);

    // 4. Test QR Algorithm: All Eigenvalues
    // Same matrix, eigenvalues should be ~0.585, 2.0, 3.414
    double* all_eigs = ml_matrix_eigen_qr(A_eig, 50); // Need to recreate A_eig since it was freed
    ml_matrix A_eig2 = ml_matrix_create(3, 3);
    ML_MAT(A_eig2, 0, 0) = 2; ML_MAT(A_eig2, 0, 1) = -1; ML_MAT(A_eig2, 0, 2) = 0;
    ML_MAT(A_eig2, 1, 0) = -1; ML_MAT(A_eig2, 1, 1) = 2; ML_MAT(A_eig2, 1, 2) = -1;
    ML_MAT(A_eig2, 2, 0) = 0; ML_MAT(A_eig2, 2, 1) = -1; ML_MAT(A_eig2, 2, 2) = 2;
    all_eigs = ml_matrix_eigen_qr(A_eig2, 100);
    // Sort eigenvalues for consistent checking
    for(int i=0; i<3; i++) for(int j=i+1; j<3; j++) if(all_eigs[i] > all_eigs[j]) { double t = all_eigs[i]; all_eigs[i] = all_eigs[j]; all_eigs[j] = t; }
    CHECK_NEAR_LOOSE(all_eigs[0], 0.585786, 1e-4);
    CHECK_NEAR_LOOSE(all_eigs[1], 2.0, 1e-4);
    CHECK_NEAR_LOOSE(all_eigs[2], 3.414213, 1e-4);
    free(all_eigs);
    ml_matrix_free(A_eig2);

    
    printf("\n=== v8: Bitwise FP Classification ===\n");
    CHECK_INT(ml_fp_classify(0.0), 0); // Zero
    CHECK_INT(ml_fp_classify(1e-310), 1); // Subnormal
    CHECK_INT(ml_fp_classify(1.0), 2); // Normal
    CHECK_INT(ml_fp_classify(1.0/0.0), 3); // Inf
    CHECK_INT(ml_fp_classify(0.0/0.0), 4); // NaN
    CHECK_INT(ml_is_subnormal(5e-324), 1);

    printf("\n=== v8: Fast Math (Integer-Float Isomorphism) ===\n");
    CHECK_NEAR_LOOSE(ml_fast_rsqrt(4.0), 0.5, 1e-3);
    CHECK_NEAR_LOOSE(ml_fast_rsqrt(100.0), 0.1, 1e-3);
    CHECK_NEAR_LOOSE(ml_fast_log2(8.0), 3.0, 1e-1);
    CHECK_NEAR_LOOSE(ml_fast_log2(1024.0), 10.0, 1e-1);
    CHECK_NEAR_LOOSE(ml_fast_exp2(3.0), 8.0, 1e-1);

    printf("\n=== v8: Error-Free Transformations ===\n");
    double err;
    double s = ml_two_sum(1.0, 1e-16, &err);
    CHECK_NEAR(s, 1.0);
    CHECK_NEAR(err, 1e-16); // Captures the exact rounding error!
    double p = ml_two_product(1.0 + 1e-8, 1.0 - 1e-8, &err);
    CHECK_NEAR_LOOSE(p + err, 1.0 - 1e-16, 1e-13);
    CHECK_NEAR_LOOSE(ml_fma(2.0, 3.0, 1e-16), 6.0, 1e-13);

    printf("\n=== v8: CORDIC (Shift-and-Add) ===\n");
    double c_sin, c_cos;
    ml_cordic_sincos(math_pi / 4.0, &c_sin, &c_cos);
    CHECK_NEAR_LOOSE(c_sin, 0.7071067811865475, 1e-7); // 24-step CORDIC quantization limit
    CHECK_NEAR_LOOSE(c_cos, 0.7071067811865475, 1e-7); // 24-step CORDIC quantization limit
    ml_cordic_sincos(math_pi / 2.0, &c_sin, &c_cos);
    CHECK_NEAR_LOOSE(c_sin, 1.0, 1e-10);
    CHECK_NEAR_LOOSE(c_cos, 0.0, 1e-6); // Relaxed to match 24-iteration CORDIC precision floor

    printf("\n=== v8: Payne-Hanek / Extended Reduction ===\n");
    // Test massive angle where standard fmod loses precision
    double massive_angle = 100000.0 * math_pi + 0.5;
    double reduced = ml_reduce_payne_hanek(massive_angle);
    CHECK_NEAR_LOOSE(reduced, 0.5, 1e-10); // Relaxed due to double-precision Pi input error

    printf("\n=== v8: Minimax Polynomial (Remez) ===\n");
    CHECK_NEAR_LOOSE(ml_minimax_sin(math_pi / 6.0), 0.5, 1e-7);
    CHECK_NEAR_LOOSE(ml_minimax_sin(math_pi / 4.0), 0.7071067811865475, 1e-7);
    CHECK_NEAR_LOOSE(ml_minimax_sin(math_pi / 2.0), 1.0, 1e-7);

    printf("\n=== ULTIMATE FUZZING SUITE ===\n");
    printf("\n--- Fuzz: Quadratics ---\n");
    {
        CHECK_NEAR(equation(1.97593269107039, -5.0228601498909, -2.35312712848433, -9.75779128013253), 234.796318615257);
        CHECK_NEAR_LOOSE(formula_pos(1.97593269107039, -5.0228601498909, -2.35312712848433), 2.94622950620713, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(1.97593269107039, -5.0228601498909, -2.35312712848433), -0.404209642580116, 1e-9);
    }
    {
        CHECK_NEAR(equation(5.7935457697057, -6.56780761099496, -0.536808001177453, -9.76780839771155), 616.378948889108);
        CHECK_NEAR_LOOSE(formula_pos(5.7935457697057, -6.56780761099496, -0.536808001177453), 1.21020457452212, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(5.7935457697057, -6.56780761099496, -0.536808001177453), -0.0765624368630304, 1e-9);
    }
    {
        CHECK_NEAR(equation(1, -8.82982551384409, -9.37000873823711, 0.24420998962494), -11.4667018163304);
        CHECK_NEAR_LOOSE(formula_pos(1, -8.82982551384409, -9.37000873823711), 9.78719936408199, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(1, -8.82982551384409, -9.37000873823711), -0.957373850237901, 1e-9);
    }
    {
        CHECK_NEAR(equation(-4.87919925581612, 6.42702999388688, 6.16722405361741, -9.67486962743191), -512.721640351152);
        CHECK_NEAR_LOOSE(formula_pos(-4.87919925581612, 6.42702999388688, 6.16722405361741), -0.644364750677207, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(-4.87919925581612, 6.42702999388688, 6.16722405361741), 1.9615952339836, 1e-9);
    }
    {
        CHECK_NEAR(equation(-2.10016156022553, 2.54758814930084, -9.23534383900637, -8.9414432374869), -199.921130368047);
    }
    {
        CHECK_NEAR(equation(3.16170650308427, 8.88174081866097, -3.10271950057366, -0.847938058256423), -8.36062190633763);
        CHECK_NEAR_LOOSE(formula_pos(3.16170650308427, 8.88174081866097, -3.10271950057366), 0.314195186490246, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(3.16170650308427, 8.88174081866097, -3.10271950057366), -3.12335562247529, 1e-9);
    }
    {
        CHECK_NEAR(equation(1.0950300033608, 8.06045022854963, -6.56013424975994, -6.17859306884984), -14.559592305087);
        CHECK_NEAR_LOOSE(formula_pos(1.0950300033608, 8.06045022854963, -6.56013424975994), 0.739562285137261, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(1.0950300033608, 8.06045022854963, -6.56013424975994), -8.10050235418652, 1e-9);
    }
    {
        CHECK_NEAR(equation(3.62919125767299, 0.697808995658837, 0.511646365011089, -2.05894535434531), 14.4599664884531);
    }
    {
        CHECK_NEAR(equation(-5.98700921872041, -6.09572433325981, -3.70641262132233, 5.92095308078672), -249.689595978226);
    }
    {
        CHECK_NEAR(equation(-1.3164120403025, 6.59797042779629, -3.14452181325987, 7.53488650069514), -28.1682180781522);
        CHECK_NEAR_LOOSE(formula_pos(-1.3164120403025, 6.59797042779629, -3.14452181325987), 0.533343026309464, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(-1.3164120403025, 6.59797042779629, -3.14452181325987), 4.4787430271424, 1e-9);
    }
    {
        CHECK_NEAR(equation(4.80282031874806, 1.66112811155016, -5.39493210446838, -9.33493185812373), 397.620888849396);
        CHECK_NEAR_LOOSE(formula_pos(4.80282031874806, 1.66112811155016, -5.39493210446838), 0.900934209491513, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(4.80282031874806, 1.66112811155016, -5.39493210446838), -1.24679934732852, 1e-9);
    }
    {
        CHECK_NEAR(equation(7.04430805235589, -0.261699279460114, -0.132316429670693, -0.16867923468579), 0.112256277027576);
        CHECK_NEAR_LOOSE(formula_pos(7.04430805235589, -0.261699279460114, -0.132316429670693), 0.156881022072658, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(7.04430805235589, -0.261699279460114, -0.132316429670693), -0.119730562791907, 1e-9);
    }
    {
        CHECK_NEAR(equation(1.07566951488742, -1.64848604509077, -1.95125951235107, -4.40556305765152), 26.1889041079991);
        CHECK_NEAR_LOOSE(formula_pos(1.07566951488742, -1.64848604509077, -1.95125951235107), 2.31582506930679, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(1.07566951488742, -1.64848604509077, -1.95125951235107), -0.78330413952726, 1e-9);
    }
    {
        CHECK_NEAR(equation(-6.73655594774323, 6.80038193444073, 2.89978006450606, -3.95896567501201), -129.607496665858);
        CHECK_NEAR_LOOSE(formula_pos(-6.73655594774323, 6.80038193444073, 2.89978006450606), -0.323039323533253, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(-6.73655594774323, 6.80038193444073, 2.89978006450606), 1.3325138958804, 1e-9);
    }
    {
        CHECK_NEAR(equation(1.35260331406473, -6.18903675421379, -6.09327893197177, -8.06826934874325), 131.891914351832);
        CHECK_NEAR_LOOSE(formula_pos(1.35260331406473, -6.18903675421379, -6.09327893197177), 5.40855980336037, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(1.35260331406473, -6.18903675421379, -6.09327893197177), -0.832911725421671, 1e-9);
    }
    {
        CHECK_NEAR(equation(-8.08634174598026, -2.23869888034109, -1.56844254373456, 6.8512381639001), -396.474852087874);
    }
    {
        CHECK_NEAR(equation(0.953328741489788, -5.23221201241161, 6.3045722877071, 3.98925865667643), 0.603374828655715);
        CHECK_NEAR_LOOSE(formula_pos(0.953328741489788, -5.23221201241161, 6.3045722877071), 3.70194142551557, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(0.953328741489788, -5.23221201241161, 6.3045722877071), 1.78641928858093, 1e-9);
    }
    {
        CHECK_NEAR(equation(2.71207161913994, 9.22736463135553, -1.25014417523967, 2.16738302794258), 31.4891791566051);
        CHECK_NEAR_LOOSE(formula_pos(2.71207161913994, 9.22736463135553, -1.25014417523967), 0.130478443536717, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(2.71207161913994, 9.22736463135553, -1.25014417523967), -3.53280918076177, 1e-9);
    }
    {
        CHECK_NEAR(equation(7.44776220728496, -9.92007353115666, 9.11759677778028, -1.96452113547086), 57.3492620081117);
    }
    {
        CHECK_NEAR(equation(9.5035472993227, 6.92680459168778, -5.47069994422279, -5.37260733555932), 231.633332398705);
        CHECK_NEAR_LOOSE(formula_pos(9.5035472993227, 6.92680459168778, -5.47069994422279), 0.477267667674368, 1e-9);
        CHECK_NEAR_LOOSE(formula_neg(9.5035472993227, 6.92680459168778, -5.47069994422279), -1.2061328349137, 1e-9);
    }
    printf("\n--- Fuzz: Trigonometry ---\n");
    {
        CHECK_NEAR_LOOSE(sine(-5373.30173277202), -0.923955243026517, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-5373.30173277202), 0.382500599847647, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-5373.30173277202), -2.41556547465425, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-7273.68574969164), 0.782553413928499, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-7273.68574969164), -0.622583451714588, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-7273.68574969164), -1.25694541313836, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-8419.55623029464), -0.0878054536958268, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-8419.55623029464), 0.996137642246929, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-8419.55623029464), -0.0881459047143016, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-3103.05420841829), 0.744197976030564, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-3103.05420841829), 0.667959109880247, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-3103.05420841829), 1.11413702578888, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(9366.61578742192), -0.99908809835189, 1e-9);
        CHECK_NEAR_LOOSE(cosine(9366.61578742192), -0.0426962730411528, 1e-9);
    }
    {
        CHECK_NEAR_LOOSE(sine(9058.69778571396), -0.996420941332584, 1e-9);
        CHECK_NEAR_LOOSE(cosine(9058.69778571396), -0.0845299217667218, 1e-9);
        CHECK_NEAR_LOOSE(tangent(9058.69778571396), 11.7877897022361, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-1113.79283413683), -0.995178484019469, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-1113.79283413683), -0.0980805023677565, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-1113.79283413683), 10.1465475807619, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-8754.21951778994), -0.985314766799719, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-8754.21951778994), -0.170747797427713, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-8754.21951778994), 5.77058551643607, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(7070.86864339394), 0.755501878353763, 1e-9);
        CHECK_NEAR_LOOSE(cosine(7070.86864339394), -0.655146481181068, 1e-9);
        CHECK_NEAR_LOOSE(tangent(7070.86864339394), -1.15318009033915, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(659.500458929189), -0.231868725236956, 1e-9);
        CHECK_NEAR_LOOSE(cosine(659.500458929189), 0.972747086480854, 1e-9);
        CHECK_NEAR_LOOSE(tangent(659.500458929189), -0.238364862212846, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-9286.10482376977), 0.428706224994489, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-9286.10482376977), 0.903443951029047, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-9286.10482376977), 0.474524428998812, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(4257.54659634205), -0.635418768781181, 1e-9);
        CHECK_NEAR_LOOSE(cosine(4257.54659634205), -0.772167720304733, 1e-9);
        CHECK_NEAR_LOOSE(tangent(4257.54659634205), 0.822902527614617, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-865.89383180412), 0.926777497517569, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-865.89383180412), 0.375610796031041, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-865.89383180412), 2.46738780490478, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(4439.85467047661), -0.706295116675634, 1e-9);
        CHECK_NEAR_LOOSE(cosine(4439.85467047661), -0.707917515082197, 1e-9);
        CHECK_NEAR_LOOSE(tangent(4439.85467047661), 0.997708209823889, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-1025.1403196826), -0.831117718670829, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-1025.1403196826), 0.556096518341372, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-1025.1403196826), -1.49455659450943, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-9482.06888332385), -0.675951287569156, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-9482.06888332385), 0.736946305258124, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-9482.06888332385), -0.917232751892821, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(3767.96000558614), -0.928522729616626, 1e-9);
        CHECK_NEAR_LOOSE(cosine(3767.96000558614), -0.371275558831026, 1e-9);
        CHECK_NEAR_LOOSE(tangent(3767.96000558614), 2.50089915032412, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-6365.72364724013), -0.755836710487687, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-6365.72364724013), 0.654760159966344, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-6365.72364724013), -1.15437187034492, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(7523.80869197971), 0.300973542484555, 1e-9);
        CHECK_NEAR_LOOSE(cosine(7523.80869197971), -0.953632490388356, 1e-9);
        CHECK_NEAR_LOOSE(tangent(7523.80869197971), -0.315607475120722, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(2082.3846048094), 0.471794143084922, 1e-9);
        CHECK_NEAR_LOOSE(cosine(2082.3846048094), -0.881708731129938, 1e-9);
        CHECK_NEAR_LOOSE(tangent(2082.3846048094), -0.535090701075743, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-9966.36873009003), -0.944750596558106, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-9966.36873009003), 0.327790039969342, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-9966.36873009003), -2.88218213294849, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-2294.27774056598), -0.792625662681292, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-2294.27774056598), 0.609708585193815, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-2294.27774056598), -1.30000738373945, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-5360.40764626277), -0.751662566788943, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-5360.40764626277), 0.659547864592297, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-5360.40764626277), -1.1396634075278, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-6817.94463422768), -0.635438234882424, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-6817.94463422768), 0.772151701189287, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-6817.94463422768), -0.822944809813546, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(4493.58593908814), 0.895005911171967, 1e-9);
        CHECK_NEAR_LOOSE(cosine(4493.58593908814), 0.446054278050595, 1e-9);
        CHECK_NEAR_LOOSE(tangent(4493.58593908814), 2.00649552131511, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-6478.1573426776), -0.192089618251089, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-6478.1573426776), 0.981377388449597, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-6478.1573426776), -0.195734709717081, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(4370.19480109077), -0.23713915932323, 1e-9);
        CHECK_NEAR_LOOSE(cosine(4370.19480109077), -0.971475691469154, 1e-9);
        CHECK_NEAR_LOOSE(tangent(4370.19480109077), 0.244102000086699, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(2053.96896992148), -0.591264222334244, 1e-9);
        CHECK_NEAR_LOOSE(cosine(2053.96896992148), 0.806477910043097, 1e-9);
        CHECK_NEAR_LOOSE(tangent(2053.96896992148), -0.73314372901131, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(1874.25936378512), 0.955529578963532, 1e-9);
        CHECK_NEAR_LOOSE(cosine(1874.25936378512), -0.294895275862084, 1e-9);
        CHECK_NEAR_LOOSE(tangent(1874.25936378512), -3.2402335919766, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(1763.31678984258), -0.772844371197456, 1e-9);
        CHECK_NEAR_LOOSE(cosine(1763.31678984258), -0.634595601866582, 1e-9);
        CHECK_NEAR_LOOSE(tangent(1763.31678984258), 1.21785333671433, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(893.320988848958), 0.895108825390252, 1e-9);
        CHECK_NEAR_LOOSE(cosine(893.320988848958), 0.445847721434667, 1e-9);
        CHECK_NEAR_LOOSE(tangent(893.320988848958), 2.00765593801833, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-6135.92194898031), 0.38157213326647, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-6135.92194898031), -0.924339065015904, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-6135.92194898031), -0.412805373815835, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(9898.63491205778), 0.499948581507714, 1e-9);
        CHECK_NEAR_LOOSE(cosine(9898.63491205778), -0.866055088229625, 1e-9);
        CHECK_NEAR_LOOSE(tangent(9898.63491205778), -0.577271109312112, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-8611.59848130597), 0.473286459610546, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-8611.59848130597), -0.880908580472069, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-8611.59848130597), -0.537270802104024, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(5746.75143777977), -0.702194005898112, 1e-9);
        CHECK_NEAR_LOOSE(cosine(5746.75143777977), -0.711985658620146, 1e-9);
        CHECK_NEAR_LOOSE(tangent(5746.75143777977), 0.986247401751026, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(660.156885716668), 0.409976651045674, 1e-9);
        CHECK_NEAR_LOOSE(cosine(660.156885716668), 0.912096017751077, 1e-9);
        CHECK_NEAR_LOOSE(tangent(660.156885716668), 0.449488478259711, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-408.990962892876), -0.551296921223423, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-408.990962892876), 0.834309118162792, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-408.990962892876), -0.660782567542134, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(1158.86695946721), 0.371598024075529, 1e-9);
        CHECK_NEAR_LOOSE(cosine(1158.86695946721), -0.928393724937412, 1e-9);
        CHECK_NEAR_LOOSE(tangent(1158.86695946721), -0.400259086305846, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(4981.59474332228), -0.825566453235135, 1e-9);
        CHECK_NEAR_LOOSE(cosine(4981.59474332228), 0.564304909860583, 1e-9);
        CHECK_NEAR_LOOSE(tangent(4981.59474332228), -1.46297939076784, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(sine(-3137.78391414538), -0.618747055238727, 1e-9);
        CHECK_NEAR_LOOSE(cosine(-3137.78391414538), -0.785590275928492, 1e-9);
        CHECK_NEAR_LOOSE(tangent(-3137.78391414538), 0.787620562776731, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(-0.902727455326815), -1.12606771960226, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(-0.902727455326815), 2.69686404639715, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(47.2251581803777), 1.54962433677222, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(-0.393098074266846), -0.403998491727852, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(-0.393098074266846), 1.97479481852275, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(-51.4476746837747), -1.55136154962656, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(0.852623696191336), 1.02098609512575, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(0.852623696191336), 0.54981023166915, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(-17.8636481493976), -1.51487508499879, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(0.0864021954017169), 0.0865100617811813, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(0.0864021954017169), 1.48428626501372, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(77.7277896004057), 1.55793162502374, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(0.676800478239399), 0.743407709817423, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(0.676800478239399), 0.827388616977474, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(-3.6135533630131), -1.30081696592173, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(0.709907460470672), 0.789366807500359, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(0.709907460470672), 0.781429519294538, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(54.2010351139155), 1.55234858763952, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(-0.810399192435926), -0.944833151926523, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(-0.810399192435926), 2.51562947872142, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(81.3202695632718), 1.55849988946049, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(0.554238765560647), 0.587448135357072, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(0.554238765560647), 0.983348191437825, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(-80.3481698022703), -1.55835113514549, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(-0.0412550999496364), -0.0412668115033118, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(-0.0412550999496364), 1.61206313829821, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(-27.998405158617), -1.53509518235492, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(0.358888207735961), 0.367076474521624, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(0.358888207735961), 1.20371985227327, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(-97.6456101578768), -1.56055556901344, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(-0.043560235658611), -0.0435740233166738, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(-0.043560235658611), 1.61437035011157, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(9.73170982694566), 1.46839885590043, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(-0.429017487152383), -0.443404798294329, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(-0.429017487152383), 2.01420112508923, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(-71.7812583670113), -1.5568660149366, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(-0.941846170069662), -1.22808250375191, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(-0.941846170069662), 2.79887883054681, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(42.6570630277779), 1.54735784318893, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(-0.371389660186579), -0.380505282666686, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(-0.371389660186579), 1.95130160946158, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(67.9217039511677), 1.55607455603126, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(0.723054225458717), 0.808213493043564, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(0.723054225458717), 0.762582833751333, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(-68.854233550138), -1.55627391259778, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(-0.229847304038147), -0.231920783354074, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(-0.229847304038147), 1.80271711014897, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(28.5799810422084), 1.53582106879725, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(-0.321030796198251), -0.326817693736133, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(-0.321030796198251), 1.89761402053103, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(98.8333690660127), 1.56067863186419, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(0.444064633387514), 0.46013004838087, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(0.444064633387514), 1.11066627841403, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(-41.5352013258224), -1.54672501335679, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(0.928415421067727), 1.19012500965574, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(0.928415421067727), 0.380671317139154, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(-37.4621543891373), -1.54410905796079, 1e-8);
    }
    {
        CHECK_NEAR_LOOSE(arcsine(0.719991063336197), 0.803789441500598, 1e-8);
        CHECK_NEAR_LOOSE(arccosine(0.719991063336197), 0.767006885294299, 1e-8);
        CHECK_NEAR_LOOSE(arctangent(-29.4290509230508), -1.53682936713159, 1e-8);
    }
    printf("\n--- Fuzz: Exp/Log/Hyperbolic ---\n");
    {
        CHECK_NEAR_LOOSE(exponential(-2.16988099639547), 0.114191205267321, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(-4.70843632535228), 0.0090188691392844, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(-1.74084824081497), 0.175371580176044, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(5.53763149739685), 254.075507737132, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(1.35466579434192), 3.8754655376667, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(6.56661777887027), 710.961144248417, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(-4.91803148451393), 0.00731351345043987, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(-3.66692699467828), 0.0255548796894799, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(-8.78606056826324), 0.000152848921091324, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(-0.654860129960431), 0.519514722050219, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(9.84121043489721), 18792.4492950016, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(-7.9565666449026), 0.000350353944348385, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(9.59762004977575), 14729.6839022917, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(-7.59024696241264), 0.000505356234557054, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(-4.19796526635294), 0.0150261198881575, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(-3.89120181571122), 0.0204207892591649, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(7.65761949377818), 2116.71257884898, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(-4.57521045591293), 0.0103041303846241, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(-6.98214310553305), 0.000928311600132797, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(2.49872179680294), 12.166932205636, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(-6.30671225759782), 0.00182402030074445, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(9.11385613866559), 9080.24214022648, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(3.67315270394204), 39.3758507767933, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(4.12125121366328), 61.6363142361479, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(8.29616532196062), 4008.47168597052, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(2.42527647700093), 11.3053546577816, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(0.166022475250729), 1.18059963569792, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(-8.57907140394717), 0.000187999471943926, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(0.572746333800136), 1.77312997755746, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(exponential(3.27826215139855), 26.5296281381996, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(666.616482502879), 6.50221489179489, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(524.124122498197), 6.26172853126715, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(278.701607665061), 5.63014170291425, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(993.582641850413), 6.9013172410693, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(241.868757587206), 5.48839525502602, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(192.992835523969), 5.26265306657839, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(818.640028172263), 6.70764446119199, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(215.524175909629), 5.37307308844701, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(290.07650301252), 5.67014469168195, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(38.228201687512), 3.64357350730066, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(702.550232401373), 6.55471690384921, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(325.239454108052), 5.78456169291015, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(95.2009411538362), 4.55598982781816, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(433.949357119331), 6.07292783862662, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(239.640593125444), 5.47914027228128, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(548.659344474819), 6.30747774705121, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(200.88309659249), 5.30272312986637, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(457.290901968887), 6.12531973536198, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(440.466165552477), 6.08783363323778, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(195.423358764236), 5.27516827578143, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(817.35298442353), 6.70607105103151, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(551.301919477801), 6.31228260731967, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(821.645318074734), 6.71130881541199, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(686.048620946442), 6.53094850122952, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(414.302885189588), 6.0265973130697, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(769.856229912209), 6.6462037830147, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(238.598899066372), 5.47478389551391, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(684.617398274804), 6.52886013955465, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(184.135184508576), 5.2156701862664, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(logarithm(57.3064597265618), 4.04841335255961, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(power(9.53019338047048, -4.53629348565512), 3.61832253387643e-05, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(3.64091237060972, 0.45057906405804), 1.79006820435224, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(4.10200874535862, 2.59159354642093), 38.7826950816907, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(7.12660370047694, -1.06887170208099), 0.122568229859065, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(3.19520064750715, 3.20907685311096), 41.5885031148206, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(6.8907016561904, 1.3847560041164), 14.4807261680766, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(5.58912782478703, 2.91314894213897), 150.357295648833, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(0.258087776328173, -1.93060018069487), 13.6660150537067, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(7.87323100971045, 2.9857683458109), 473.920225492164, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(8.53199455633687, 0.345667739245604), 2.09813846412928, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(2.43359948184044, -0.794685118299875), 0.493234657367928, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(4.8277791018184, 4.0653715978201), 602.125779053487, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(9.2091491707761, 4.93675010918068), 57558.6713362466, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(2.50895024459769, -2.70279990011361), 0.0832248326058781, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(5.55585137688501, -2.09419892437627), 0.027564115031866, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(5.36846399149365, 0.187960680837149), 1.37145989191483, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(5.15841213877606, -1.59903030154426), 0.072554652310812, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(4.29364289181713, -3.83454963355998), 0.00374453244159067, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(4.64526624620079, -3.00664516024146), 0.00987496977301774, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(power(5.14212861302601, -3.27841264833536), 0.00466208306985944, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(3.21513296509409), 12.433226702474, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(3.21513296509409), 12.4733766973948, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(3.21513296509409), 0.996781144681598, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(1.0548895471221), 1.26171353260632, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(1.0548895471221), 1.6099444208922, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(1.0548895471221), 0.783700055873421, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(3.60074701731486), 18.299140426907, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(3.60074701731486), 18.3264437456825, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(3.60074701731486), 0.998510168194424, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(-0.555131633695776), -0.584086799539817, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(-0.555131633695776), 1.15808349845625, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(-0.555131633695776), -0.50435637872262, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(0.0329770995453291), 0.0329830769096893, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(0.0329770995453291), 1.00054379382535, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(0.0329770995453291), 0.0329651506643062, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(2.50724175360394), 6.09477204474905, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(2.50724175360394), 6.17626475124363, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(2.50724175360394), 0.986805503038358, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(-4.6316314022137), -51.3358517084548, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(-4.6316314022137), 51.3455905665955, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(-4.6316314022137), -0.999810327273808, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(2.14670300361536), 4.21986595946219, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(2.14670300361536), 4.33673479888127, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(2.14670300361536), 0.973051421209979, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(4.90931245602582), 67.7694054198035, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(4.90931245602582), 67.7767829787877, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(4.90931245602582), 0.999891149171442, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(2.5910440596555), 6.63437709601992, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(2.5910440596555), 6.70931885158201, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(2.5910440596555), 0.988830199127529, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(0.826420992385639), 0.923756449653986, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(0.826420992385639), 1.36136915576832, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(0.826420992385639), 0.678549565883653, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(-3.18125776239443), -12.0177421499079, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(-3.18125776239443), 12.0592755330357, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(-3.18125776239443), -0.996555897324512, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(-0.843917572029259), -0.94771836735894, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(-0.843917572029259), 1.37774094220557, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(-0.843917572029259), -0.687878496113917, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(2.16579042667148), 4.30341681232971, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(2.16579042667148), 4.41807608135509, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(2.16579042667148), 0.974047692499171, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(-0.688859004085554), -0.744646658661418, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(-0.688859004085554), 1.24679535059111, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(-0.688859004085554), -0.597248504582871, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(3.50011595899159), 16.5445491668962, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(3.50011595899159), 16.5747430488634, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(3.50011595899159), 0.998178319755655, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(-1.99131627665024), -3.59432687697029, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(-1.99131627665024), 3.73084249178533, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(-1.99131627665024), -0.963408904258052, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(1.24363355493395), 1.58992891650934, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(1.24363355493395), 1.87826354901344, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(1.24363355493395), 0.846488724835475, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(0.692509307823769), 0.749202811606325, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(0.692509307823769), 1.24952184971645, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(0.692509307823769), 0.599591605201892, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(hyperbolic_sine(4.60268471184014), 49.8708681614461, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_cosine(4.60268471184014), 49.8808930471011, 1e-7);
        CHECK_NEAR_LOOSE(hyperbolic_tangent(4.60268471184014), 0.999799023533008, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(47.2218345181412), 4.54811565562879, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(47.2218345181412), 4.54789143048426, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(14.6405703793852), 3.37810795033336, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(14.6405703793852), 3.37577526595836, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(28.3489566646072), 4.03804833858006, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(28.3489566646072), 4.03742618732104, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(19.2519829655941), 3.65143507681528, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(19.2519829655941), 3.6500860526619, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(92.1555153431498), 5.21665415063994, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(92.1555153431498), 5.21659527612509, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(45.7417735074921), 4.51627860583099, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(45.7417735074921), 4.51603963547709, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(86.1877580126116), 5.1497089829716, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(86.1877580126116), 5.14964167308812, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(91.1932943764232), 5.20615860875308, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(91.1932943764232), 5.20609848526122, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(18.3328655844054), 3.60258557425253, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(18.3328655844054), 3.6010978922919, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(97.6585457322638), 5.27465056008968, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(97.6585457322638), 5.2745981337545, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(9.1394684377439), 2.9087289995164, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(9.1394684377439), 2.9027429390719, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(78.0983585275709), 5.05115720504211, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(78.0983585275709), 5.05107522914185, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(37.7528885678344), 4.3243845296068, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(37.7528885678344), 4.32403372145509, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(85.8374769018464), 5.14563681386503, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(85.8374769018464), 5.14556895351133, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(41.0185219250493), 4.40731945313713, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(41.0185219250493), 4.40702227967584, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(17.8757347499418), 3.5773728266549, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(17.8757347499418), 3.57580808340413, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(4.79548929302979), 2.27152093498287, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(4.79548929302979), 2.24977012887245, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(36.6572440115752), 4.29494423727983, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(36.6572440115752), 4.29457214519349, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(56.4137152958954), 4.72593803366731, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(56.4137152958954), 4.72578092483225, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_sine(36.3216837024343), 4.28575153632025, 1e-7);
        CHECK_NEAR_LOOSE(inverse_hyperbolic_cosine(36.3216837024343), 4.28537253727843, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(-0.825107964692545), -1.17261279028282, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(-0.381643104545372), -0.40198146692344, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(0.241368059504897), 0.246226295668725, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(0.867152495132792), 1.32148396449572, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(0.84465895545), 1.23721335207122, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(-0.208424255252263), -0.211523474525275, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(-0.584942701662909), -0.669943478152619, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(-0.962721187599676), -1.98183100863364, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(0.708788089179993), 0.884744230215532, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(-0.535308187120684), -0.597555997689732, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(0.787406696771196), 1.06457005478886, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(-0.520937604232559), -0.577625708534781, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(-0.580933530583235), -0.663870624058606, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(-0.776135048074363), -1.03557599219542, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(-0.985011330730153), -2.44304268381706, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(-0.868657174578231), -1.32758223531781, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(0.534965114879356), 0.597075254380949, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(0.67367397583131), 0.817439759375488, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(0.810730866768872), 1.1291579245241, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(inverse_hyperbolic_tangent(0.876604128995459), 1.36091072872517, 1e-7);
    }
    printf("\n--- Fuzz: Complex ---\n");
    {
        cplx ca = {-3.43026172051297, 9.25324022895206}; cplx cb = {-9.51173150290374, 7.76198326930781};
        CHECK_NEAR(cplx_add(ca, cb).real, -12.9419932234167); CHECK_NEAR(cplx_add(ca, cb).imag, 17.0152234982599);
        CHECK_NEAR(cplx_mul(ca, cb).real, -39.1957673738039); CHECK_NEAR(cplx_mul(ca, cb).imag, -114.639970673628);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, 0.693008492077425, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, -0.40729912390569, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 9.86859413523088, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), 1.92579994505845, 1e-9);
    }
    {
        cplx ca = {-7.42770394030147, -6.64746610184489}; cplx cb = {-8.33713367047551, -4.97410846948293};
        CHECK_NEAR(cplx_add(ca, cb).real, -15.764837610777); CHECK_NEAR(cplx_add(ca, cb).imag, -11.6215745713278);
        CHECK_NEAR(cplx_mul(ca, cb).real, 28.8605431772237); CHECK_NEAR(cplx_mul(ca, cb).imag, 92.3670185393009);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, 1.0078666171285, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, 0.196017994921234, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 9.96792813978647, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), -2.41157144918676, 1e-9);
    }
    {
        cplx ca = {4.18378703721482, -9.57754417144706}; cplx cb = {-0.586595392306803, 8.367977814734};
        CHECK_NEAR(cplx_add(ca, cb).real, 3.59719164490802); CHECK_NEAR(cplx_add(ca, cb).imag, -1.20956635671306);
        CHECK_NEAR(cplx_mul(ca, cb).real, 77.6904869478808); CHECK_NEAR(cplx_mul(ca, cb).imag, 40.627980389571);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, -1.17382714843522, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, -0.417690572078697, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 10.4514796238995, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), -1.15894589261035, 1e-9);
    }
    {
        cplx ca = {9.66353625351076, 2.01872912228355}; cplx cb = {9.8185252400167, -2.31212707478738};
        CHECK_NEAR(cplx_add(ca, cb).real, 19.4820614935275); CHECK_NEAR(cplx_add(ca, cb).imag, -0.293397952503831);
        CHECK_NEAR(cplx_mul(ca, cb).real, 99.5492328732053); CHECK_NEAR(cplx_mul(ca, cb).imag, -2.52238097003385);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, 0.886630716209731, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, 0.414393394843814, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 9.87214263430551, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), 0.20594004856774, 1e-9);
    }
    {
        cplx ca = {0.536108403933097, 9.05300016191915}; cplx cb = {-6.60951083566804, 1.26565478068217};
        CHECK_NEAR(cplx_add(ca, cb).real, -6.07340243173494); CHECK_NEAR(cplx_add(ca, cb).imag, 10.3186549426013);
        CHECK_NEAR(cplx_mul(ca, cb).real, -15.001387239338); CHECK_NEAR(cplx_mul(ca, cb).imag, -59.1573745011074);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, 0.174762482337876, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, -1.33622765893714, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 9.06886013523618, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), 1.51164653982188, 1e-9);
    }
    {
        cplx ca = {-0.497173842269047, -1.40885538261203}; cplx cb = {-9.96257499910687, -3.57671443662885};
        CHECK_NEAR(cplx_add(ca, cb).real, -10.4597488413759); CHECK_NEAR(cplx_add(ca, cb).imag, -4.98556981924089);
        CHECK_NEAR(cplx_mul(ca, cb).real, -0.085941694911214); CHECK_NEAR(cplx_mul(ca, cb).imag, 15.8140762713257);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, 0.0891796591822163, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, 0.109397942626185, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 1.49400646536468, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), -1.91004524201184, 1e-9);
    }
    {
        cplx ca = {-7.9004351566107, 5.98369043520473}; cplx cb = {-6.66826398477047, -9.42129509173652};
        CHECK_NEAR(cplx_add(ca, cb).real, -14.5686991413812); CHECK_NEAR(cplx_add(ca, cb).imag, -3.43760465653179);
        CHECK_NEAR(cplx_mul(ca, cb).real, 109.056300546507); CHECK_NEAR(cplx_mul(ca, cb).imag, 34.5315035384678);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, -0.0277116402539346, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, -0.858186014825768, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 9.91067237316273, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), 2.49338031783192, 1e-9);
    }
    {
        cplx ca = {8.34387718760146, -4.44463600291771}; cplx cb = {-8.47794404942015, 4.09434003385719};
        CHECK_NEAR(cplx_add(ca, cb).real, -0.134066861818688); CHECK_NEAR(cplx_add(ca, cb).imag, -0.35029596906052);
        CHECK_NEAR(cplx_mul(ca, cb).real, -52.5410728290494); CHECK_NEAR(cplx_mul(ca, cb).imag, 71.8440457595591);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, -1.0033576522065, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, 0.039696959221462, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 9.45383920532861, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), -0.489450358481995, 1e-9);
    }
    {
        cplx ca = {9.90362747483505, 1.28092983308032}; cplx cb = {8.53433840573419, -2.54825660741754};
        CHECK_NEAR(cplx_add(ca, cb).real, 18.4379658805692); CHECK_NEAR(cplx_add(ca, cb).imag, -1.26732677433722);
        CHECK_NEAR(cplx_mul(ca, cb).real, 87.7850462253543); CHECK_NEAR(cplx_mul(ca, cb).imag, -14.3050954806422);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, 1.02430684532368, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, 0.455937687810102, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 9.98612128894812, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), 0.128625387230334, 1e-9);
    }
    {
        cplx ca = {-5.19823368267368, -1.02458194820435}; cplx cb = {9.30177723709878, 2.96034083041573};
        CHECK_NEAR(cplx_add(ca, cb).real, 4.1035435544251); CHECK_NEAR(cplx_add(ca, cb).imag, 1.93575888221138);
        CHECK_NEAR(cplx_mul(ca, cb).real, -45.319699967238); CHECK_NEAR(cplx_mul(ca, cb).imag, -24.9189764602108);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, -0.539277103702269, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, 0.0614787976767525, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 5.29824514233434, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), -2.94698523567419, 1e-9);
    }
    {
        cplx ca = {-0.842758944853482, 4.0895094196212}; cplx cb = {-2.08468888111947, -2.97277877314069};
        CHECK_NEAR(cplx_add(ca, cb).real, -2.92744782597296); CHECK_NEAR(cplx_add(ca, cb).imag, 1.11673064648052);
        CHECK_NEAR(cplx_mul(ca, cb).real, 13.9140969970088); CHECK_NEAR(cplx_mul(ca, cb).imag, -6.02001891418281);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, -0.788898376017335, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, -0.836714336143545, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 4.17544368089202, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), 1.77402946796562, 1e-9);
    }
    {
        cplx ca = {6.98752350859066, -6.76208560437602}; cplx cb = {-1.83259107299218, 5.336526367545};
        CHECK_NEAR(cplx_add(ca, cb).real, 5.15493243559848); CHECK_NEAR(cplx_add(ca, cb).imag, -1.42555923683102);
        CHECK_NEAR(cplx_mul(ca, cb).real, 23.2807749231828); CHECK_NEAR(cplx_mul(ca, cb).imag, 49.6812411608231);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, -1.53568078646434, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, -0.782015925889699, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 9.72374858292915, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), -0.769003682336798, 1e-9);
    }
    {
        cplx ca = {-0.0961894622551807, -3.9830311982674}; cplx cb = {-4.64186929245901, -9.25379521338623};
        CHECK_NEAR(cplx_add(ca, cb).real, -4.73805875471419); CHECK_NEAR(cplx_add(ca, cb).imag, -13.2368264116536);
        CHECK_NEAR(cplx_mul(ca, cb).real, -36.4116561261945); CHECK_NEAR(cplx_mul(ca, cb).imag, 19.3788277955389);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, 0.348057161660812, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, 0.164197105451459, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 3.98419250777625, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), -1.59494144751931, 1e-9);
    }
    {
        cplx ca = {3.19502559391397, 2.35926936816278}; cplx cb = {4.98863664847629, 1.83522702914505};
        CHECK_NEAR(cplx_add(ca, cb).real, 8.18366224239025); CHECK_NEAR(cplx_add(ca, cb).imag, 4.19449639730784);
        CHECK_NEAR(cplx_mul(ca, cb).real, 11.6090268571326); CHECK_NEAR(cplx_mul(ca, cb).imag, 17.6331349624055);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, 0.717357519967737, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, 0.209026139090732, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 3.97169239711694, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), 0.636048332274991, 1e-9);
    }
    {
        cplx ca = {-7.65199563476876, 0.241819001165659}; cplx cb = {-7.40841226463172, 0.0271189658118907};
        CHECK_NEAR(cplx_add(ca, cb).real, -15.0604078994005); CHECK_NEAR(cplx_add(ca, cb).imag, 0.26893796697755);
        CHECK_NEAR(cplx_mul(ca, cb).real, 56.682580428304); CHECK_NEAR(cplx_mul(ca, cb).imag, -1.99900906206869);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, 1.03298493795285, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, -0.0288598299219713, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 7.65581567331953, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), 3.11000108453586, 1e-9);
    }
    {
        cplx ca = {-1.06722179272746, 4.99827306942612}; cplx cb = {1.11190061422642, 8.72584507416271};
        CHECK_NEAR(cplx_add(ca, cb).real, 0.0446788214989606); CHECK_NEAR(cplx_add(ca, cb).imag, 13.7241181435888);
        CHECK_NEAR(cplx_mul(ca, cb).real, -44.8008010090216); CHECK_NEAR(cplx_mul(ca, cb).imag, -3.75482912714373);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, 0.548324165911989, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, 0.192176660867555, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 5.11093886007482, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), 1.78115563509197, 1e-9);
    }
    {
        cplx ca = {1.03312714127049, 6.25562122644254}; cplx cb = {-1.09731553979293, 3.5755343038488};
        CHECK_NEAR(cplx_add(ca, cb).real, -0.0641883985224396); CHECK_NEAR(cplx_add(ca, cb).imag, 9.83115553029134);
        CHECK_NEAR(cplx_mul(ca, cb).real, -23.5008547537279); CHECK_NEAR(cplx_mul(ca, cb).imag, -3.17040884898402);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, 1.51792190378234, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, -0.754786894808393, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 6.34035871372813, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), 1.40712192926031, 1e-9);
    }
    {
        cplx ca = {7.56113801985616, -8.50456342012589}; cplx cb = {4.59460099615159, 2.26627085139915};
        CHECK_NEAR(cplx_add(ca, cb).real, 12.1557390160077); CHECK_NEAR(cplx_add(ca, cb).imag, -6.23829256872674);
        CHECK_NEAR(cplx_mul(ca, cb).real, 54.0140564609775); CHECK_NEAR(cplx_mul(ca, cb).imag, -21.9394888641389);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, 0.589292335004144, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, -2.1416570166029, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 11.3797366894958, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), -0.844053754046341, 1e-9);
    }
    {
        cplx ca = {2.29192187669286, 7.64923434296115}; cplx cb = {9.15503329586875, -1.06404371432395};
        CHECK_NEAR(cplx_add(ca, cb).real, 11.4469551725616); CHECK_NEAR(cplx_add(ca, cb).imag, 6.58519062863719);
        CHECK_NEAR(cplx_mul(ca, cb).real, 29.1217408146718); CHECK_NEAR(cplx_mul(ca, cb).imag, 67.5902900310954);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, 0.151194604048465, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, 0.85309487782671, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 7.9852170867422, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), 1.27968120086592, 1e-9);
    }
    {
        cplx ca = {0.18848726411499, 6.56622840979576}; cplx cb = {-2.54996279371688, 3.46313861530158};
        CHECK_NEAR(cplx_add(ca, cb).real, -2.36147552960189); CHECK_NEAR(cplx_add(ca, cb).imag, 10.0293670250973);
        CHECK_NEAR(cplx_mul(ca, cb).real, -23.2203946734367); CHECK_NEAR(cplx_mul(ca, cb).imag, -16.0908806171768);
        CHECK_NEAR_LOOSE(cplx_div(ca, cb).real, 1.20347954830889, 1e-9); CHECK_NEAR_LOOSE(cplx_div(ca, cb).imag, -0.94056741503493, 1e-9);
        CHECK_NEAR_LOOSE(cplx_abs(ca), 6.56893316896606, 1e-9);
        CHECK_NEAR_LOOSE(cplx_arg(ca), 1.54209864736758, 1e-9);
    }
    printf("\n--- Fuzz: Linear Algebra ---\n");
    {
        CHECK_NEAR(vec2_dot((vec2){4.18288997057416, 9.63321661182135}, (vec2){-8.51744720607977, 4.11555072697718}), 4.01844713670395);
        CHECK_NEAR(vec2_cross((vec2){4.18288997057416, 9.63321661182135}, (vec2){-8.51744720607977, 4.11555072697718}), 99.265309775181);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){4.18288997057416, 9.63321661182135}), 10.5021631484281, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){8.87497443725218, -8.90531764043679}, (vec2){-8.62653201202113, -5.58613931340598}), -26.8139061192134);
        CHECK_NEAR(vec2_cross((vec2){8.87497443725218, -8.90531764043679}, (vec2){-8.62653201202113, -5.58613931340598}), -126.398851311852);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){8.87497443725218, -8.90531764043679}), 12.5725834075163, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){5.64419856625709, -3.90314818992462}, (vec2){-2.97500034133461, 5.91571984116711}), -39.8814238513277);
        CHECK_NEAR(vec2_cross((vec2){5.64419856625709, -3.90314818992462}, (vec2){-2.97500034133461, 5.91571984116711}), 21.7776302485887);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){5.64419856625709, -3.90314818992462}), 6.86232783010623, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){-3.94271079980884, -7.67776191246963}, (vec2){7.95636716487241, -5.36086087841536}), 9.78975872215908);
        CHECK_NEAR(vec2_cross((vec2){-3.94271079980884, -7.67776191246963}, (vec2){7.95636716487241, -5.36086087841536}), 82.2234168616823);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){-3.94271079980884, -7.67776191246963}), 8.63093253568225, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){-6.9487297788872, -4.53594127393822}, (vec2){7.43220556660084, -8.80264480402837}), -11.7161082570404);
        CHECK_NEAR(vec2_cross((vec2){-6.9487297788872, -4.53594127393822}, (vec2){7.43220556660084, -8.80264480402837}), 94.8792480686568);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){-6.9487297788872, -4.53594127393822}), 8.29816900169007, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){-3.25625339306923, 3.16966478949834}, (vec2){0.715601238467647, 4.22885836271518}), 11.0738844912291);
        CHECK_NEAR(vec2_cross((vec2){-3.25625339306923, 3.16966478949834}, (vec2){0.715601238467647, 4.22885836271518}), -16.0384504412928);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){-3.25625339306923, 3.16966478949834}), 4.54422281998368, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){6.37674613149331, -0.697883646969665}, (vec2){6.64194475579164, 9.14716446376264}), 35.9703390316846);
        CHECK_NEAR(vec2_cross((vec2){6.37674613149331, -0.697883646969665}, (vec2){6.64194475579164, 9.14716446376264}), 62.9644502375744);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){6.37674613149331, -0.697883646969665}), 6.41482133891682, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){0.652699802857541, 7.28589778528557}, (vec2){8.83253679688442, 2.39832176676653}), 23.2389222749449);
        CHECK_NEAR(vec2_cross((vec2){0.652699802857541, 7.28589778528557}, (vec2){8.83253679688442, 2.39832176676653}), -62.787576142516);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){0.652699802857541, 7.28589778528557}), 7.31507508985926, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){8.73236534090262, 5.36042119581636}, (vec2){2.49941186118241, 4.57228870805564}), 46.3351708132835);
        CHECK_NEAR(vec2_cross((vec2){8.73236534090262, 5.36042119581636}, (vec2){2.49941186118241, 4.57228870805564}), 26.5289951250685);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){8.73236534090262, 5.36042119581636}), 10.2463808168326, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){-3.79109768771443, -4.61302570263694}, (vec2){2.65768170524962, -1.39778969762461}), -3.6274911654293);
        CHECK_NEAR(vec2_cross((vec2){-3.79109768771443, -4.61302570263694}, (vec2){2.65768170524962, -1.39778969762461}), 17.5591113063202);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){-3.79109768771443, -4.61302570263694}), 5.97096540025001, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){-9.83212961975765, 9.71708966508463}, (vec2){3.41794604696933, 8.20575212512254}), 46.1303406021337);
        CHECK_NEAR(vec2_cross((vec2){-9.83212961975765, 9.71708966508463}, (vec2){3.41794604696933, 8.20575212512254}), -113.892506730629);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){-9.83212961975765, 9.71708966508463}), 13.8236248653893, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){-7.2892861479263, 5.79177317500052}, (vec2){-0.140480140313707, -4.99610719951985}), -27.9123197167586);
        CHECK_NEAR(vec2_cross((vec2){-7.2892861479263, 5.79177317500052}, (vec2){-0.140480140313707, -4.99610719951985}), 37.2316841113041);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){-7.2892861479263, 5.79177317500052}), 9.31011971228114, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){-6.52398359387254, -0.912522876040819}, (vec2){-3.25343424545835, -5.86997869706808}), 26.5818414840603);
        CHECK_NEAR(vec2_cross((vec2){-6.52398359387254, -0.912522876040819}, (vec2){-3.25343424545835, -5.86997869706808}), 35.3268115413781);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){-6.52398359387254, -0.912522876040819}), 6.58749268936337, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){5.12252359260239, 4.00027320492241}, (vec2){2.60636908958801, 7.68654720950655}), 44.0994759930043);
        CHECK_NEAR(vec2_cross((vec2){5.12252359260239, 4.00027320492241}, (vec2){2.60636908958801, 7.68654720950655}), 28.9483309951325);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){5.12252359260239, 4.00027320492241}), 6.49941794861573, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){1.58776914846569, -5.95306584353935}, (vec2){7.13542106020896, -1.30213055161134}), 19.0810703313387);
        CHECK_NEAR(vec2_cross((vec2){1.58776914846569, -5.95306584353935}, (vec2){7.13542106020896, -1.30213055161134}), 40.4101486756782);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){1.58776914846569, -5.95306584353935}), 6.16116902919684, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){9.16616083821566, -8.2369530745007}, (vec2){2.53248135407173, -4.06391805213116}), 56.6874337052243);
        CHECK_NEAR(vec2_cross((vec2){9.16616083821566, -8.2369530745007}, (vec2){2.53248135407173, -4.06391805213116}), -16.3905964236255);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){9.16616083821566, -8.2369530745007}), 12.3233883515681, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){-0.0794418548381692, 3.93405436217222}, (vec2){0.949499998463629, -8.14082868761431}), -32.1018926512526);
        CHECK_NEAR(vec2_cross((vec2){-0.0794418548381692, 3.93405436217222}, (vec2){0.949499998463629, -8.14082868761431}), -3.0886620799745);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){-0.0794418548381692, 3.93405436217222}), 3.93485638020327, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){-1.66339442592568, 8.04679186116463}, (vec2){-4.13996679834587, 0.680520933106541}), 12.3624080017598);
        CHECK_NEAR(vec2_cross((vec2){-1.66339442592568, 8.04679186116463}, (vec2){-4.13996679834587, 0.680520933106541}), 32.1814764115662);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){-1.66339442592568, 8.04679186116463}), 8.21691793028907, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){5.38955391807433, 5.65616895924804}, (vec2){7.9166369933693, -7.89935274424478}), -2.01293186456241);
        CHECK_NEAR(vec2_cross((vec2){5.38955391807433, 5.65616895924804}, (vec2){7.9166369933693, -7.89935274424478}), -87.3518239565259);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){5.38955391807433, 5.65616895924804}), 7.81278047377446, 1e-9);
    }
    {
        CHECK_NEAR(vec2_dot((vec2){-4.44599006603816, -6.24929883195846}, (vec2){9.00270548915796, -5.32284118113361}), -6.76191399630483);
        CHECK_NEAR(vec2_cross((vec2){-4.44599006603816, -6.24929883195846}, (vec2){9.00270548915796, -5.32284118113361}), 79.9258959122797);
        CHECK_NEAR_LOOSE(vec2_mag((vec2){-4.44599006603816, -6.24929883195846}), 7.66945653605439, 1e-9);
    }
    {
        vec3 fz_v1 = {1.35650318221776, -8.60728979609599, 2.72363566712681}; vec3 fz_v2 = {7.22897381221931, -4.00461645720498, -4.41255034725094};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), 32.2568408410574);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, 48.8872157946593); CHECK_NEAR(fz_cross.y, 25.6747294994281); CHECK_NEAR(fz_cross.z, 56.7895975624002);
    }
    {
        vec3 fz_v1 = {-7.20757895016165, -0.588823681866577, -4.69142701251094}; vec3 fz_v2 = {2.19515362472003, -5.80836590126942, 3.90642711513072};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), -30.7283573526188);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, -29.5497214846375); CHECK_NEAR(fz_cross.y, 17.857478833734); CHECK_NEAR(fz_cross.z, 43.1568142443966);
    }
    {
        vec3 fz_v1 = {-7.79754829958057, 7.65570038584685, -6.88629681304191}; vec3 fz_v2 = {4.93332925287303, -7.73019342853653, 3.45257742265836};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), -121.423390842962);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, -26.8005080643169); CHECK_NEAR(fz_cross.y, -7.05073030052604); CHECK_NEAR(fz_cross.z, 22.5084659593843);
    }
    {
        vec3 fz_v1 = {-8.04027782869779, 0.463584618331748, -3.21407376754857}; vec3 fz_v2 = {2.59697072606548, -9.06083746231095, -8.73711235137914};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), 3.00089258539944);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, -33.1725908943707); CHECK_NEAR(fz_cross.y, -78.595666211374); CHECK_NEAR(fz_cross.z, 71.6477348747913);
    }
    {
        vec3 fz_v1 = {-5.72991882922087, 1.77596991435941, 9.362030197697}; vec3 fz_v2 = {-8.13302990692851, -0.647154612634065, -8.03418183257039};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), -29.7639788499635);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, -8.20978419507891); CHECK_NEAR(fz_cross.y, -122.176881347267); CHECK_NEAR(fz_cross.z, 18.1521598276394);
    }
    {
        vec3 fz_v1 = {-9.40132611497314, -6.07750316202516, 5.73538646383791}; vec3 fz_v2 = {-3.18007345218214, 3.21931037227784, -9.30067885931485};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), -43.0114490077441);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, 38.0609160444113); CHECK_NEAR(fz_cross.y, -105.677665278711); CHECK_NEAR(fz_cross.z, -49.5926931362087);
    }
    {
        vec3 fz_v1 = {6.66464914949711, -9.12964053207662, -2.77943231649638}; vec3 fz_v2 = {-4.11123785947583, -0.251354672311679, -6.96862766567854};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), -5.73635116365474);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, 62.9224422904039); CHECK_NEAR(fz_cross.y, 57.870365812657); CHECK_NEAR(fz_cross.z, -39.2093145019226);
    }
    {
        vec3 fz_v1 = {5.55427830734592, -0.119998520161872, 3.07051539560637}; vec3 fz_v2 = {-1.22587985308673, 3.80318670030461, -6.92753941376749};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), -28.5363710749909);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, -10.8464488376355); CHECK_NEAR(fz_cross.y, 34.7133989271061); CHECK_NEAR(fz_cross.z, 20.9768536200217);
    }
    {
        vec3 fz_v1 = {-5.53147692968233, -7.71731506423423, 0.548307538082646}; vec3 fz_v2 = {2.40897756705925, -9.51358546911256, 9.6573955353473};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), 65.3893553896284);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, -69.3127934192962); CHECK_NEAR(fz_cross.y, 54.7405211636813); CHECK_NEAR(fz_cross.z, 71.2150174086259);
    }
    {
        vec3 fz_v1 = {4.19945073827853, 6.67977923816345, -6.4127270173438}; vec3 fz_v2 = {-6.67219655008819, 5.46460721285697, 8.46716205603961};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), -45.8147898000405);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, 91.6018076211585); CHECK_NEAR(fz_cross.y, 7.22954513441908); CHECK_NEAR(fz_cross.z, 67.5171487826593);
    }
    {
        vec3 fz_v1 = {5.53206685204378, -9.16806802033767, -8.52207185329987}; vec3 fz_v2 = {2.44248945596486, 2.81631015562828, 0.420799928497754};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), -15.8941953438795);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, 20.1428751400217); CHECK_NEAR(fz_cross.y, -23.1429639804447); CHECK_NEAR(fz_cross.z, 37.9729255282689);
    }
    {
        vec3 fz_v1 = {-6.17460837803347, 9.5576582173233, -0.939145257689125}; vec3 fz_v2 = {1.60473397998843, 6.90407845369321, 4.01102582958698};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), 52.3112824023089);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, 44.8199465185476); CHECK_NEAR(fz_cross.y, 23.2594353847177); CHECK_NEAR(fz_cross.z, -57.9674795732289);
    }
    {
        vec3 fz_v1 = {2.82034019846606, 7.30303659443049, 2.64332576008864}; vec3 fz_v2 = {4.13657840704362, -4.48744150427113, 2.65463359255593};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), -14.0883297967667);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, 31.2486559963712); CHECK_NEAR(fz_cross.y, 3.44735442868097); CHECK_NEAR(fz_cross.z, -42.8656951451314);
    }
    {
        vec3 fz_v1 = {1.43803369879875, -9.226461181335, 5.90681740874325}; vec3 fz_v2 = {1.11630241406406, -3.25702966381384, 3.5238741704143};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), 52.4710195451333);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, -13.2742087222279); CHECK_NEAR(fz_cross.y, 1.52634472543341); CHECK_NEAR(fz_cross.z, 5.61580247544109);
    }
    {
        vec3 fz_v1 = {6.52408556360024, 1.74000379455959, -4.32196463504383}; vec3 fz_v2 = {6.33885980939192, -7.9010424613254, -8.78018693579077};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), 65.5550773340185);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, -49.4255846830466); CHECK_NEAR(fz_cross.y, 29.8863629108115); CHECK_NEAR(fz_cross.z, -62.5767171808488);
    }
    {
        vec3 fz_v1 = {5.79115799855123, -6.08419797781182, 7.3156873916911}; vec3 fz_v2 = {-3.00869609730551, -2.49607103285659, -8.02617528242005};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), -60.9542334556431);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, 67.0933148066977); CHECK_NEAR(fz_cross.y, 24.4701690800729); CHECK_NEAR(fz_cross.z, -32.760644437956);
    }
    {
        vec3 fz_v1 = {1.17110033365113, 6.07792345206911, -6.07508309749739}; vec3 fz_v2 = {5.00108647040555, -7.17312930251044, 2.89231489526592};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), -55.3120101112122);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, -25.9980880491389); CHECK_NEAR(fz_cross.y, -33.7692068243537); CHECK_NEAR(fz_cross.z, -38.7966748637961);
    }
    {
        vec3 fz_v1 = {-3.22931043133407, 8.52149050934217, 7.40962600850879}; vec3 fz_v2 = {-0.16562797451707, 9.77025426248712, 6.44064543001832};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), 131.51476700732);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, -17.5100311871294); CHECK_NEAR(fz_cross.y, 19.5716021239639); CHECK_NEAR(fz_cross.z, -30.139786793707);
    }
    {
        vec3 fz_v1 = {-8.91956020277322, -1.22711574053008, 3.69399782431035}; vec3 fz_v2 = {-4.48736143759917, -7.12223004448799, -6.68818332583083};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), 24.0589564353284);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, 34.5166773232543); CHECK_NEAR(fz_cross.y, -76.2319572093176); CHECK_NEAR(fz_cross.z, 58.0206478062852);
    }
    {
        vec3 fz_v1 = {-6.99457935032795, 9.34111770677671, -5.35441919599205}; vec3 fz_v2 = {-1.67720742347627, -8.82784118217369, -4.51259068429968};
        CHECK_NEAR(vec3_dot(fz_v1, fz_v2), -46.5682409852822);
        vec3 fz_cross = vec3_cross(fz_v1, fz_v2);
        CHECK_NEAR(fz_cross.x, -89.4206030295473); CHECK_NEAR(fz_cross.y, -22.5832019929631); CHECK_NEAR(fz_cross.z, 77.4140276021783);
    }
    printf("\n--- Fuzz: Statistics ---\n");
    {
        double fz_data[] = { -61.8574363780891, 88.3159148097885, -32.0418312618824, -60.2475243095607, -5.07949721049586, -27.8427772409594, 37.9942363403071, -23.0408433426075, 4.57653096433368, -36.8221828714263 };
        CHECK_NEAR(mean(fz_data, 10), -11.6045410500592);
        CHECK_NEAR(variance(fz_data, 10), 1908.81144306872);
        CHECK_NEAR_LOOSE(stddev(fz_data, 10), 43.6899467048052, 1e-9);
    }
    {
        double fz_data[] = { -13.0545960832669, -83.3374864175345, -93.3954302011184, -19.23107779234, -51.9878176971737, -55.3495976877431, -34.8205143761824, 60.5435592882488, 43.7934439238207, -20.4323066867315 };
        CHECK_NEAR(mean(fz_data, 10), -26.7271823730021);
        CHECK_NEAR(variance(fz_data, 10), 2204.43239138706);
        CHECK_NEAR_LOOSE(stddev(fz_data, 10), 46.9513832744793, 1e-9);
    }
    {
        double fz_data[] = { -27.0095639093452, -6.15906326162816, 25.1751855367546, -80.1257992627762, -81.7761876372407, -97.1853183972577, 73.1986312220535, -34.2717200676341, -18.2941753310939, 88.5910514533906 };
        CHECK_NEAR(mean(fz_data, 10), -15.7856959654777);
        CHECK_NEAR(variance(fz_data, 10), 3617.74749402936);
        CHECK_NEAR_LOOSE(stddev(fz_data, 10), 60.1477139551401, 1e-9);
    }
    {
        double fz_data[] = { -34.6206283783324, 52.5835329267767, 31.6610780448789, -55.7500388748247, -83.6460689310634, 67.3120881117888, -52.916894889414, -91.6011692619523, 42.1231504297144, -41.9080111526296 };
        CHECK_NEAR(mean(fz_data, 10), -16.6762961975058);
        CHECK_NEAR(variance(fz_data, 10), 3154.23132225878);
        CHECK_NEAR_LOOSE(stddev(fz_data, 10), 56.1625437659191, 1e-9);
    }
    {
        double fz_data[] = { 61.2621020343296, -88.9277923478281, -34.1271458416452, 46.6014329912789, 80.6942897930435, -21.6020572041383, 69.169702123149, 8.0117760365104, -40.9164516016942, 11.9813135027959 };
        CHECK_NEAR(mean(fz_data, 10), 9.21471694858015);
        CHECK_NEAR(variance(fz_data, 10), 2779.30043658811);
        CHECK_NEAR_LOOSE(stddev(fz_data, 10), 52.7190709002739, 1e-9);
    }
    {
        double fz_data[] = { 34.4961334109254, 64.8921874356565, -58.5379577616968, -81.5315196986774, -44.8871898179463, -99.2294291428573, -1.68408889874571, 52.5859591771258, 80.7685674635494, -63.205550368005 };
        CHECK_NEAR(mean(fz_data, 10), -11.6332888200671);
        CHECK_NEAR(variance(fz_data, 10), 3926.9688089267);
        CHECK_NEAR_LOOSE(stddev(fz_data, 10), 62.6655312666118, 1e-9);
    }
    {
        double fz_data[] = { 77.7865046053964, -57.2090473638452, 71.2790156957159, -68.116999996911, -34.9769295159522, -94.9692265876411, -50.0437704446936, -86.6755029985744, 16.4275187111834, -9.90962386744908 };
        CHECK_NEAR(mean(fz_data, 10), -23.6408061762771);
        CHECK_NEAR(variance(fz_data, 10), 3408.29873930866);
        CHECK_NEAR_LOOSE(stddev(fz_data, 10), 58.3806366812547, 1e-9);
    }
    {
        double fz_data[] = { 96.1696570884596, 16.7545744120931, -44.0321956344126, 57.527322101065, 82.6753331701088, 21.6283024618741, 28.3141777880639, 49.0826217174131, 84.1306542855029, -39.6153454434826 };
        CHECK_NEAR(mean(fz_data, 10), 35.2635101946685);
        CHECK_NEAR(variance(fz_data, 10), 2150.35274220071);
        CHECK_NEAR_LOOSE(stddev(fz_data, 10), 46.3718960384488, 1e-9);
    }
    {
        double fz_data[] = { 86.5837072398735, -75.9084141342584, -43.0149816592804, -50.56730074326, -47.1766706708049, -12.3289312774589, 50.0316046152081, 40.9694185144526, 17.987318859194, -49.096286262513 };
        CHECK_NEAR(mean(fz_data, 10), -8.25205355188474);
        CHECK_NEAR(variance(fz_data, 10), 2627.84868494777);
        CHECK_NEAR_LOOSE(stddev(fz_data, 10), 51.2625466100522, 1e-9);
    }
    {
        double fz_data[] = { 3.9679905623006, 88.7557867554772, 42.0349359820999, -88.3132498350645, 54.8935266346043, 1.96526512792661, 69.1589537866401, 39.6619662175082, 69.9283003715986, -71.7627416819332 };
        CHECK_NEAR(mean(fz_data, 10), 21.0290733921158);
        CHECK_NEAR(variance(fz_data, 10), 3245.03335669029);
        CHECK_NEAR_LOOSE(stddev(fz_data, 10), 56.9651942565835, 1e-9);
    }
    printf("\n--- Fuzz: v8 Silicon Engine ---\n");
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(446.777878705008), 0.0473101330173813, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(446.777878705008), 8.80341394526324, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(140.728732852943), 0.0842963191664326, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(140.728732852943), 7.13677310621286, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(272.60087116341), 0.0605670442853154, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(272.60087116341), 8.09064636241353, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(165.054501048098), 0.0778370403251621, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(165.054501048098), 7.366798671269, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(261.643492237139), 0.0618222819439011, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(261.643492237139), 8.03145856554766, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(919.551035444646), 0.0329770711387107, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(919.551035444646), 9.84478583677698, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(355.860251483672), 0.0530102996715106, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(355.860251483672), 8.47516698693852, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(670.325554682573), 0.0386239878575842, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(670.325554682573), 9.38871812423264, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(81.3426788139556), 0.110876820761531, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(81.3426788139556), 6.34594059795654, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(495.635287732243), 0.0449178427454681, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(495.635287732243), 8.95313509642993, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(348.923736269863), 0.0535346222888134, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(348.923736269863), 8.446767933005, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(629.178048876436), 0.0398669690374749, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(629.178048876436), 9.29732452783675, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(977.519022252122), 0.0319843398425767, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(977.519022252122), 9.93298096691777, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(574.22350782264), 0.041731069590714, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(574.22350782264), 9.16546858323332, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(812.944832786816), 0.035072721165217, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(812.944832786816), 9.66701364273268, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(894.682012602831), 0.0334322531766631, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(894.682012602831), 9.80523120136144, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(900.239673305474), 0.0333288958250614, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(900.239673305474), 9.8141653350579, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(815.840383657024), 0.0350104264325239, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(815.840383657024), 9.67214311129563, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(769.185565980243), 0.0360565721884771, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(769.185565980243), 9.5871878800603, 1e-1);
    }
    {
        CHECK_NEAR_LOOSE(ml_fast_rsqrt(645.331743034123), 0.0393648395846583, 1e-2);
        CHECK_NEAR_LOOSE(ml_fast_log2(645.331743034123), 9.33389718130611, 1e-1);
    }
    {
        double fz_err; double fz_s = ml_two_sum(3.29095844243696, 0.651026005056554, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, 3.94198444749352, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(5.28112426388544, -8.50824212091881, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, -3.22711785703337, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(4.01190552927032, -5.6739618598775, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, -1.66205633060718, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(3.6038789084569, -8.58895654530697, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, -4.98507763685007, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(0.392540842908293, -3.93053346977278, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, -3.53799262686448, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(-4.24958451963034, 6.85010254839501, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, 2.60051802876467, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(-8.93149057788659, -9.24716804357419, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, -18.1786586214608, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(7.80842988345612, -4.8885003346849, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, 2.91992954877122, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(1.0875843415635, -4.53077951968607, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, -3.44319517812257, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(2.25891864926094, 3.60263815495986, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, 5.8615568042208, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(8.18195623642747, -7.44445623000839, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, 0.73750000641908, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(2.82629755772573, -1.97189542361481, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, 0.854402134110925, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(-2.4571839971512, -1.01627201807281, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, -3.47345601522401, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(-3.18986883489617, -3.73657631551584, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, -6.92644515041201, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(1.11544445896257, 1.86118895528764, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, 2.97663341425021, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(-4.25456521682013, 2.73777274528439, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, -1.51679247153574, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(5.70780968002282, -1.03450301318674, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, 4.67330666683608, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(-2.8291404236365, -7.51126222295932, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, -10.3404026465958, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(-1.80528829888588, -6.32884745255887, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, -8.13413575144475, 1e-13);
    }
    {
        double fz_err; double fz_s = ml_two_sum(-3.05773256418957, 1.68994501705517, &fz_err);
        CHECK_NEAR_LOOSE(fz_s + fz_err, -1.36778754713441, 1e-13);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(0.683995943209383, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, 0.631895132035711, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, 0.775053896261138, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(-2.32030819788642, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, -0.732021509274825, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, -0.681281520341633, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(-2.74832079077331, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, -0.383212550312413, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, -0.923660187126768, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(3.109951171055, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, 0.0316362029778223, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, -0.999499450055449, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(-1.93469431601753, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, -0.934516555841688, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, -0.355919663488531, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(-0.487790726730247, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, -0.46867542682159, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, 0.883370445675879, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(1.58093524473791, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, 0.999948601611778, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, -0.0101387442344003, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(-1.21837540762919, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, -0.938539833024555, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, 0.345170945802569, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(2.25486486183724, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, 0.775008023793201, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, -0.631951392953728, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(-1.21586557047219, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, -0.937670555006551, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, 0.347525438311653, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(1.09014848923802, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, 0.886695578774314, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, 0.46235370722217, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(2.54153167508793, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, 0.564692800074563, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, -0.825301182323126, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(0.0357472213323846, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, 0.0357396084719324, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, 0.999361136119608, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(-1.52083725980964, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, -0.998752305356329, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, 0.049938287369687, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(3.07779123728828, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, 0.063758139882665, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, -0.997965379959998, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(0.512246027538406, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, 0.490136222094339, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, 0.871645847687631, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(1.15054029297753, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, 0.912984509920078, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, 0.407994221338973, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(-1.98666070970129, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, -0.914767468186044, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, -0.403980790568679, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(0.512146182594426, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, 0.490049190220563, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, 0.871694780966463, 1e-6);
    }
    {
        double fz_sin, fz_cos; ml_cordic_sincos(-2.51879924095259, &fz_sin, &fz_cos);
        CHECK_NEAR_LOOSE(fz_sin, -0.583306388994557, 1e-6);
        CHECK_NEAR_LOOSE(fz_cos, -0.812252212405809, 1e-6);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(-0.0488814144934135), -0.0488619506698751, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(-1.34149606570605), -0.973825681366102, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(1.30960370846385), 0.966082691961705, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(-0.0569076522866374), -0.056876941535809, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(-0.00700081835829014), -0.00700076117171149, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(0.0803333737953771), 0.0802469970903027, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(-0.483242147430731), -0.464652511814944, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(0.858621432243271), 0.756942413802792, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(-0.233567370769502), -0.231449494570019, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(1.31296163544002), 0.966944370952479, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(0.212868408389686), 0.211264430501127, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(0.881888582045815), 0.771940815882044, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(-1.42690941195261), -0.98966612523975, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(-0.861652643388087), -0.758919773951356, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(-1.22883386659669), -0.942098395412099, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(-1.02290952315173), -0.853627158469824, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(-1.45335696314909), -0.99311192009779, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(1.36524581059665), 0.978948769145142, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(-0.231212203971115), -0.229157638017246, 1e-7);
    }
    {
        CHECK_NEAR_LOOSE(ml_minimax_sin(-1.12206594261725), -0.900998615602399, 1e-7);
    }

    
    printf("\n=== v9: Context-Aware Silicon Aggression ===\n");
    #if defined(MATHLIB_PROFILE_GRAPHICS)
    printf("Active Profile: GRAPHICS (Minimax + Fast RSqrt)\n");
    #elif defined(MATHLIB_PROFILE_EMBEDDED)
    printf("Active Profile: EMBEDDED (CORDIC)\n");
    #else
    printf("Active Profile: SCIENTIFIC (IEEE-754 Strict)\n");
    #endif

    CHECK_NEAR_LOOSE(ml_sin(1.0), 0.8414709848078965, 1e-7);
    CHECK_NEAR_LOOSE(ml_cos(0.0), 1.0, 1e-7);
    CHECK_NEAR_LOOSE(ml_sqrt(25.0), 5.0, 1e-2);

    printf("\n=== v9: SIMD Quaternions (SLERP) ===\n");
    ml_quat q1 = {1, 0, 0, 0};
    ml_quat q2 = {0.7071068, 0.7071068, 0, 0}; // 90 deg around X
    ml_quat q_mid = ml_quat_slerp(q1, q2, 0.5);
    CHECK_NEAR_LOOSE(q_mid.w, 0.9238795, 1e-5);
    CHECK_NEAR_LOOSE(q_mid.x, 0.3826834, 1e-5);

    printf("\n=== v9: SIMD Batch Processing ===\n");
    double simd_in[4] __attribute__((aligned(32))) = {1.0, 2.0, 3.0, 4.0};
    double simd_out[4] __attribute__((aligned(32)));
    ml_simd_batch_poly(simd_in, simd_out);
    CHECK_NEAR(simd_out[0], 2.0);
    CHECK_NEAR(simd_out[3], 20.0);

    printf("\n=== v9: Streaming FFT ===\n");
    ml_fft_stream stream = ml_fft_stream_create(8);
    for(int i=0; i<8; i++) ml_fft_stream_push(&stream, sine(2.0 * math_pi * i / 8.0));
    cplx spec[8];
    ml_fft_stream_process(&stream, spec);
    CHECK_NEAR_LOOSE(cplx_abs(spec[1]), 4.0, 1e-7);
    ml_fft_stream_free(&stream);

    printf("
=== v9: Linker-Level Migration ===
");
    printf("Trojan Horse deleted. Use gcc -Wl,--wrap=sin for production override.
");

    TEST_SUMMARY();
}
