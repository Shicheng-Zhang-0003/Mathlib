#ifndef LIBMATHC_TEST_H
#define LIBMATHC_TEST_H

#include <stdio.h>
#include <math.h>
#include <stdint.h>

#define TEST_EPSILON 1e-9

extern int tests_passed;
extern int tests_failed;

#define CHECK_NEAR(got,expected) {double diff=fabs((double)((got)-(expected)));if(diff<TEST_EPSILON){tests_passed++;}else{tests_failed++;printf("FAIL %s:%d  got %.15g  expected %.15g  diff %.15g\n",__FILE__,__LINE__,(double)(got),(double)(expected),diff);}}
#define CHECK_NEAR_LOOSE(got,expected,eps) {double diff=fabs((double)((got)-(expected)));if(diff<(eps)){tests_passed++;}else{tests_failed++;printf("FAIL %s:%d  got %.15g  expected %.15g  diff %.15g\n",__FILE__,__LINE__,(double)(got),(double)(expected),diff);}}
#define CHECK_INT(got,expected) {if((long long)(got)==(long long)(expected)){tests_passed++;}else{tests_failed++;printf("FAIL %s:%d  got %lld  expected %lld\n",__FILE__,__LINE__,(long long)(got),(long long)(expected));}}
#define CHECK_NAN(got) {if((got)!=(got)){tests_passed++;}else{tests_failed++;printf("FAIL %s:%d  got %.15g  expected NaN\n",__FILE__,__LINE__,(double)(got));}}
#define CHECK_INF(got) {if(isinf(got)){tests_passed++;}else{tests_failed++;printf("FAIL %s:%d  got %.15g  expected INF\n",__FILE__,__LINE__,(double)(got));}}
#define CHECK_NEG_INF(got) {if(isinf(got) && (got)<0){tests_passed++;}else{tests_failed++;printf("FAIL %s:%d  got %.15g  expected -INF\n",__FILE__,__LINE__,(double)(got));}}
#define TEST_SUMMARY() {printf("\n=== SUMMARY ===\n");printf("passed: %d  failed: %d\n",tests_passed,tests_failed);return tests_failed>0?1:0;}

#endif
