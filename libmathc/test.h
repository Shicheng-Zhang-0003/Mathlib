#ifndef LIBMATHC_TEST_H
#define LIBMATHC_TEST_H

#include <stdio.h>
#include <math.h>

#define TEST_EPSILON 0.002

extern int tests_passed;
extern int tests_failed;

#define CHECK_NEAR(got,expected) {float diff=fabsf((float)((got)-(expected)));if(diff<TEST_EPSILON){tests_passed++;}else{tests_failed++;printf("FAIL %s:%d  got %f  expected %f  diff %f\n",__FILE__,__LINE__,(float)(got),(float)(expected),diff);}}

#define CHECK_NEAR_LOOSE(got,expected,eps) {float diff=fabsf((float)((got)-(expected)));if(diff<(eps)){tests_passed++;}else{tests_failed++;printf("FAIL %s:%d  got %f  expected %f  diff %f\n",__FILE__,__LINE__,(float)(got),(float)(expected),diff);}}

#define CHECK_INT(got,expected) {if((got)==(expected)){tests_passed++;}else{tests_failed++;printf("FAIL %s:%d  got %d  expected %d\n",__FILE__,__LINE__,(got),(expected));}}

#define TEST_SUMMARY() {printf("\n=== SUMMARY ===\n");printf("passed: %d  failed: %d\n",tests_passed,tests_failed);return tests_failed>0?1:0;}

#endif
