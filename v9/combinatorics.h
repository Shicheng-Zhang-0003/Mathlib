#ifndef LIBMATHC_COMBINATORICS_H
#define LIBMATHC_COMBINATORICS_H

//Library header file for combinatorics
#include <stdio.h>
inline int factorial (int x);
inline int factorial (int x) {
    if (x < 0) {return (int) (0);}
    if (x == 0) {return (int) (1);}
    int result = (int) (x);
    while ((x - 1) > 0) {result *= (x - 1);
    x -= 1;} return result;
} inline int npr (int n, int r);
inline int npr (int n, int r) {
    if (r < 0 || r > n) {return (int) (0);}
    return factorial (n) / factorial (n - r);
}
inline int ncr (int n, int r);
inline int ncr (int n, int r) {
    if (r < 0 || r > n) {return (int) (0);}
    return factorial (n) / (factorial (n - r) * factorial (r));
}



#endif
