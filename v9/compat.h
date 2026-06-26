#ifndef MATHLIB_COMPAT_H
#define MATHLIB_COMPAT_H
#include "profiles.h"

// The Trojan Horse: Override standard C math functions
#define sin(x) ml_sin(x)
#define cos(x) ml_cos(x)
#define sqrt(x) ml_sqrt(x)

#endif
