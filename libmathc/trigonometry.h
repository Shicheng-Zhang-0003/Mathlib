#ifndef LIBMATHC_TRIGONOMETRY_H
#define LIBMATHC_TRIGONOMETRY_H

#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define math_pi M_PI
inline float sine (float x);
inline float cosine (float x);
inline float tangent (float x);
inline float sine (float x) {
    while (x > math_pi) {x -= 2 * math_pi;}
    while (x < -math_pi) {x += 2 * math_pi;}
    float result = x;
    for (int step = 3; step < 100; step += 2) {
        float additive = pow (x, (float) (step)) / (tgamma (step + 1));
        if (((step - 1) / 2) % 2 == 1) {result -= additive;}
        else {result += additive;}
    } return result;
    //Domain: (-infinity, infinity)
    //Range: [-1, 1]
    //Period: 2pi
} inline float cosine (float x) {
    while (x > math_pi) {x -= 2 * math_pi;}
    while (x < -math_pi) {x += 2 * math_pi;}
    float result = 1;
    for (int step = 2; step < 100; step += 2) {
        float additive = pow (x, (float) (step)) / (tgamma (step + 1));
        if ((step / 2) % 2 == 1) {result -= additive;}
        else {result += additive;}
    } return result;
    //Domain: (-infinity, infinity)
    //Range: [-1, 1]
    //Period: 2pi
} inline float tangent (float x) {
    while (x > math_pi / 2) {x -= math_pi;}
    while (x < -math_pi / 2) {x += math_pi;}
    float divisor_now = 99;
    float divisor_coefficient_now = 97;
    while (1) {
        if (divisor_coefficient_now == -1) {break;}
        divisor_now = divisor_coefficient_now - (pow (x, 2.0) / divisor_now);
        divisor_coefficient_now -= 2;
    } return x / divisor_now;
    //Domain: (-infinity, infinity), constraint (x != 0.5pi + kpi, k ->= integer)
    //Range (-infinity, infinity)
    //Period: pi
    //Restrict: theta = 90, 270
} inline float cosecant (float x);
inline float secant (float x);
inline float cotangent (float x);
inline float cosecant (float x) {return 1 / sine (x);
    //Domain: (-infinity, infinity), constraint (x != kpi, k ->= integer)
    //Range: (-infinity, -1] <-> [1, infinity)
    //Restriction: Undefined at kpi, k ->= integer
} inline float secant (float x) {return 1 / cosine (x);
    //Domain: (-infinity, infinity), constraint (x != 0.5pi + kpi, k ->= integer)
    //Range: (-infinity, -1] <-> [1, infinity)
    //Restriction: Undefined at 0.5pi + kpi, k ->= integer
} inline float cotangent (float x) {return 1 / tangent (x);
    //Domain: (-infinity, infinity), constraint (x != kpi, k ->= integer)
    //Range (-infinity, infinity)
    //Period: pi
    //Restriction: Undefined at kpi, k ->= integer
} inline float arcsine (float x);
inline float arccosine (float x);
inline float arctangent (float x);
inline float arcsine (float x) {return arctangent (x / sqrt (1 - x * x));
    //Domain: [-1, 1]
    //Range: [-0.5pi, 0.5pi]
    //Restriction: -1 <= x <= 1
} inline float arccosine (float x) {return (math_pi / 2) - arcsine (x);
    //Domain: [-1, 1]
    //Range: [0, pi]
    //Restriction -1 <= x <= 1
} inline float arctangent (float x) {
    if (x > 1) {return (math_pi / 2) - arctangent (1 / x);}
    if (x < -1) {return -(math_pi / 2) - arctangent (1 / x);}
    if (x > 0.5) {return (math_pi / 4) + arctangent ((x - 1) / (x + 1));}
    if (x < -0.5) {return -(math_pi / 4) + arctangent ((x + 1) / (1 - x));}
    float result = x;
    for (int step = 3; step < 100; step += 2) {
        float additive = pow (x, (float) (step)) / step;
        if (((step - 1) / 2) % 2 == 1) {result -= additive;}
        else {result += additive;}
    } return result;
    //Domain: (-infinity, infinity)
    //Range: (-0.5pi, 0.5pi)
    //No Restriction
} inline float arccosecant (float x);
inline float arcsecant (float x);
inline float arccotangent (float x);
inline float arccosecant (float x) {return arcsine (1 / x);
    //Domain: (-infinity, -1] <-> [1, infinity)
    //Range: [-0.5pi, 0] <-> [0, 0.5pi]
    //Restriction: |x| >= 1
} inline float arcsecant (float x) {return arccosine (1 / x);
    //Domain: (-infinity, -1] <-> [1, infinity)
    //Range: [0, 0.5pi) <-> (0.5pi, pi]
    //Restriction: |x| >= 1
} inline float arccotangent (float x) {
    if (x == 0) {return math_pi / 2;}
    if (x > 0) {return arctangent (1 / x);}
    return math_pi + arctangent (1 / x);
    //Domain: (-infinity, infinity)
    //Range: (0, pi)
    //No Restriction
}

#endif
