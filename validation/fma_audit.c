#include "ml_exp_log.h"
#include "internal/minimax.h"
double audit_exp(double x) { return ml_exp(x); }
double audit_sin(double x) { return ml_minimax_sin(x); }
