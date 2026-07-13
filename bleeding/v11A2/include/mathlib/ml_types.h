#ifndef MATHLIB_TYPES_H
#define MATHLIB_TYPES_H

// v11A1 API Standardization: Unified Status Codes
typedef enum {
    ML_SUCCESS          =  0,
    ML_ERR_SINGULAR     = -1,
    ML_ERR_WORKSPACE    = -2,
    ML_ERR_INVALID_ARG  = -3,
    ML_ERR_NAN_INPUT    = -4,
    ML_ERR_INTERNAL     = -99
} ml_status_t;

#endif
