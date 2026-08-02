# v12A1 Minimax Coefficient Validation Report

Generated: 2026-08-02
Method: np.polyfit on Chebyshev-distributed nodes
Ground truth: mpmath at 80 decimal places

## Results

| Function | Max ULP | Worst Input | Status |
|----------|---------|-------------|--------|
| sin | 17 | -7.85398163397448279e-01 | WARN |
| cos | 11 | -6.12436191620152393e-01 | WARN |
| exp | 54 | -3.46365625329309590e-01 | WARN |
| log | 7 | 1.15909068140762989e-01 | WARN |

## Acceptance Criteria

All functions must achieve <= 5 ULP on their reduced domain.
