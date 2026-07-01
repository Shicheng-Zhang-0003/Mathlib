# MathLib

> A modern, zero-dependency C99 mathematics and numerical computing library focused on performance, determinism, and portability.

## Overview

MathLib is an ambitious attempt to build a complete numerical computing library in modern C.

Rather than being a simple wrapper around `<math.h>`, MathLib aims to provide a unified collection of mathematical routines ranging from elementary functions to advanced numerical algorithms used in scientific computing, graphics, simulation, signal processing, optimization, and embedded systems.

### Long-term goals

- Deterministic numerical behavior
- Zero runtime heap allocations in hot paths
- Modular architecture
- Configurable performance profiles
- Portable C99 implementation
- Optional SIMD acceleration

> **Status:** Work in Progress. APIs and implementations may change.

## Current Features

- Elementary mathematics
- Trigonometric and hyperbolic functions
- Exponential and logarithmic functions
- Fixed-point arithmetic
- IEEE-754 utilities
- Polynomial evaluation
- CORDIC
- Payne–Hanek argument reduction
- Linear algebra
- Matrix and vector operations
- Quaternions
- Tensor utilities
- FFT
- Statistics
- Numerical integration
- ODE solvers
- Optimization algorithms
- SIMD acceleration
- Testing and benchmarking

## Building

```bash
git clone <repository-url>
cd mathlib
cmake -B build
cmake --build build
```

## Project History

MathLib began as a personal learning project.

- **Version 1** was written entirely by the maintainer.
- **Versions 2–10** relied heavily on AI-assisted development, primarily using **Qwen**, to accelerate implementation, experimentation, and refactoring.

The project's architecture, feature selection, testing, integration, and long-term direction remain directed by the maintainer. This history is documented intentionally for transparency.

## License

Licensed under the **GNU GPL v3.0**.

See `LICENSE` for details.
