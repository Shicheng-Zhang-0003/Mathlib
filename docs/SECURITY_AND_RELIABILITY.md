# Security and Reliability Model

v11S treats numerical correctness as a reliability problem.

Primary protections:
- compiler warnings promoted to errors
- sanitizer builds
- fuzz testing
- boundary testing
- regression preservation

The library avoids undefined behavior and documents unsupported numerical regions.
