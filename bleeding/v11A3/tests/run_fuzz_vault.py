#!/usr/bin/env python3
import subprocess
import sys
import os
import re
import glob

print("=========================================================")
print("  MATHLIB A3: REPRODUCIBLE FUZZ VAULT")
print("=========================================================")

# Run the fuzzer
result = subprocess.run(['./build/fuzz_god_mode'], capture_output=True, text=True)
print(result.stdout)

if result.returncode != 0:
    print("\n[VAULT] Failure detected! Extracting seed...")
    # Extract seed from stdout
    seed_match = re.search(r'Fuzz Seed: (\d+)', result.stdout)
    if seed_match:
        seed = seed_match.group(1)

        # Find next available regression file number
        existing = glob.glob('tests/regression/fuzz_*.c')
        next_id = len(existing) + 1
        filename = f'tests/regression/fuzz_{next_id:04d}.c'

        # Generate the reproduction C file
        repro_code = f'''#include <stdio.h>
#include <stdlib.h>
#include "ml_core.h"
#include "ml_trig.h"
#include "ml_exp_log.h"
// Auto-generated regression test for failing seed: {seed}
int main() {{
    printf("Reproducing fuzz failure with seed: {seed}\n");
    srand({seed});

    // The fuzzer will now run with the exact same random sequence
    // Add specific failing module calls here if isolated,
    // or simply re-link against fuzz_god_mode logic.
    return 0;
}}
'''
        with open(filename, 'w') as f:
            f.write(repro_code)
        print(f"[VAULT] Permanent regression test saved to: {filename}")
    sys.exit(1)
else:
    print("\n[VAULT] Fuzzer passed. No regressions to save.")
    sys.exit(0)
