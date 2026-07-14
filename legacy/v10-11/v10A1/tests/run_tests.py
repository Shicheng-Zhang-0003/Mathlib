#!/usr/bin/env python3
import subprocess
import sys

tests = ["test_core", "test_trig", "test_linalg", "test_dsp"]
total_failed = 0

print("\n" + "="*50)
print("  MATHLIB V1.0.5 MODULAR CI/CD TEST RUNNER")
print("="*50)

for t in tests:
    print(f"\n--- Running {t} ---")
    try:
        res = subprocess.run([f"./{t}"], capture_output=True, text=True)
        print(res.stdout.strip())
        if res.returncode != 0:
            total_failed += 1
            if res.stderr: print(res.stderr.strip())
    except FileNotFoundError:
        print(f"  [ERROR] Binary ./{t} not found. Did you run 'make test_modular'?")
        total_failed += 1

print("\n" + "="*50)
if total_failed == 0:
    print("  🎉 ALL MODULAR TEST SUITES PASSED 🎉")
else:
    print(f"  ❌ {total_failed} SUITE(S) FAILED")
print("="*50 + "\n")
sys.exit(total_failed)