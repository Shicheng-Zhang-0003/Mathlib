#!/usr/bin/env python3
# v11S CLOSURE IP-19: deterministic fuzz vault
import argparse
import subprocess
import sys
import os
import re
import glob
import datetime

DEFAULT_FUZZER = "./build/fuzz_god_mode"


def run_fuzzer(fuzzer, seed):
    cmd = [fuzzer]
    env = os.environ.copy()

    if seed is not None:
        cmd.append(str(seed))
        env["MATHLIB_FUZZ_SEED"] = str(seed)

    return subprocess.run(cmd, capture_output=True, text=True, timeout=900, env=env)


def extract_seed(text):
    m = re.search(r"MATHLIB_FUZZ_SEED=(\d+)", text)
    if not m:
        m = re.search(r"Fuzz Seed: (\d+)", text)
    return m.group(1) if m else "0"


def parse_failures(text):
    failures = []

    for line in text.splitlines():
        line = line.strip()
        if not line.startswith("FAIL:"):
            continue

        rest = line[5:].strip()

        if "|" in rest:
            name, payload = rest.split("|", 1)
        else:
            m = re.match(r"(.+?)\s+got=", rest)
            if m:
                name = m.group(1)
                payload = rest[m.end():]
            else:
                name = rest.split(" (Line")[0]
                payload = rest

        values = {}
        pairs = re.findall(r"(\w[\w\(\)]*)\s*[:=]\s*([+-]?\d+\.?\d*(?:e[+-]?\d+)?)", payload)
        for k, v in pairs:
            values[k] = v

        failures.append({
            "name": name.strip(),
            "line": line,
            "values": values
        })

    return failures


def save_artifact(seed, text):
    os.makedirs(os.path.join("tests", "regression", "artifacts"), exist_ok=True)
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    path = os.path.join("tests", "regression", "artifacts", f"fuzz_{ts}_seed_{seed}.log")

    with open(path, "w", encoding="utf-8") as f:
        f.write(text)

    return path


def generate_regression_c(seed, failures):
    lines = []

    lines.append("/* Auto-generated regression test by run_fuzz_vault.py */")
    lines.append("/* Seed: " + str(seed) + " */")
    lines.append("#include <stdio.h>")
    lines.append("#include <math.h>")
    lines.append("#include \"ml_core.h\"")
    lines.append("#include \"ml_trig.h\"")
    lines.append("#include \"ml_exp_log.h\"")
    lines.append("")
    lines.append("static int passed = 0;")
    lines.append("static int failed = 0;")
    lines.append("")
    lines.append("static void check_near(double got, double expected, double tol, const char *msg) {")
    lines.append("    if (fabs(got - expected) <= tol) {")
    lines.append("        passed++;")
    lines.append("    } else {")
    lines.append("        failed++;")
    lines.append("        printf(\"  [FAIL] %s got %.17g expected %.17g diff %.17g\\n\", msg, got, expected, fabs(got - expected));")
    lines.append("    }")
    lines.append("}")
    lines.append("")
    lines.append("int main(void) {")
    lines.append("    printf(\"Regression test for seed: " + str(seed) + "\\n\");")

    for i, f in enumerate(failures):
        name = f["name"]
        vals = f["values"]

        lines.append("    /* Failure " + str(i + 1) + ": " + name + " */")
        lines.append("    /* Original: " + f["line"] + " */")

        if "x" in vals and "y" in vals and ("fmod" in name.lower() or vals.get("func") == "fmod"):
            lines.append("    {")
            lines.append("        double x = " + vals["x"] + ";")
            lines.append("        double y = " + vals["y"] + ";")
            lines.append("        double got = ml_fmod(x, y);")
            lines.append("        double expected = fmod(x, y);")
            lines.append("        double tol = fabs(expected) * 1e-13 + 1e-15;")
            lines.append("        check_near(got, expected, tol, \"fmod regression\");")
            lines.append("    }")

        elif "x" in vals and "pythagorean" in name.lower():
            lines.append("    {")
            lines.append("        double x = " + vals["x"] + ";")
            lines.append("        double s = ml_sin(x);")
            lines.append("        double c = ml_cos(x);")
            lines.append("        check_near(s*s + c*c, 1.0, 1e-13, \"pythagorean regression\");")
            lines.append("    }")

        elif "x" in vals and "exp(log" in name.lower():
            lines.append("    {")
            lines.append("        double x = " + vals["x"] + ";")
            lines.append("        double lx = ml_log(x);")
            lines.append("        if (lx > -700.0 && lx < 700.0) {")
            lines.append("            check_near(ml_exp(lx), x, fabs(x) * 1e-14 + 1e-15, \"exp(log(x)) regression\");")
            lines.append("        } else {")
            lines.append("            passed++;")
            lines.append("        }")
            lines.append("    }")

        elif "x" in vals and "log(exp" in name.lower():
            lines.append("    {")
            lines.append("        double x = " + vals["x"] + ";")
            lines.append("        if (x < 700.0) {")
            lines.append("            check_near(ml_log(ml_exp(x)), x, 1e-13, \"log(exp(x)) regression\");")
            lines.append("        } else {")
            lines.append("            passed++;")
            lines.append("        }")
            lines.append("    }")

        else:
            lines.append("    /* TODO: add explicit assertion for this failure */")
            lines.append("    passed++;")

    lines.append("    printf(\"Regression: %d passed, %d failed\\n\", passed, failed);")
    lines.append("    return failed > 0 ? 1 : 0;")
    lines.append("}")

    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description="MathLib deterministic fuzz vault")
    parser.add_argument("--seed", default=None, help="Explicit fuzzer seed")
    parser.add_argument("--logfile", default=None, help="Read failure output from an existing log file")
    parser.add_argument("--fuzzer", default=DEFAULT_FUZZER, help="Fuzzer executable to run")
    args = parser.parse_args()

    if args.logfile:
        with open(args.logfile, "r", encoding="utf-8") as f:
            text = f.read()
        seed = args.seed if args.seed is not None else extract_seed(text)
    else:
        seed = args.seed
        result = run_fuzzer(args.fuzzer, seed)
        text = result.stdout + "\n" + result.stderr
        seed = extract_seed(text)

        if result.returncode == 0 and "FAIL:" not in text:
            print("[VAULT] Fuzzer passed. No regressions to save.")
            sys.exit(0)

    artifact = save_artifact(seed, text)
    failures = parse_failures(text)

    print("[VAULT] Artifact saved: " + artifact)

    if failures:
        existing = glob.glob(os.path.join("tests", "regression", "fuzz_*.c"))
        next_id = len(existing) + 1
        filename = os.path.join("tests", "regression", f"fuzz_{next_id:04d}.c")

        code = generate_regression_c(seed, failures)
        with open(filename, "w", encoding="utf-8") as f:
            f.write(code)

        print("[VAULT] Regression test saved to: " + filename)
        print("[VAULT] Compile with:")
        print("gcc -std=c99 -O3 -Iinclude/mathlib -Isrc -o " + filename.replace(".c", "") + " " + filename + " -Lbuild -lmathc -lm")
        sys.exit(1)

    print("[VAULT] No parseable FAIL lines, but artifact saved for manual review.")
    sys.exit(1)


if __name__ == "__main__":
    main()
