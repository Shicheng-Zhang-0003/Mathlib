# Math Diary: Exponentials & Logarithms (Horner's & Cody-Waite)

> **v11S Closure Status:** Exp/log use Cody-Waite reduction and Horner evaluation. Pow, hyperbolics, and special-case behavior were hardened during v11S closure.
> **Orchestrator's Mission:** You cannot ship what you cannot explain.
> Use this document to reverse-engineer the AI-generated code in `src/exp_log.c`.
> Consult the Red Team (Gemini/ChatGPT) or Green Team (Copilot) to fill in the blanks.

---

## 1. The Mathematical Identity
*What fundamental mathematical truth or definition does this module enforce?*
- **Core Equation:** [Fill in the mathematical definition, e.g., Taylor series, Euler's formula]
- **Algebraic Invariants:** [e.g., sin²(θ) + cos²(θ) = 1, det(AB) = det(A)det(B)]

## 2. The Algorithm (Execution Path)
*Trace the data flow from the raw `double` input to the final `double` output. What are the distinct phases?*
1. **Phase 1: Range Reduction / Pre-processing**
   - *Code location:* `[Function Name]`
   - *Why we do it:* [Explain why we can't just plug the raw input into the polynomial]
2. **Phase 2: The Core Approximation**
   - *Code location:* `[Function Name]`
   - *Why we do it:* [Explain the polynomial or iterative method used]
3. **Phase 3: Reconstruction / Post-processing**
   - *Code location:* `[Function Name]`
   - *Why we do it:* [Explain how we map the reduced result back to the original domain]

## 3. The "Dark Arts" (Hardware & IEEE-754 Sympathy)
*Why did the AI write it this way instead of the naive textbook way?*
- **Bitwise / Memory Tricks:** [e.g., Why use `memcpy` instead of `union`? Why use Horner's Method?]
- **Branchless Logic:** [Are there masks or bitwise ops used to avoid CPU pipeline stalls?]
- **Constants:** [Where do the magic numbers (e.g., `0x5fe6ec85...`) come from?]

## 4. Error Budget & Failure Modes
*Where does the physical limit of the silicon break the mathematical ideal?*
- **Maximum ULP Error:** [ ] ULPs (Units in the Last Place) vs MPFR ground-truth.
- **Catastrophic Cancellation:** [Identify inputs where subtracting two nearly equal numbers destroys precision]
- **Overflow/Underflow Boundaries:** [At what exact input value does the result become `Inf` or `0.0`?]
- **NaN Propagation:** [How does the code handle `NaN` and `Inf` inputs?]

## 5. Orchestrator's Reverse-Engineering Log
*Paste insights, proofs, and "Aha!" moments from your Red/Green team consultations here.*
- [Date] - Insight 1: ...
- [Date] - Insight 2: ...
