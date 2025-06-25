# Critical Debugging Clues for Xic Model

## Priority 1: Parameter Mismatches

### Coupling Factors to Check
These coupling factors appear in the Julia implementation - verify they match Python:

| Wave | Julia Coupling Factor | Expected Value |
|------|---------------------|----------------|
| K2(1430)1 | 0.4472135954999579 | 1/√5 ≈ 0.447 |
| K(700)1 | 1.0 | 1.0 |
| L(1405)1 | 0.7071067811865475 | 1/√2 ≈ 0.707 |
| D(1232)1 | 0.5 | 0.5 |

**Action**: Check if Python uses the same coupling factors for these waves.

### Phase Space Point
Crosscheck uses this specific phase space point:
- `m2kpi`: 0.7980703453578917
- `m2pk`: 3.6486261122281745
- `m2ppi`: 2.7875826337931926

**Action**: Ensure Python uses exactly these values.

## Priority 2: Off-diagonal Element Issues

### Pattern: Diagonal Elements Match, Off-diagonal Don't
Most waves show this pattern:
- A++ and A-- elements match between Julia and Python
- A+- and A-+ elements differ significantly

**Examples**:
- **K2(1430)1**: Diagonal (0,0) and (1,1) match, off-diagonal (0,1) differs by factor ~50
- **L(1405)1**: 3/4 elements match exactly, only A+- differs

**Hypothesis**: Issue with helicity mixing or angular momentum coupling in off-diagonal elements.

**Action**: Check how Python constructs the A+- and A-+ matrix elements vs Julia.

## Priority 3: Specific Wave Analysis

### K2(1430)1 - Most Problematic
```
Julia: [0.0+0.0im  0.0933271+0.00052055im]
Python: [0.0+0.0im  -0.401402-5.37907im]
```
- **Ratio**: -0.00138 + 0.0172i (completely different)
- **Coupling factor**: 0.447 (1/√5)

**Action**: Debug this specific wave in Python - check coupling application and matrix element construction.

### K(700)1 - Red Flag Case
```
Julia: [-0.0+0.0im    -0.0+0.0im]
Python: [0.0+0.0im      -0.494729-5.37959im]
```
- **Red flag**: One zero, one nonzero in (0,0) element
- **Coupling factor**: 1.0

**Action**: Check why Python produces nonzero value where Julia produces zero.

### L(1405)1 - Good Baseline
```
Julia: [-0.180229+0.0440099im  0.504805-0.123268im]
Python: [-0.180229+0.0440099im  0.0100762-5.50286im]
```
- **3/4 elements match exactly**
- **Only A+- differs**: 0.504805-0.123268i vs 0.0100762-5.50286i

**Action**: Use this as a baseline - the matching elements show the calculation is mostly correct.

## Priority 4: Systematic Issues

### No Helicity Flip Problems
- All waves show "original with 3/4 matches" or similar
- No evidence of systematic matrix ordering issues
- **Conclusion**: Matrix element ordering is correct

### No Selection Rule Issues
- No systematic red flags (one-zero/one-nonzero patterns)
- When elements are zero in one calculation, they're typically zero in both
- **Conclusion**: Selection rules are consistent

### Magnitude vs Phase Issues
- Some mismatches: similar magnitude, different phase
- Others: completely different magnitude
- **Conclusion**: Both phase convention and normalization issues possible

## Debugging Strategy

### Step 1: Parameter Verification
1. Print all resonance parameters in Python (mass, width, couplings)
2. Compare with Julia model definitions
3. Check if crosscheck data uses different parameters

### Step 2: Amplitude Construction
1. For K2(1430)1, print intermediate calculation steps:
   - Breit-Wigner value
   - Angular factors
   - Coupling application
   - Final matrix element

### Step 3: Matrix Element Analysis
1. Focus on A+- and A-+ elements (most problematic)
2. Check helicity mixing calculations
3. Verify angular momentum coupling

### Step 4: Phase Conventions
1. Check global phase factors
2. Verify phase conventions in amplitude construction
3. Look for any normalization factors

## Key Questions for Python Code

1. **Are the coupling factors the same?** (Check the values above)
2. **How are A+- and A-+ matrix elements constructed?**
3. **Are there any global phase factors applied?**
4. **What normalization conventions are used?**
5. **Are the resonance parameters exactly the same?**

## Expected Outcome

Based on the 73.84% match rate and 3/4 elements matching pattern, the issue is likely:
- **Parameter mismatches** (most likely)
- **Phase convention differences** (possible)
- **Coupling normalization** (possible)

The fundamental amplitude construction appears correct, suggesting the issue is in the details rather than the core algorithm. 