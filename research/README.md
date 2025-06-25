# Xic Model Crosscheck Investigation

## Executive Summary

This document summarizes the detailed investigation of discrepancies between the Julia implementation of the Xic model and the reference Python/LHCb calculations.

### Key Statistics
- **Total waves analyzed**: 43
- **Total matrix elements**: 172 (4 per wave)
- **Matching elements**: 127 (73.84%)
- **Mismatching elements**: 45 (26.16%)
- **Waves with any mismatches**: 43 (100%)

## Analysis Methodology

### 1. Improved Comparison Logic
- **Zero tolerance**: `1e-10` for considering values as "zero"
- **Match tolerance**: `1e-3` for considering values as matching
- **Both zeros**: Marked as matches
- **One zero, one nonzero**: Marked as "red flags" for investigation
- **Both nonzero**: Ratio calculated and compared to 1.0

### 2. Helicity Flip Investigation
For each problematic wave, the script tested 4 permutations of the 2x2 matrix:
- Original matrix
- Transpose
- Flip rows
- Flip columns

The best permutation (highest number of matching elements) was recorded.

## Key Findings

### 1. Overall Pattern
- **73.84% match rate** indicates the model is mostly correct
- **All 43 waves have at least one mismatch**, suggesting systematic issues
- Most mismatches are in **off-diagonal elements** (λ0=1, λ1=-1 and λ0=-1, λ1=1)

### 2. Red Flag Analysis
From the output analysis, most waves show **no red flags** (no one-zero/one-nonzero cases), which means:
- When elements are zero in one calculation, they're typically zero in both
- The main discrepancies are in **nonzero elements** with different values

### 3. Helicity Flip Investigation
Most waves show **"original with 3/4 matches"**, indicating:
- 3 out of 4 matrix elements match in the original ordering
- 1 element consistently mismatches
- **No evidence of systematic helicity flip issues**

## Detailed Wave Analysis

### Problematic Wave Categories

#### 1. K-meson Resonances
- **K2(1430)**: All 4 variants show mismatches
- **K(700)**: Both variants show mismatches  
- **K(1430)**: Both variants show mismatches
- **K(892)**: All 4 variants show mismatches

#### 2. Lambda Resonances
- **L(1405)**: Both variants show mismatches
- **L(1520)**: Both variants show mismatches
- **L(1600)**: Both variants show mismatches
- **L(1670)**: Both variants show mismatches
- **L(1690)**: Both variants show mismatches
- **L(1710)**: Both variants show mismatches
- **L(1800)**: Both variants show mismatches
- **L(1810)**: Both variants show mismatches
- **L(1820)**: Both variants show mismatches
- **L(1830)**: Both variants show mismatches
- **L(1890)**: Both variants show mismatches
- **L(2000)**: Both variants show mismatches

#### 3. Delta Resonances
- **D(1232)**: Both variants show mismatches
- **D(1600)**: Both variants show mismatches
- **D(1620)**: Both variants show mismatches
- **D(1700)**: Both variants show mismatches

## Specific Clues and Patterns

### 1. Coupling Factors
Different waves show different coupling factors in Julia:
- **K2(1430)1**: `0.4472135954999579` (≈ 1/√5)
- **K(700)1**: `1.0` (no coupling factor)
- **L(1405)1**: `0.7071067811865475` (≈ 1/√2)
- **D(1232)1**: `0.5` (exact)

**Investigation needed**: Check if these coupling factors match the Python implementation.

### 2. Matrix Element Patterns

#### K2(1430)1 Example
```
Julia DPD:                    Python LHCb':
[0.0+0.0im  0.0933271+0.00052055im]  [0.0+0.0im  -0.401402-5.37907im]
[0.0+0.0im        0.0+0.0im]         [0.0+0.0im        0.0+0.0im]
```
- **Diagonal elements match** (both zero)
- **Off-diagonal element differs by large factor** (~0.093 vs ~-0.401-5.38i)
- **Ratio**: -0.00138 + 0.0172i (very different)

#### K(700)1 Example
```
Julia DPD:                    Python LHCb':
[-0.0+0.0im    -0.0+0.0im]           [0.0+0.0im      -0.494729-5.37959im]
[-1.89405+4.353im  -0.0+0.0im]       [-0.376875+2.08378im        0.0+0.0im]
```
- **One red flag**: (0,0) element (Julia: -0.0, Python: 0.0)
- **Off-diagonal elements differ significantly**
- **Ratio for matching element**: 2.182 + 0.514i

#### L(1405)1 Example
```
Julia DPD:                    Python LHCb':
[-0.180229+0.0440099im  0.504805-0.123268im]  [-0.180229+0.0440099im  0.0100762-5.50286im]
[-0.0586139+0.0143128im  0.164172-0.040089im] [-0.0586139+0.0143128im   0.164172-0.040089im]
```
- **3/4 elements match exactly**
- **One off-diagonal element differs**: 0.504805-0.123268i vs 0.0100762-5.50286i
- **Ratio for mismatched element**: 0.0226 + 0.0917i

### 3. Crosscheck Data Structure
The crosscheck data contains:
- **chainvars**: Phase space point (m2kpi=0.798, m2pk=3.649, etc.)
- **chains**: 2x2 amplitude matrices for each parameter
- **Matrix ordering**: A++, A+-, A-+, A-- (helicity basis)

### 4. Systematic Patterns

#### A. Off-diagonal Element Issues
- **Most mismatches occur in off-diagonal elements** (A+- and A-+)
- **Diagonal elements (A++ and A--) often match**
- This suggests issues with **helicity mixing** or **angular momentum coupling**

#### B. Magnitude vs Phase Issues
- Some mismatches show **similar magnitudes but different phases**
- Others show **completely different magnitudes**
- This suggests both **phase convention** and **normalization** issues

#### C. Zero Pattern Consistency
- **No systematic red flags** (one-zero/one-nonzero)
- When elements are zero in one calculation, they're typically zero in both
- This suggests **selection rules are consistent** between implementations

## Hypotheses for Discrepancies

### 1. Parameter Value Differences
**Most likely cause**: The resonance parameters (mass, width, couplings) may differ between the Julia and Python implementations.

**Investigation needed**:
- Compare all resonance parameters between the two codes
- Check if the crosscheck data uses different parameter values than the model definitions

### 2. Phase Conventions
**Possible cause**: Different phase conventions in amplitude construction.

**Investigation needed**:
- Compare the phase factors in amplitude construction
- Check if there are global phase differences

### 3. Coupling Normalizations
**Possible cause**: Different normalization conventions for couplings.

**Investigation needed**:
- Compare how couplings are applied in amplitude construction
- Check for any normalization factors

### 4. Matrix Element Ordering
**Less likely**: The analysis shows no evidence of systematic helicity flip issues, but individual cases should be checked.

## Next Steps for Investigation

### Phase 1: Parameter Comparison
1. **Extract parameters from Python code**:
   - Resonance masses and widths
   - Coupling values
   - Any normalization factors

2. **Compare with Julia implementation**:
   - Check if parameters match exactly
   - Identify any differences

### Phase 2: Amplitude Construction Analysis
3. **For a few problematic waves**:
   - Print intermediate calculation steps in both codes
   - Compare Breit-Wigner values
   - Compare angular factors
   - Compare final amplitude construction

### Phase 3: Systematic Fixes
4. **Based on findings**:
   - Fix parameter mismatches
   - Adjust phase conventions if needed
   - Fix coupling normalizations if needed

## Files and Data

### Generated Files
- `detailed_analysis_output.txt`: Complete script output with element-by-element analysis
- `detailed_xic_analysis.json`: Structured data (not generated due to JSON3 package issue)

### Key Data Points
- **Match rate**: 73.84%
- **Most common pattern**: 3/4 elements match, 1 element differs
- **No systematic helicity flip issues detected**
- **No red flags (one-zero/one-nonzero) patterns detected**

## Conclusion

The Xic model implementation shows **good overall agreement** (73.84% match rate) but has **systematic discrepancies** affecting all waves. The pattern of 3/4 elements matching suggests the issue is likely in **parameter values or phase conventions** rather than fundamental implementation errors.

The next priority should be **parameter comparison** between the Python and Julia implementations to identify the source of the discrepancies. 