# YAML Structure Analysis: XiC → pKπ Amplitude Model

## Overview

The XiC → pKπ amplitude model workspace uses **three main YAML files** to define the complete physical model. These files provide a human-readable, structured approach to specifying particle properties, model configurations, and fit parameters.

## YAML File Structure

### 1. **xic-particle-definitions.yaml**
**Purpose**: Defines physical properties of all particles involved in the decay

**Structure**:
```yaml
# Initial and final state particles
Lambda_c+:
  latex: \Lambda_c^+
  jp: 1/2^+           # Spin-parity quantum numbers
  mass: 2467.94       # Mass in MeV
  width: 3.25e-9      # Width in MeV

# Resonances
L(1405):
  latex: \Lambda(1405)
  jp: 1/2^-
  mass: 1405.1
  width: 50.5
  lineshape: Flatte1405  # Specific lineshape type
```

**Key Features**:
- **Particle definitions**: Masses, widths, quantum numbers (JP)
- **LaTeX names**: For presentation and plotting
- **Lineshape hints**: Some particles specify their preferred lineshape

### 2. **xic-model-definitions.yaml**
**Purpose**: Specifies model structure, lineshape assignments, and coupling parameters

**Structure**:
```yaml
Default amplitude model:
  lineshapes:
    L(1405): Flatte1405
    L(1520): BreitWignerMinL
    K(892): BreitWignerMinL
    # ... more resonances
    
  parameters:
    # Complex amplitude coefficients with uncertainties
    ArL(1405)1: "-0.736593 ± 0.074727"  # Real part
    AiL(1405)1: "0.738192 ± 0.040446"   # Imaginary part
    # ... more parameters
```

**Key Features**:
- **Lineshape mapping**: Associates each resonance with its lineshape function
- **Parameter values**: Central values with uncertainties for all floating parameters
- **Helicity structure**: Parameters indexed by helicity (1, 2) for different spin projections

### 3. **xic2pKpi.yaml** (Experimental fit results)
**Purpose**: Contains experimental fit results with systematic uncertainties

**Structure**:
```yaml
parameters:
  ArL(1405)1: "-0.88 ± 0.22 ± 0.69 ± 0.17"  # stat ± syst1 ± syst2 ± syst3
  AiL(1405)1: "0.76 ± 0.13 ± 0.6 ± 0.44"
  # ... more parameters with full uncertainty breakdown
```

## YAML Loading Workflow

### 1. **File Loading Pattern**
```julia
# Load particle and model definitions
particledict = YAML.load_file("data/xic-particle-definitions.yaml")
modelparameters = YAML.load_file("data/xic-model-definitions.yaml")
defaultparameters = modelparameters["Default amplitude model"]

# Parse into Julia data structures
(; chains, couplings, isobarnames) = parse_model_dictionaries(defaultparameters; particledict)

# Create the amplitude model
model = Lc2ppiKModel(; chains, couplings, isobarnames)
```

### 2. **Key Processing Steps**
1. **Mass extraction**: Converts MeV → GeV and creates `ThreeBodyMasses`
2. **Lineshape assignment**: Maps resonance names to lineshape functions
3. **Parameter parsing**: Converts string format "value ± uncertainty" to `MeasuredParameter`
4. **Chain construction**: Builds decay chain objects with proper quantum numbers
5. **Model assembly**: Creates final `Lc2ppiKModel` ready for calculations

### 3. **Parameter Name Convention**
- **Format**: `A{r/i}{resonance_name}{helicity_index}`
- **Examples**:
  - `ArL(1405)1`: Real part of L(1405) coupling, helicity +1/2
  - `AiK(892)2`: Imaginary part of K(892) coupling, helicity -1/2
  - `ArD(1232)1`: Real part of D(1232) coupling, helicity +1/2

## YAML vs JSON Comparison

### **YAML Advantages**:
1. **Human readability**: Easy to edit and understand
2. **Comments**: Can include documentation inline
3. **Natural structure**: Mirrors physics organization
4. **Version control friendly**: Clear diffs when parameters change
5. **Uncertainty notation**: Supports "value ± error" format naturally

### **JSON Advantages**:
1. **Tool integration**: Works with many external tools
2. **Standardized format**: Well-defined schema for amplitude models
3. **Validation**: Can be validated against JSON schemas
4. **Cross-platform**: Universal support across programming languages

### **Conversion Process**:
```julia
# YAML → Julia model → JSON serialization
YAML files → parse_model_dictionaries() → Lc2ppiKModel() → serializeToDict() → JSON
```

## Model Structure Hierarchy

```
XiC → pKπ Model
├── Particle Properties (xic-particle-definitions.yaml)
│   ├── Initial state: Lambda_c+
│   ├── Final states: p, pi+, K-
│   └── Resonances: L(1405), L(1520), K(892), D(1232), etc.
│
├── Model Configuration (xic-model-definitions.yaml)
│   ├── Lineshape assignments
│   │   ├── BreitWignerMinL (standard)
│   │   ├── Flatte1405 (coupled-channel)
│   │   ├── L1670Flatte (special case)
│   │   └── BuggBreitWignerMinL (scalar mesons)
│   └── Amplitude parameters
│       ├── Lambda resonances: L(1405), L(1520), L(1600), ...
│       ├── Delta resonances: D(1232), D(1600), D(1620), ...
│       └── Kaon resonances: K(700), K(892), K(1430), ...
│
└── Experimental Data (xic2pKpi.yaml)
    └── Fit results with full systematic uncertainties
```

## Technical Implementation Details

### **Mass Convention**:
- YAML files: Masses in MeV
- Julia model: Masses in GeV (automatically converted)
- Numbering: 0: XiC, 1: p, 2: π, 3: K

### **Quantum Numbers**:
- Initial state: XiC (JP = 1/2+)
- Final states: p (1/2+), π (0-), K (0-)
- Resonances: Various JP assignments specified in YAML

### **Lineshape Types**:
- `BreitWignerMinL`: Standard relativistic Breit-Wigner
- `Flatte1405`: Coupled-channel for Λ(1405)
- `L1670Flatte`: Special treatment for Λ(1670)
- `BuggBreitWignerMinL`: Bugg parameterization for scalars

### **Physical Validation**:
- Kinematic constraints ensure physical phase space
- Spin-parity conservation checked
- Unitarity respected in lineshape definitions

## Usage Examples

### **Basic Model Loading**:
```julia
using Lc2ppiKSemileptonicModelLHCb
using YAML

# Load and parse YAML model
particledict = YAML.load_file("data/xic-particle-definitions.yaml")
modelparameters = YAML.load_file("data/xic-model-definitions.yaml")
defaultparameters = modelparameters["Default amplitude model"]

# Create model
(; chains, couplings, isobarnames) = parse_model_dictionaries(defaultparameters; particledict)
model = Lc2ppiKModel(; chains, couplings, isobarnames)
```

### **Parameter Access**:
```julia
# Access masses
masses = masses(model)  # Gets ThreeBodyMasses object

# Access specific resonance parameters
# (Parameters are embedded in the model structure)
```

### **Amplitude Calculation**:
```julia
# Define phase space point
σs = (σ1=3.0, σ2=2.0, σ3=1.5)  # Invariant masses squared

# Calculate amplitude
amp = amplitude(model, σs, [λ1, λ2, λ3, λ0])  # Helicity configuration
```

## Conclusion

The YAML-based approach provides a **clean, human-readable way to specify complex amplitude models**. The three-file structure (particles, model, fit results) separates concerns effectively:

1. **Physical properties** remain stable across analyses
2. **Model structure** can be modified for systematic studies  
3. **Fit results** are kept separate for different datasets

This design enables easy model development, systematic studies, and collaboration between physicists who need to understand and modify the model parameters.
