# -------------------------------------------------------------
# Imports and dependencies
# -------------------------------------------------------------
using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using DelimitedFiles
using ProgressMeter
using Statistics
using Test
using YAML

# -------------------------------------------------------------
# Load model and particle definitions from YAML files
# -------------------------------------------------------------
begin
    particledict = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-particle-definitions.yaml"))
    modelparameters = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-model-definitions.yaml"))
end

defaultparameters = modelparameters["Default amplitude model"]

# -------------------------------------------------------------
# Parse model dictionaries and convert to standard convention
# -------------------------------------------------------------
(; chains, couplings, isobarnames) =
    parse_model_dictionaries(defaultparameters; particledict)

# -------------------------------------------------------------
# Set up the amplitude model with particle numbering
# 0: Xic, 1:p, 2:pi, 3:K
# -------------------------------------------------------------
model = Lc2ppiKModel(; chains, couplings, isobarnames)

# -------------------------------------------------------------
# Generate a random point in the phase space (Dalitz variables)
# -------------------------------------------------------------
σs0 = Invariants(model.chains[1].tbs.ms;
    σ1 = 0.7980703453578917,
    σ2 = 3.6486261122281745)

ms = model.chains[1].tbs.ms
σ1 = 1.4  # GeV²
σ2 = 3.2  # GeV²
σs = Invariants(ms; σ1, σ2)

# -------------------------------------------------------------
# Evaluate intensity and amplitude at the random point
# -------------------------------------------------------------
_I = unpolarized_intensity(model, σs)
_A = amplitude(model, σs)  # pars: model, mandelstam variables, helicity values

print(_I, "\n")
print(_A, "\n")
# -------------------------------------------------------------
# Basic tests to check evaluation
# -------------------------------------------------------------

