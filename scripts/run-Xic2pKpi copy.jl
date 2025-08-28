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

# -------------------------------------------------------------
# Evaluate intensity and amplitude at the random point
# -------------------------------------------------------------
_I = unpolarized_intensity(model, σs0)
_A = amplitude(model, σs0, [1, 0, 0, 1])  # pars: model, mandelstam variables, helicity values

# -------------------------------------------------------------
# Basic tests to check evaluation
# -------------------------------------------------------------
@testset "Evaluation of the meeting" begin
    @test _I isa Real
    @test _A isa Complex
    @test _A == amplitude(model, σs0, ThreeBodySpins(1, 0, 0; two_h0 = 1))
end

println("Test point", σs0)
println("Intensity at test point: $_I\n")
print("Amplitude at test point: $_A\n")

# ----------
# -------------------------------------------------------------
# Save the final projections plot
# -------------------------------------------------------------
#savefig(joinpath(@__DIR__, "..", "plots", "xic2pKpi-projections.png"))


