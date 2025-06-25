cd(joinpath(@__DIR__, ".."))
using Pkg
Pkg.activate(".")
Pkg.instantiate()
#

using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
import Plots.PlotMeasures.mm
using DelimitedFiles
using ProgressMeter
using Measurements
using DataFrames
using Statistics
using QuadGK
using Plots
using Test
using YAML

theme(:boxed)


# -------------------------------------------------------------
# Load model and particle definitions from YAML files
# -------------------------------------------------------------
particledict = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-particle-definitions.yaml"))
modelparameters = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-model-definitions.yaml"))

default_model = begin
    defaultparameters = modelparameters["Default amplitude model"]
    parced_settings = parse_model_dictionaries(defaultparameters; particledict)
    Lc2ppiKModel(; parced_settings...)
end


alt_model = begin
    parameters = modelparameters["Alternative amplitude model"]
    parced_settings = parse_model_dictionaries(parameters; particledict)
    Lc2ppiKModel(; parced_settings...)
end


using Setfield
alt_model.names[7]
alt_model.names[8]

alt_model1 = @set alt_model.couplings[7] = alt_model.couplings[7] * 0.8
alt_model2 = @set alt_model1.couplings[8] = alt_model1.couplings[8] * 0.8

let
    x_values = range(1.62, 1.72, 51) .^ 2
    k = 2
    plot(xlab = "m²(Kp) [GeV²]")
    plot(x_values; xlab = "m²(Kp) [GeV²]") do σ
        I = Base.Fix1(unpolarized_intensity, alt_model2)
        integrand = projection_integrand(I, masses(alt_model2), σ; k)
        quadgk(integrand, 0, 1)[1]
    end
    _model_1670 = alt_model2[alt_model2.names.=="L(1670)"]
    plot!(x_values; xlab = "m²(Kp) [GeV²]", fill = 0, fillalpha = 0.5) do σ
        I = Base.Fix1(unpolarized_intensity, _model_1670)
        integrand = projection_integrand(I, masses(alt_model2), σ; k)
        quadgk(integrand, 0, 1)[1]
    end
end


# -------------------------------------------------------------
# Plot projections for L1670 components
# -------------------------------------------------------------
let
    x_values = range(1.62, 1.72, 51) .^ 2
    k = 2
    plot(xlab = "m²(Kp) [GeV²]")
    plot(x_values; xlab = "m²(Kp) [GeV²]") do σ
        I = Base.Fix1(unpolarized_intensity, default_model)
        integrand = projection_integrand(I, masses(default_model), σ; k)
        quadgk(integrand, 0, 1)[1]
    end
    _model_1670 = default_model[default_model.names.=="L(1670)"]
    plot!(x_values; xlab = "m²(Kp) [GeV²]", fill = 0, fillalpha = 0.5) do σ
        I = Base.Fix1(unpolarized_intensity, _model_1670)
        integrand = projection_integrand(I, masses(default_model), σ; k)
        quadgk(integrand, 0, 1)[1]
    end
end
