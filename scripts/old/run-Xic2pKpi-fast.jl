# -------------------------------------------------------------
# FAST TEST VERSION - Reduced computational load for quick testing
# -------------------------------------------------------------
# Imports and dependencies
# -------------------------------------------------------------
using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using DelimitedFiles
using ProgressMeter
using Statistics
using Test
using JSON

# -------------------------------------------------------------
# Load model from JSON file (faster and avoids YAML parsing issues)
# -------------------------------------------------------------
model_file = joinpath(@__DIR__, "..", "data", "xic2pKpi-model_new.json")
if !isfile(model_file)
    error("JSON model file not found: $model_file. Please run xic_yaml_to_json_new.jl first.")
end

# Read the JSON model directly
model_json = read(model_file, String)
parsed_model = JSON.parse(model_json)

# Extract chains, couplings, and isobarnames from JSON
chains = parsed_model["chains"]
couplings = parsed_model["couplings"] 
isobarnames = parsed_model["isobarnames"]

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

println("✓ Basic model evaluation tests passed")
println("Model has $(length(model.chains)) decay chains")
println("Intensity at test point: $_I")
println("Amplitude at test point: $_A")

# -------------------------------------------------------------
# FAST VERSION: Reduced computational load
# -------------------------------------------------------------
using Plots
using QuadGK
using DataFrames
using Measurements
import Plots.PlotMeasures.mm

theme(:boxed)

# Quick Dalitz plot (reduced resolution)
println("\n--- Creating Dalitz plot (reduced resolution) ---")
@time plot(
    masses(model), σs -> unpolarized_intensity(model, σs);
    iσx = 2, iσy = 1, title = "Dalitz plot (fast)",
    xlab = "m²(pK⁻) [GeV²]", ylab = "m²(K⁻π⁺) [GeV²]",
    resolution = (200, 200))  # Reduced resolution for speed

# -------------------------------------------------------------
# 1D projections (reduced bins and faster computation)
# -------------------------------------------------------------
labels = [(2, "pK⁻"), (1, "K⁻π⁺"), (3, "pπ⁺")]
println("\n--- Creating 1D projections (reduced bins) ---")
let
    plot(layout = grid(1, 3), size = (1200, 400), bottom_margin = 5mm)
    n_bins = 50  # Reduced from 100 for speed
    @showprogress map(enumerate(labels)) do (sp, (k, part_lab))
        x_ranges = lims(masses(model); k)
        plot!(range(x_ranges..., n_bins); sp, xlab = "m²($part_lab) [GeV²]") do σ
            I = Base.Fix1(unpolarized_intensity, model)
            integrand = projection_integrand(I, masses(model), σ; k)
            quadgk(integrand, 0, 1)[1]
        end
    end
    plot!()
end

# -------------------------------------------------------------
# Generate random phase space points (REDUCED SAMPLE SIZE)
# -------------------------------------------------------------
println("\n--- Generating random samples (10k instead of 100k) ---")
data = map(eachrow(rand(10_000, 2))) do y  # Reduced from 100_000 to 10_000
    y2σs(y, masses(model))
end |> v -> DataFrame(σs = v)
subset!(data, :σs => ByRow(σs -> isphysical(σs, masses(model))))

println("Generated $(nrow(data)) physical phase space points")

@time transform!(data,
    :σs => ByRow(σs -> unpolarized_intensity(model, σs)) => :weights,
    :σs => ByRow(identity) => AsTable)

# -------------------------------------------------------------
# Plot histograms (reduced bins)
# -------------------------------------------------------------
println("\n--- Creating histograms (reduced bins) ---")
let
    plot(layout = grid(1, 3), size = (1200, 400),
        ylabel = "Candidate", bottom_margin = 5mm)
    map(enumerate(labels)) do (sp, (k, part_lab))
        σ = Symbol("σ$k")
        stephist!(data[!, σ]; data.weights, sp,
            xlabel = "m²($part_lab) [GeV]",
            bins = 50, normalize = true)  # Reduced from 100 bins
    end
    plot!()
end

# -------------------------------------------------------------
# Calculate weights for each amplitude chain (component)
# -------------------------------------------------------------
println("\n--- Computing component weights ---")
@showprogress for name in unique(model.names)
    _model = model[model.names.==name]
    _weight = Symbol("weights_$name")
    transform!(data,
        :σs => ByRow(σs -> unpolarized_intensity(_model, σs)) => _weight)
end

# -------------------------------------------------------------
# Reduce the weights to sum and compute fit fractions
# -------------------------------------------------------------
println("\n--- Computing fit fractions ---")
stacked_weights = let
    _table = leftjoin(
        stack(
            combine(data,
                Not(:σs, :σ1, :σ2, :σ3) .=> sum; renamecols = false),
            value_name = :mean,
        ),
        stack(
            combine(data,
                Not(:σs, :σ1, :σ2, :σ3) .=> std; renamecols = false),
            value_name = :std,
        ); on = :variable)
    # create value as  m ± σ / sqrt(N)
    transform!(_table, [:mean, :std] => ByRow((m, σ) -> m ± σ / sqrt(size(data, 1))) => :value)
    # normalize the weights
    select!(_table, :variable, :value => (x -> x ./ _table[1, :value] * 100) =>
        :fraction)
    # sort by fraction
    sort!(_table, [:fraction]; rev = true)
    subset!(_table, :variable => ByRow(x -> occursin("weights_", x)))
    _table.variable .= map(x -> x[9:end], _table.variable)
    _table
end

# -------------------------------------------------------------
# Load reference fit fractions and compare
# -------------------------------------------------------------
println("\n--- Comparing with reference fit fractions ---")
fractions_ref = let
    _data = readdlm(joinpath(@__DIR__, "..", "data", "fit-fractions-ref.txt"))[:, 1:3]
    DataFrame(
        variable = _data[:, 1] |> collect,
        ref_fraction = _data[:, 2] .± _data[:, 3])
end
fit_fractions = leftjoin(stacked_weights, fractions_ref; on = :variable)

println("Fit fractions comparison:")
print(DataFrames.pretty_table(fit_fractions; compact_printing = false))

# Normalize fractions
let
    _, i = findmax(fit_fractions.fraction)
    fit_fractions[:, 2] ./= fit_fractions[i, 2] / 100
    fit_fractions[:, 3] ./= fit_fractions[i, 3] / 100
    fit_fractions
end

# Calculate ratios
ratios = getproperty.(fit_fractions.fraction, :val) ./ fit_fractions.ref_fraction
println("\nRatios (computed/reference):")
for (i, row) in enumerate(eachrow(fit_fractions))
    println("$(row.variable): $(round(ratios[i], digits=3))")
end

# -------------------------------------------------------------
# Test individual resonances
# -------------------------------------------------------------
println("\n--- Testing individual resonances ---")

# Filter individual resonances by name
lambda_1520_model = model[model.names .== "L(1520)"]
k892_model = model[model.names .== "K(892)"]

# Test with Mandelstam variables
σs0 = Invariants(model.chains[1].tbs.ms; σ1 = 0.7980703453578917, σ2 = 3.2)

# Calculate amplitude for individual resonances
amplitude_L1520 = amplitude(lambda_1520_model, σs0, [1, 0, 0, 1])
intensity_L1520 = unpolarized_intensity(lambda_1520_model, σs0)

amplitude_K892 = amplitude(k892_model, σs0, [1, 0, 0, 1])
intensity_K892 = unpolarized_intensity(k892_model, σs0)

println("L(1520) - Amplitude: $amplitude_L1520, Intensity: $intensity_L1520")
println("K(892) - Amplitude: $amplitude_K892, Intensity: $intensity_K892")

# -------------------------------------------------------------
# Summary
# -------------------------------------------------------------
println("\n" * "="^60)
println("FAST TEST COMPLETED SUCCESSFULLY")
println("="^60)
println("Model chains: $(length(model.chains))")
println("Data points: $(nrow(data))")
println("Total intensity: $(sum(data.weights))")
println("Number of resonances: $(length(unique(model.names)))")

# Test that the model evaluates correctly
@test _I > 0
@test abs(_A) > 0
@test nrow(data) > 1000  # Should have enough data points

println("✓ All fast tests passed!")
