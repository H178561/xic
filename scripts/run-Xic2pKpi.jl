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

# -------------------------------------------------------------
# Plotting Dalitz plot and projections
# -------------------------------------------------------------
using Plots
using QuadGK
using DataFrames
using Measurements
import Plots.PlotMeasures.mm

theme(:boxed)

# Dalitz plot: 2D distribution of invariant masses
@time plot(
    masses(model), σs -> unpolarized_intensity(model, σs);
    iσx = 2, iσy = 1, title = "Dalitz plot",
    xlab = "m²(pK⁻) [GeV²]", ylab = "m²(K⁻π⁺) [GeV²]")

# -------------------------------------------------------------
# 1D projections of invariant mass distributions
# -------------------------------------------------------------
labels = [(2, "pK⁻"), (1, "K⁻π⁺"), (3, "pπ⁺")]
let
    plot(layout = grid(1, 3), size = (1200, 400), bottom_margin = 5mm)
    n_bins = 100
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
# Generate random phase space points and compute weights
# -------------------------------------------------------------
# Histogram of the invariant mass distributions
# -------------------------------------------------------------
data = map(eachrow(rand(100_000, 2))) do y
    y2σs(y, masses(model))
end |> v -> DataFrame(σs = v)
subset!(data, :σs => ByRow(σs -> isphysical(σs, masses(model))))
@time transform!(data,
    :σs => ByRow(σs -> unpolarized_intensity(model, σs)) => :weights,
    :σs => ByRow(identity) => AsTable)

# -------------------------------------------------------------
# Plot histograms of invariant mass distributions
# -------------------------------------------------------------
let
    plot(layout = grid(1, 3), size = (1200, 400),
        ylabel = "Candidate", bottom_margin = 5mm)
    map(enumerate(labels)) do (sp, (k, part_lab))
        σ = Symbol("σ$k")
        stephist!(data[!, σ]; data.weights, sp,
            xlabel = "m²($part_lab) [GeV]",
            bins = 100, normalize = true)
    end
    plot!()
end

# -------------------------------------------------------------
# Calculate weights for each amplitude chain (component)
# -------------------------------------------------------------
@showprogress for name in unique(model.names)
    _model = model[model.names.==name]
    _weight = Symbol("weights_$name")
    transform!(data,
        :σs => ByRow(σs -> unpolarized_intensity(_model, σs)) => _weight)
end

# -------------------------------------------------------------
# Reduce the weights to sum and compute fit fractions
# -------------------------------------------------------------
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
# Load reference fit fractions from file and compare
# -------------------------------------------------------------
fractions_ref = let
    _data = readdlm(joinpath(@__DIR__, "..", "data", "fit-fractions-ref.txt"))[:, 1:3]
    DataFrame(
        variable = _data[:, 1] |> collect,
        ref_fraction = _data[:, 2] .± _data[:, 3])
end
fit_fractions = leftjoin(stacked_weights, fractions_ref; on = :variable)

print(DataFrames.pretty_table(fit_fractions; compact_printing = false))
let
    _, i = findmax(fit_fractions.fraction)
    fit_fractions[:, 2] ./= fit_fractions[i, 2] / 100
    fit_fractions[:, 3] ./= fit_fractions[i, 3] / 100
    fit_fractions
end

# -------------------------------------------------------------
# Plot projections for most significant amplitude components
# -------------------------------------------------------------
let
    most_significant = stacked_weights[2:end, :].variable
    plot(layout = grid(1, 3), size = (1200, 400), bottom_margin = 5mm, yaxis = nothing)
    map(enumerate(labels)) do (sp, (k, part_lab))
        σ = Symbol("σ$k")
        stephist!(data[!, σ]; data.weights, sp,
            xlabel = "m²($part_lab) [GeV]",
            bins = 100)
        map(most_significant) do branch
            @show branch
            data_branch_name = "weights_" * branch
            stephist!(data[!, σ]; weights = getproperty(data, data_branch_name),
                sp, lab = branch,
                bins = 100)
        end
    end
    plot!()
end
# -------------------------------------------------------------
# Save the final projections plot
# -------------------------------------------------------------
savefig(joinpath(@__DIR__, "..", "plots", "xic2pKpi-projections.png"))

