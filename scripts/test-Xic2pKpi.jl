using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using ProgressMeter
using Statistics
using Test
using YAML

begin
    particledict = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-particle-definitions.yaml"))
    modelparameters = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-model-definitions.yaml"))
end

defaultparameters = modelparameters["Default amplitude model"]

# get parameters from json files
# convert to the standard convention
(; chains, couplings, isobarnames) =
    parse_model_dictionaries(defaultparameters; particledict)


# set up the model with a numbering
# 0: Xic, 1:p, 2:pi, 3:K
model = Lc2ppiKModel(; chains, couplings, isobarnames)

# get a random point in the phase space
σs0 = Invariants(model.chains[1].tbs.ms;
    σ1=0.7980703453578917,
    σ2=3.6486261122281745)

# call intensity
_I = unpolarized_intensity(model, σs0)

# call the amplitude
_A = amplitude(model, σs0, [1, 0, 0, 1])  # pars: model, mandelstam variables, helicity values

@testset "Evaluation of the meeting" begin
    @test _I isa Real
    @test _A isa Complex
    @test _A == amplitude(model, σs0, ThreeBodySpins(1, 0, 0; two_h0=1))
end


using Plots
using QuadGK
using DataFrames
using Measurements
import Plots.PlotMeasures.mm


theme(:boxed)
@time plot(
    masses(model), σs -> unpolarized_intensity(model, σs);
    iσx=2, iσy=1, title="Dalitz plot",
    xlab="m²(pK⁻) [GeV²]", ylab="m²(K⁻π⁺) [GeV²]")

labels = [(2, "pK⁻"), (1, "K⁻π⁺"), (3, "pπ⁺")]
let
    plot(layout=grid(1, 3), size=(1200, 400), bottom_margin=5mm)
    n_bins = 100
    @showprogress map(enumerate(labels)) do (sp, (k, part_lab))
        x_ranges = lims(masses(model); k)
        plot!(range(x_ranges..., n_bins); sp, xlab="m²($part_lab) [GeV²]") do σ
            I = Base.Fix1(unpolarized_intensity, model)
            integrand = projection_integrand(I, masses(model), σ; k)
            quadgk(integrand, 0, 1)[1]
        end
    end
    plot!()
end

# histogram of the invariant mass distributions
data = map(eachrow(rand(1_000_000, 2))) do y
    y2σs(y, masses(model))
end |> v -> DataFrame(σs=v)
subset!(data, :σs => ByRow(σs -> isphysical(σs, masses(model))))
@time transform!(data,
    :σs => ByRow(σs -> unpolarized_intensity(model, σs)) => :weights,
    :σs => ByRow(identity) => AsTable)

let
    plot(layout=grid(1, 3), size=(1200, 400),
        ylabel="Candidate", bottom_margin=5mm)
    map(enumerate(labels)) do (sp, (k, part_lab))
        σ = Symbol("σ$k")
        stephist!(data[!, σ]; data.weights, sp,
            xlabel="m²($part_lab) [GeV]",
            bins=100, normalize=true)
    end
    plot!()
end

@showprogress for name in unique(model.names)
    _model = model[model.names.==name]
    _weight = Symbol("weights_$name")
    transform!(data,
        :σs => ByRow(σs -> unpolarized_intensity(_model, σs)) => _weight)
end


# reduce the weights to sum
stacked_weights = leftjoin(
    stack(
        combine(data,
            Not(:σs, :σ1, :σ2, :σ3) .=> sum; renamecols=false),
        value_name=:mean,
    ),
    stack(
        combine(data,
            Not(:σs, :σ1, :σ2, :σ3) .=> std; renamecols=false),
        value_name=:std,
    ); on=:variable)
# create value as  m ± σ / sqrt(N)
transform!(stacked_weights, [:mean, :std] => ByRow((m, σ) -> m ± σ / sqrt(size(data, 1))) => :value)
# normalize the weights
select!(stacked_weights, :variable, :value => (x -> x ./ stacked_weights[1, :value] * 100) =>
    :fraction)
# sort by fraction
sort!(stacked_weights, [:fraction]; rev=true)

print(DataFrames.pretty_table(stacked_weights; compact_printing=false))


let
    most_significant = stacked_weights[2:end, :].variable
    plot(layout=grid(1, 3), size=(1200, 400), bottom_margin=5mm, yaxis=nothing)
    map(enumerate(labels)) do (sp, (k, part_lab))
        σ = Symbol("σ$k")
        stephist!(data[!, σ]; data.weights, sp,
            xlabel="m²($part_lab) [GeV]",
            bins=100)
        map(most_significant) do branch
            lab = replace(branch, "weights_" => "")
            stephist!(data[!, σ]; weights=getproperty(data, branch),
                sp, lab,
                bins=100)
        end
    end
    plot!()
end
savefig(joinpath(@__DIR__, "..", "plots", "xic2pKpi-projections.png"))