using Test
using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
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



@testset "Exact values of amplitude and intensity" begin
    @test _I ≈ 9345.853380852352
    @test _A ≈ -45.1323269502508 + 54.85942516648639im
    # 
    @test model.chains[1].Xlineshape(σs0.σ2) ≈
          model.chains[2].Xlineshape(σs0.σ2) ≈ -0.5636481410171861 + 0.13763637759224928im
    # 
    @test model.chains[21].Xlineshape(σs0.σ1) ≈
          model.chains[22].Xlineshape(σs0.σ1) ≈
          model.chains[23].Xlineshape(σs0.σ1) ≈
          model.chains[24].Xlineshape(σs0.σ1) ≈ 2.1687201455088894 + 23.58225917009096im
end

using DelimitedFiles
using YAML
using OrderedCollections

data = readdlm(joinpath(@__DIR__, "..", "data", "raw-xic-values.dlm"))

function consistent_label(name, coupling)
    if startswith(name, "K*")
        coupling == "1/2,0" && return 1
        coupling == "1/2,-1" && return 2
        coupling == "-1/2,1" && return 3
        coupling == "-1/2,0" && return 4
    end
    if startswith(name, "K0")
        coupling == "1/2,0" && return 2
        coupling == "-1/2,0" && return 1
    end
    if startswith(name, "Λ") || startswith(name, "∆")
        coupling == "1/2,0" && return 1
        coupling == "-1/2,0" && return 2
    end
    return -1
end

processed_data = map(eachrow(data)) do row
    type, name, couplings = row[1], row[2], row[3]
    if last(type) == 'H'
        index = consistent_label(name, couplings)
        _reim = (row[1] == "ReH") ? "r" : ((row[1] == "ImH") ? "i" : "")
        _name = replace(name, "Λ" => "L", "∆" => "D", "*" => "", "+" => "", ")0" => ")", "0(" => "(")
        return "A$(_reim)$(_name)$(index)" => "$(row[4]) ± $(row[5]) ± $(row[6]) ± $(row[7])"
    end
    if first(type) == 'P'
        _name = replace(name, "γ" => "gamma", "*" => "", "+" => "", ")0" => ")", "0(" => "(")
        return "$(_name)" => "$(row[4]) ± $(row[5]) ± $(row[6]) ± $(row[7])"
    end
    error("Unknown name format: $type")
end

YAML.write_file(joinpath(@__DIR__, "..", "data", "xic2pKpi.yaml"), Dict("parameters" => LittleDict(processed_data...)))
# processed_data