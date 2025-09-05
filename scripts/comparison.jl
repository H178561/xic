import Pkg
Pkg.activate(dirname(@__DIR__))  # Geht ein Verzeichnis hoch zum Hauptprojekt
Pkg.instantiate()



using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
#using ThreeBodyDecaysIO.HadronicLineshapes
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.DataFrames
using ThreeBodyDecaysIO.JSON
using Measurements
using Statistics
using QuadGK
using Plots
using StatsBase
using FHist
using CSV
using YAML

#include("/home/hvdsmagt/ThreeBodyDecays.jl/IO.jl/HadronicLineshapes/HadronicLineshapes.jl/src/shapes.jl")


# Load the JSON model
input = open(
    #joinpath("/home/hvdsmagt/ThreeBodyDecays.jl/IO.jl/", "lc2ppik-lhcb-2683025.json"),
    joinpath(dirname(@__DIR__), "data", "xic2pKpi-model_test.json")
) do io
    JSON.parse(io)
end

workspace = Dict{String,Any}()

@unpack functions = input
for fn in functions
    @unpack name, type = fn
    instance_type = eval(Symbol(type))
    workspace[name] = dict2instance(instance_type, fn)
end

@unpack distributions = input
for dist in distributions
    @unpack name, type = dist
    instance_type = eval(Symbol(type))
    workspace[name] = dict2instance(instance_type, distributions[1]; workspace)
end

# Get the model
json_model_dist = [v for (k, v) in workspace if v isa HadronicUnpolarizedIntensity] |> first


# Load the YAML model following the same pattern as run-Xic2pKpi.jl
println("Loading YAML model...")

# Load YAML files
particledict = YAML.load_file(joinpath(dirname(@__DIR__), "data", "xic-particle-definitions.yaml"))
modelparameters = YAML.load_file(joinpath(dirname(@__DIR__), "data", "xic-model-definitions.yaml"))

defaultparameters = modelparameters["Default amplitude model"]

# Parse model dictionaries and convert to standard convention
(; chains, couplings, isobarnames) = parse_model_dictionaries(defaultparameters; particledict)

# Set up the amplitude model with particle numbering
# 0: Xic, 1:p, 2:pi, 3:K
yaml_model = Lc2ppiKModel(; chains, couplings, isobarnames)

println("✓ YAML model loaded successfully")
println("  - Chains: $(length(chains))")
println("  - Couplings: $(length(couplings))")
println("  - Isobars: $(length(isobarnames))")

# Test 1: Simple amplitude and intensity calculation
function test_simple_amplitude()
    json_model = json_model_dist.model
    ms = masses(json_model)
    println("Masses: ", ms)

    # Define a specific point
    σ1 = 1.4  # GeV²
    σ2 = 3.2  # GeV²
    σs = Invariants(ms; σ1, σ2)

    # Check if the point is physical
    if Kibble(σs, ms^2) < 0
        # Calculate amplitude for this point
        amplitude_value_j = amplitude(json_model, σs; refζs = (1, 1, 1, 1))
        intensity_value_j = unpolarized_intensity(json_model, σs; refζs = (1, 1, 1, 1))
        amplitude_value_y = amplitude(yaml_model, σs)  # helicity values
        intensity_value_y = unpolarized_intensity(yaml_model, σs)

        println("σ1 = $σ1, σ2 = $σ2, σs = $σs")
        
        println("JSON Amplitude: $amplitude_value_j")
        println("YAML Amplitude: $amplitude_value_y")

        if(amplitude_value_j ≈ amplitude_value_y)
            println("✓ Amplitudes match")
        else
            println("✗ Amplitudes do not match")
        end
    else
        println("Point is not in physical phase space")
    end
end


function test_complex_lineshapes()
    json_model = json_model_dist.model

    ms = masses(json_model)

    # Define a specific point
    σ1 = 1.4  # GeV²
    σ2 = 3.2  # GeV²
    σs = Invariants(ms; σ1, σ2)

    
    print("Testing with σ1 = $σ1, σ2 = $σ2 and σs = $σs\n")



    println("=== Complex Lineshape values =")

    # Direkte Iteration über chains und names
    for (i, name) in enumerate(json_model.names)
        # Zugriff auf die DecayChain über model.chains[i]
        dcj = json_model.chains[i]
        dcy = yaml_model.chains[i]
        kj = dcj.k
        ky = dcy.k

        s_test = σs[kj]
        
        amplitudeindj = amplitude(json_model.chains[i], σs; refζs = (1, 1, 1, 1))
        amplitudeindy = amplitude(yaml_model.chains[i], σs)

        lineshape_valuej = dcj.Xlineshape(s_test)
        lineshape_valuey = dcy.Xlineshape(s_test)

        if(lineshape_valuej ≈ lineshape_valuey)
            println("✓ Lineshapes match for $name")
            println("  Lineshape(s=$s_test): $lineshape_valuej")
        else
            println("✗ Lineshapes do not match for $name")
            println("$name:", kj, ky)
            println("  Lineshape JSON(s=$s_test): $lineshape_valuej")
            println("  Lineshape YAML(s=$s_test): $lineshape_valuey")
            println("  |Lineshape JSON|²: $(abs2(lineshape_valuej))")
            println("  |Lineshape YAML|²: $(abs2(lineshape_valuey))")
            println("  Amplitude JSON: $amplitudeindj")
            println("  Amplitude YAML: $amplitudeindy")
            println()
        end
        
    end
end

# Run the tests
function main()
    println("=== Test 1: Simple Amplitude ===")
    test_simple_amplitude()
    println("\n" * "="^50 * "\n")

    

    println("=== Test 3: Complex Lineshapes ===")
    test_complex_lineshapes()

    
end

# Run if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

main()