import Pkg
Pkg.activate(dirname(@__DIR__))  # Geht ein Verzeichnis hoch zum Hauptprojekt
Pkg.instantiate()

Pkg.develop(path="/home/hvdsmagt/ThreeBodyDecays.jl/IO.jl/ThreeBodyDecaysIO.jl/")


using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.HadronicLineshapes
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
model_dist = [v for (k, v) in workspace if v isa HadronicUnpolarizedIntensity] |> first

# Test 1: Simple amplitude and intensity calculation
function test_simple_amplitude()
    model = model_dist.model
    ms = masses(model)
    println("Masses: ", ms)

    # Define a specific point
    σ1 = 1.4  # GeV²
    σ2 = 3.2  # GeV²
    σs = Invariants(ms; σ1, σ2)

    # Check if the point is physical
    if Kibble(σs, ms^2) < 0
        # Calculate amplitude for this point
        amplitude_value = amplitude(model, σs; refζs = (1, 1, 1, 1))
        intensity_value = unpolarized_intensity(model, σs; refζs = (1, 1, 1, 1))

        println("σ1 = $σ1, σ2 = $σ2, σs = $σs")
        println("Amplitude: $amplitude_value")
        println("Intensity: $intensity_value")
    else
        println("Point is not in physical phase space")
    end
end

# Test 2: Individual resonance contributions
function test_individual_resonances()
    model = model_dist.model
    ms = masses(model)

    # Define a specific point
    σ1 = 1.4  # GeV²
    σ2 = 3.2  # GeV²
    σs = Invariants(ms; σ1, σ2)

  
    print("Testing with σ1 = $σ1, σ2 = $σ2 and σs = $σs\n")

    # Calculate σ3 explicitly
    m0² = ms.m0^2
    m1² = ms.m1^2
    m2² = ms.m2^2
    m3² = ms.m3^2
    σ3 = m0² + m1² + m2² + m3² - σ1 - σ2

    if Kibble(σs, ms^2) < 0
        println("=== Test Point ===")
        println("σ1 = $σ1, σ2 = $σ2, σ3 = $σ3")
        println("Masses: m0=$(ms.m0), m1=$(ms.m1), m2=$(ms.m2), m3=$(ms.m3)")
        println()

        # Test individual resonances
        chain_names = Set(model.names) |> collect |> sort

        println("=== Individual Resonance Contributions ===")
        for name in chain_names
            _model = model[model.names.==name]

            # Amplitude for this resonance
            amplitude_val = amplitude(_model, σs; refζs = (1, 1, 1, 1))


            # Intensity for this resonance
            intensity_val = unpolarized_intensity(_model, σs; refζs = (1, 1, 1, 1))

            println("$name:")
            println("  Amplitude: $amplitude_val")
            println("  Intensity: $intensity_val")
            println()
        end

        # Total amplitude and intensity
        total_amplitude = amplitude(model, σs; refζs = (1, 1, 1, 1))
        total_intensity = unpolarized_intensity(model, σs; refζs = (1, 1, 1, 1))

        println("=== Total ===")
        println("Total Amplitude: $total_amplitude")
        println("Total Intensity: $total_intensity")

    else
        println("Point is not in physical phase space")
        println("Kibble = $(Kibble(σs, ms^2))")
    end
end

function test_complex_lineshapes()
    model = model_dist.model

    ms = masses(model)

    # Define a specific point
    σ1 = 1.4  # GeV²
    σ2 = 3.2  # GeV²
    σs = Invariants(ms; σ1, σ2)

    
    print("Testing with σ1 = $σ1, σ2 = $σ2 and σs = $σs\n")



    println("=== Complex Lineshape values =")

    # Direkte Iteration über chains und names
    for (i, name) in enumerate(model.names)
        # Zugriff auf die DecayChain über model.chains[i]
        dc = model.chains[i]
        k = dc.k

        s_test = σs[k]
        amplitudeind = amplitude(model.chains[i], σs; refζs = (1, 1, 1, 1))

        lineshape_value = dc.Xlineshape(s_test)

        println("$name:", k)
        println("  Lineshape(s=$s_test): $lineshape_value")
        println("  |Lineshape|²: $(abs2(lineshape_value))")
        println("  Real part: $(real(lineshape_value))")
        println("  Imag part: $(imag(lineshape_value))")
        println("  Amplitude: $amplitudeind")
    end
end

function BreitWignerTest()

    m0 = 0.8955
    Γ = 0.047299

    s_test = 1.5
    bw = BreitWigner(m0, Γ)
    bw_value = bw(s_test)
    println("Breit-Wigner value at s = $s_test GeV²: $bw_value")

end

# Run the tests
function main()
    println("=== Test 1: Simple Amplitude ===")
    test_simple_amplitude()
    println("\n" * "="^50 * "\n")

    println("=== Test 2: Individual Resonances ===")
    test_individual_resonances()
    println("\n" * "="^50 * "\n")

    println("=== Test 3: Complex Lineshapes ===")
    test_complex_lineshapes()

    println("=== Test 4: Breit-Wigner ===")
    BreitWignerTest()
end

# Run if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

main()