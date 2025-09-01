import Pkg
Pkg.activate(dirname(@__DIR__))  # Go one directory up to the main project
Pkg.instantiate()

using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.Parameters
using YAML
using Measurements
using Statistics
using QuadGK
using Plots
using StatsBase

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
model = Lc2ppiKModel(; chains, couplings, isobarnames)

println("✓ YAML model loaded successfully")
println("  - Chains: $(length(chains))")
println("  - Couplings: $(length(couplings))")
println("  - Isobars: $(length(isobarnames))")

# Test 1: Simple amplitude and intensity calculation
function test_simple_amplitude()
    # Check the model structure
    println("Model type: $(typeof(model))")
    println("Model fields: $(fieldnames(typeof(model)))")
    
    # Get masses from the first chain
    ms = masses(model)
    println("Masses: ", ms)

    # Define a specific point
    σ1 = 1.4  # GeV²
    σ2 = 3.2  # GeV²
    σs = Invariants(ms; σ1, σ2)

    # Check if the point is physical
    if isphysical(σs, ms)
        # Calculate amplitude for this point
        amplitude_value = amplitude(model, σs)  # helicity values
        intensity_value = unpolarized_intensity(model, σs)

        println("σ1 = $σ1, σ2 = $σ2, σs = $σs")
        println("Amplitude: $amplitude_value")
        println("Intensity: $intensity_value")
    else
        println("Point is not in physical phase space")
    end
end

# Test 2: Individual resonance contributions
function test_individual_resonances()
    ms = masses(model)

    # Define a physical point (use the same as the simple test)
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

    if isphysical(σs, ms)
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
            amplitude_val = amplitude(_model, σs)

            # Intensity for this resonance
            intensity_val = unpolarized_intensity(_model, σs)

            println("$name:")
            println("  Amplitude: $amplitude_val")
            println("  Intensity: $intensity_val")
            println()
        end

        # Total amplitude and intensity
        total_amplitude = amplitude(model, σs)
        total_intensity = unpolarized_intensity(model, σs)

        println("=== Total ===")
        println("Total Amplitude: $total_amplitude")
        println("Total Intensity: $total_intensity")

    else
        println("Point is not in physical phase space")
        println("Kibble test failed for this point")
    end
end

function test_complex_lineshapes()
    ms = masses(model)

    # Define a physical point (use the same as the simple test)
    σ1 = 1.4  # GeV²
    σ2 = 3.2  # GeV²
    σs = Invariants(ms; σ1, σ2)
    
    print("Testing with σ1 = $σ1, σ2 = $σ2 and σs = $σs\n")

    println("=== Complex Lineshape values ===")

    # Direct iteration over chains and names
    for (i, name) in enumerate(model.names)
        # Access the DecayChain via model.chains[i]
        dc = model.chains[i]
        k = dc.k

        s_test = σs[k]
        amplitude_ind = amplitude(model.chains[i], σs)

        # Get lineshape value
        lineshape_value = dc.Xlineshape(s_test)

        println("$name (k=$k):")
        println("  Lineshape(s=$s_test): $lineshape_value")
        println("  |Lineshape|²: $(abs2(lineshape_value))")
        println("  Real part: $(real(lineshape_value))")
        println("  Imag part: $(imag(lineshape_value))")
        println("  Amplitude: $amplitude_ind")
        println()
    end
end

function BreitWignerTest()
    # Test a Breit-Wigner lineshape directly
    m0 = 0.8955
    Γ = 0.047299

    s_test = 1.5
    
    # Use the BreitWigner from the loaded model workspace if available
    # Otherwise create a simple test
    println("=== Breit-Wigner Test ===")
    println("Testing Breit-Wigner with m0=$m0, Γ=$Γ at s=$s_test")
    
    # Simple Breit-Wigner calculation
    bw_value = 1 / (m0^2 - s_test - 1im * m0 * Γ)
    println("Simple Breit-Wigner value: $bw_value")
    println("  |BW|²: $(abs2(bw_value))")
    println("  Real part: $(real(bw_value))")
    println("  Imag part: $(imag(bw_value))")
end

# Test 3: Compare with specific point from run-Xic2pKpi.jl
function test_reference_point()
    println("=== Reference Point Test ===")
    
    # Use the same point as in run-Xic2pKpi.jl
    ms = masses(model)
    σs0 = Invariants(ms;
        σ1 = 0.7980703453578917,
        σ2 = 3.6486261122281745)
    
    println("Reference point from run-Xic2pKpi.jl:")
    println("  σ1 = $(σs0.σ1)")
    println("  σ2 = $(σs0.σ2)")
    
    if isphysical(σs0, ms)
        _I = unpolarized_intensity(model, σs0)
        _A = amplitude(model, σs0, [1, 0, 0, 1])
        
        println("  Intensity: $_I")
        println("  Amplitude: $_A")
        
        # Test alternative helicity specification
        _A_alt = amplitude(model, σs0, ThreeBodySpins(1, 0, 0; two_h0 = 1))
        println("  Amplitude (alt): $_A_alt")
        println("  Match: $(_A == _A_alt)")
    else
        println("  Point is not in physical phase space")
    end
end

# Test 4: Test specific resonances mentioned in the model
function test_specific_resonances()
    println("=== Specific Resonance Tests ===")
    
    ms = masses(model)
    σs0 = Invariants(ms; σ1 = 3.2, σ2 = 1.4)
    
    # Test specific resonances that appear in both models
    resonances_to_test = ["L(1520)", "D(1232)", "K(892)", "L(1405)"]
    
    for res_name in resonances_to_test
        # Find chains that contain this resonance
        matching_chains = [name for name in model.names if occursin(res_name, name)]
        
        if !isempty(matching_chains)
            println("Testing $res_name resonance:")
            for chain_name in matching_chains
                _model = model[model.names.==chain_name]
                
                if isphysical(σs0, ms)
                    amplitude_val = amplitude(_model, σs0, [1, 0, 0, 1])
                    intensity_val = unpolarized_intensity(_model, σs0)
                    
                    println("  $chain_name:")
                    println("    Amplitude: $amplitude_val")
                    println("    Intensity: $intensity_val")
                else
                    println("  Point not physical for $chain_name")
                end
            end
        else
            println("Resonance $res_name not found in model")
        end
        println()
    end
end

# Run the tests
function main()
    println("=== YAML Model Test Suite ===")
    println("Model loaded with $(length(model.chains)) chains")
    println()
    
    println("=== Test 1: Simple Amplitude ===")
    test_simple_amplitude()
    println("\n" * "="^50 * "\n")

    println("=== Test 2: Individual Resonances ===")
    test_individual_resonances()
    println("\n" * "="^50 * "\n")

    println("=== Test 3: Complex Lineshapes ===")
    test_complex_lineshapes()
    println("\n" * "="^50 * "\n")

    println("=== Test 4: Breit-Wigner ===")
    BreitWignerTest()
    println("\n" * "="^50 * "\n")
    
    println("=== Test 5: Reference Point ===")
    test_reference_point()
    println("\n" * "="^50 * "\n")
    
    println("=== Test 6: Specific Resonances ===")
    test_specific_resonances()
end

# Run if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

main()