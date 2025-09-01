import Pkg
Pkg.activate(dirname(@__DIR__))  # Geht ein Verzeichnis hoch zum Hauptprojekt
Pkg.instantiate()

using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
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
using YAML

# Load the JSON model
input = open(
    #joinpath("/home/hvdsmagt/ThreeBodyDecays.jl/IO.jl/", "lc2ppik-lhcb-2683025.json"),
    joinpath(dirname(@__DIR__), "data", "test_json.json")
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
        if(name!="L(1600)")
            continue
        end
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
        # 0.045492
        # ✅ BERECHNE Form-Faktoren bei der AKTUELLEN Energie (nicht nur am Pol!)
    # Extract actual parameters from YAML model
        yaml_lineshape = dcy.Xlineshape
        m_res = yaml_lineshape.pars[1]  # Mass from YAML
        Γ₀ = yaml_lineshape.pars[2]     # Width from YAML  
        l = yaml_lineshape.l            # Angular momentum
        minL = yaml_lineshape.minL      # Minimum L
        
        # Extract parameters from JSON model  
        json_lineshape = dcj.Xlineshape
        json_m = json_lineshape.m
        json_channel = json_lineshape.channels[1]
        json_gsq = json_channel.gsq
        
        sqrt_s_test = sqrt(s_test)
        
        # Current momenta at s_test
        p_current = breakup(sqrt_s_test, ms[1], ms[3])  # p-K breakup
        q_current = breakup(ms[4], sqrt_s_test, ms[2])  # Λc decay
        
        # Reference momenta at resonance pole
        p0 = breakup(m_res, ms[1], ms[3])  # p-K at pole
        q0 = breakup(ms[4], m_res, ms[2])  # Λc at pole
        
        # Form factors (only for l > 0)
        dR = 1.5    # resonance radius
        dΛc = 5.0   # decay radius
        FF_resonance = HadronicLineshapes.BlattWeisskopf{l}(dR)
        FF_decay = HadronicLineshapes.BlattWeisskopf{minL}(dΛc)
        
        # YAML formula breakdown:
        # 1. Momentum-dependent width
        Γ_current = Γ₀ * (p_current/p0)^(2*l + 1) * m_res/sqrt_s_test * 
                    FF_resonance(p_current)^2 / FF_resonance(p0)^2
        
        # 2. Basic Breit-Wigner
        basic_BW = 1 / (m_res^2 - s_test - 1im * m_res * Γ_current)
        
        # 3. Additional kinematic factors
        momentum_factors = (p_current/p0)^l * (q_current/q0)^minL
        
        # 4. Form factor corrections
        ff_factors = sqrt(FF_resonance(p_current)^2 / FF_resonance(p0)^2 * 
                         FF_decay(q_current)^2 / FF_decay(q0)^2)
        
        # Complete YAML reconstruction
        yaml_reconstructed = basic_BW * momentum_factors * ff_factors
        
        # JSON should implement the same physics but with effective coupling
        # Calculate what the JSON should give with correct gsq
        p_json = breakup(sqrt_s_test, json_channel.ma, json_channel.mb)
        p0_json = breakup(json_m, json_channel.ma, json_channel.mb)
        
        # JSON width calculation (from MultichannelBreitWigner)
        FF_json = HadronicLineshapes.BlattWeisskopf{json_channel.l}(json_channel.d)
        momentum_factor_json = (p_json / p0_json)^(2*json_channel.l + 1)
        ff_ratio_json = FF_json(p_json)^2 / FF_json(p0_json)^2
        
        mΓ_json = json_gsq * 2*p_json / sqrt_s_test * momentum_factor_json * ff_ratio_json
        json_reconstructed = 1 / (json_m^2 - s_test - 1im * json_m * mΓ_json / json_m)
        
        println("$name analysis:")
        println("  s_test = $s_test GeV²")
        println("  YAML parameters: m = $m_res, Γ₀ = $Γ₀, l = $l, minL = $minL")
        println("  JSON parameters: m = $json_m, gsq = $json_gsq")
        println("  Current momenta: p = $p_current, q = $q_current")
        println("  Reference momenta: p₀ = $p0, q₀ = $q0")
        println("  Γ_current = $Γ_current")
        println("  momentum_factors = $momentum_factors")
        println("  ff_factors = $ff_factors")
        println("  Basic BW = $basic_BW")
        println("  YAML reconstructed = $yaml_reconstructed")
        println("  YAML actual = $lineshape_valuey")
        println("  JSON reconstructed = $json_reconstructed")
        println("  JSON actual = $lineshape_valuej")
        
        # Check reconstructions
        yaml_error = abs(yaml_reconstructed - lineshape_valuey) / abs(lineshape_valuey) * 100
        json_error = abs(json_reconstructed - lineshape_valuej) / abs(lineshape_valuej) * 100
        
        println("  YAML reconstruction error: $(yaml_error)%")
        println("  JSON reconstruction error: $(json_error)%")
        

       
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