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

my_breakup(m², m1², m2²) = sqrtKallenFact(sqrt(m²), sqrt(m1²), sqrt(m2²)) / (2 * sqrt(m²))


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
        #s_test = 2.8
        
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
        


        #σ = sqrt(s_test)                
        σ = s_test

        m1, m2, mk = ms[3], ms[1], ms[2]
        m0 = ms[4]
        m, Γ₀ = yaml_lineshape.pars
        print(m, " ", Γ₀, " ", m1, " ", m2, " ", mk, " ", m0, "\n")
        l, minL = yaml_lineshape.l, yaml_lineshape.minL
        println(l, minL)
        p, p0 = my_breakup(σ, m1^2, m2^2), my_breakup(m^2, m1^2, m2^2)
        q, q0 = my_breakup(m0^2, σ, mk^2), my_breakup(m0^2, m^2, mk^2)
        Γ = Γ₀ * (p / p0)^(2l + 1) * m / sqrt(σ) * F²(l, p, p0, dR)
        res = 1 / (m^2 - σ - 1im * m * Γ) * (p / p0)^l * (q / q0)^minL *sqrt(F²(l, p, p0, dR) * F²(minL, q, q0, dΛc))

        println(res)

        bwminl = bwminl = BreitWignerMinL(
            pars = (m, Γ₀),
            l = l,
            minL = minL,
            name = "L(1600)",
            m1 = m1,
            m2 = m2,
            mk = mk,
            m0 = m0
        )
        lineshape_valuej2 = bwminl(σ)
        println("BreitWignerMinL result: ", lineshape_valuej2)
        #print all minl parameters
        println("BreitWignerMinL parameters: ", bwminl.pars, " ", bwminl.l, " ", bwminl.minL, " ", bwminl.m1, " ", bwminl.m2, " ", bwminl.mk, " ", bwminl.m0)
        
       
       # Ersetzen Sie die inkonsistenten breakup-Aufrufe:

        # Check reconstructions
        p = breakup(sqrt(σ), m1, m2)      # ← sqrt(σ) für Masse!
        p0 = breakup(m, m1, m2)           # ← m für Masse!
        q = breakup(m0, sqrt(σ), mk)      # ← sqrt(σ) für Masse!
        q0 = breakup(m0, m, mk)           # ← m für Masse!

        # VOLLSTÄNDIGE Width-Korrektur (wie BreitWignerMinL):
        FFp = BlattWeisskopf{l}(dR)
        FFq = BlattWeisskopf{minL}(dΛc)

        # Resonance form factor correction (für die Width)
        resonance_ff_correction = FFp(p)^2 / FFp(p0)^2

        # Decay form factor correction  
        decay_ff_correction = FFq(q)^2 / FFq(q0)^2

        # Momentum-abhängige Width (exakt wie BreitWignerMinL):
        Γ_full = Γ₀ * (p / p0)^(2*l + 1) * m / sqrt(σ) * resonance_ff_correction

        # Erstelle MultichannelBreitWigner mit der VOLLEN width:
        bw = MultichannelBreitWigner(m, Γ_full, m1, m2, l, dR)
        ls = bw(σ)

        # Amplitude correction (OHNE die FF, die schon in Γ_full sind):
        amp_correction = (p / p0)^l * (q / q0)^minL * sqrt(decay_ff_correction)

        final_result = amp_correction * ls
        println("Final reconstructed result: ", final_result)
  
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
        

        #=
        {
            "mk": 0.938272046,
            "m2": 0.49367700000000003,
            "x": "m_31_sq",
            "minL": 0,
            "m1": 0.13957018,
            "width": 0.27,
            "l": 0,
            "name": "K(1430)_BuggBWminL",
            "gamma": 0.051629,
            "mass": 1.425,
            "m0": 2.46794,
            "type": "BuggBreitWignerMinL"
        },
        =#

        buggbwminl = BuggBreitWignerMinL(
            pars = (1.425, 0.27, 0.051629),
            l = 0,
            minL = 0,
            m1 = 0.13957018,
            m2 = 0.49367700000000003,
            mk = 0.938272046,
            m0 = 2.46794,
            name = "K(1430)"
        )
        s = 1.4
        buggbwminl_value = buggbwminl(s)
        println("BuggBreitWignerMinL result: ", buggbwminl_value)


        #=
          {
            "mk": 0.13957018,
            "m2": 0.938272046,
            "x": "m_31_sq",
            "minL": 0,
            "m1": 0.49367700000000003,
            "width": 0.0505,
            "l": 0,
            "name": "L(1405)_Flatte1405",
            "mass": 1.4051,
            "m0": 2.46794,
            "type": "Flatte1405"
        },
        =#

        flatte1405 = Flatte1405(
            pars = (1.4051, 0.0505),
            l = 0,
            minL = 0,
            m1 = 0.49367700000000003,
            m2 = 0.938272046,
            mk = 0.13957018,
            m0 = 2.46794,
            name = "L(1405)"
        )
        s = 3.2
        flatte1405_value = flatte1405(s)
        println("Flatte1405 result: ", flatte1405_value)
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