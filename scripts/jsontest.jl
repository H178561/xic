import Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.Parameters
using YAML
using Measurements
using QuadGK
using Printf
using Statistics
using DataFrames

# ============================================================
# Load model from YAML
# ============================================================
function load_model_from_yaml(model_name::String)
    particledict = YAML.load_file(joinpath(dirname(@__DIR__), "data", "xic-particle-definitions.yaml"))
    modelparameters = YAML.load_file(joinpath(dirname(@__DIR__), "data", "xic-model-definitions.yaml"))
    
    defaultparameters = modelparameters[model_name]
    
    # Fix fixed parameters by adding dummy uncertainties
    if haskey(defaultparameters, "parameters")
        params = defaultparameters["parameters"]
        for (param_name, param_value) in params
            if isa(param_value, String)
                if !contains(param_value, "±") && !contains(param_value, "+/-")
                    try
                        val = parse(Float64, strip(param_value))
                        params[param_name] = "$val ± 0.0001"
                    catch
                        # If parsing fails, leave as is
                    end
                end
            end
        end
    end
    
    # Parse model
    (; chains, couplings, isobarnames) = parse_model_dictionaries(defaultparameters; particledict)
    
    # Create model (0: Xic, 1:p, 2:pi, 3:K)
    model = Lc2ppiKModel(; chains, couplings, isobarnames)
    
    return model
end

# ============================================================
# Calculate Fit Fractions for each resonance
# ============================================================
function calculate_fit_fractions(model; num_points=10000)
    # Initialize
    ms = masses(model)
    chain_names = Set(model.names) |> collect |> sort
    
    println("Generating $num_points random points...")
    
    # Pre-create component models (this is the optimization!)
    println("Creating component models for each resonance...")
    component_models = Dict(name => model[model.names.==name] for name in chain_names)
    println("  Created $(length(component_models)) component models")
    
    # Storage for statistics
    all_intensities = Dict(name => Float64[] for name in chain_names)
    all_model_intensities = Float64[]
    
    # Generate random points across Dalitz plot
    x2 = rand(num_points, 2)
    println("Generated random points, creating invariants...")
    
    data = map(eachslice(x2; dims=1)) do (x, y)
        σ1 = lims1(ms)[1] + x * diff(lims1(ms) |> collect)[1]
        σ2 = lims2(ms)[1] + y * diff(lims2(ms) |> collect)[1]
        σs = Invariants(ms; σ1, σ2)
    end
    
    println("Created $(length(data)) invariant points, filtering...")
    
    # Filter physically allowed points
    filter!(data) do σs
        Kibble(σs, ms^2) < 0
    end
    
    println("Filtered to $(length(data)) physical points, calculating intensities...")
    
    # Calculate intensities for each event
    n_processed = 0
    for σs in data
        # Total model intensity
        model_intensity = unpolarized_intensity(model, σs; refζs = (1, 1, 1, 1))
        push!(all_model_intensities, model_intensity)
        
        # Component intensities (using pre-created models)
        for name in chain_names
            component_intensity = unpolarized_intensity(component_models[name], σs; refζs = (1, 1, 1, 1))
            push!(all_intensities[name], component_intensity)
        end
        
        n_processed += 1
        if n_processed % 1000 == 0
            println("  Processed $n_processed / $(length(data)) points...")
        end
    end
    
    # Print final statistics in same format as the C++ version
    println("Component intensities: N = $(length(all_model_intensities))")
    
    # Calculate and print statistics for each component
    for name in chain_names
        # Mean of absolute intensities
        mean_value = mean(all_intensities[name])
        model_mean = mean(all_model_intensities)
        
        # Standard error calculation
        variance = sum((val - mean_value)^2 for val in all_intensities[name]) / length(all_intensities[name])
        std_error = sqrt(variance / length(all_intensities[name]))
        
        # Print formatted output
        println("$name Mean: $(round((mean_value/model_mean)*100, digits=3)) ± $(round((std_error/model_mean)*100, digits=3)))")
    end
    
    # Calculate final fit fractions for return value
    _int_i = map(chain_names) do name
        _intensities = all_intensities[name]
        _value = mean(_intensities)
        _err = sqrt(cov(_intensities, _intensities) / length(_intensities))
        _value ± _err
    end
    _int0 = mean(all_model_intensities)
    ff = round.(_int_i ./ _int0 .* 100; digits=2)
    
    return DataFrame("Resonance" => chain_names, "Fit Fraction [%]" => ff)
end

# ============================================================
# Calculate Fit Fractions using Grid method (200x200)
# ============================================================
function calculate_fit_fractions_grid(model; n_bins=200)
    # Initialize
    ms = masses(model)
    chain_names = Set(model.names) |> collect |> sort
    
    println("Creating $n_bins × $n_bins grid over Dalitz plot...")
    
    # Pre-create component models
    println("Creating component models for each resonance...")
    component_models = Dict(name => model[model.names.==name] for name in chain_names)
    println("  Created $(length(component_models)) component models")
    
    # Storage for statistics
    all_intensities = Dict(name => Float64[] for name in chain_names)
    all_model_intensities = Float64[]
    
    # Get kinematic limits
    σ1_min, σ1_max = lims1(ms)
    σ2_min, σ2_max = lims2(ms)
    
    println("Kinematic limits:")
    println("  σ1: [$σ1_min, $σ1_max] GeV²")
    println("  σ2: [$σ2_min, $σ2_max] GeV²")
    
    # Create grid
    σ1_range = range(σ1_min, σ1_max, length=n_bins)
    σ2_range = range(σ2_min, σ2_max, length=n_bins)
    
    println("Calculating intensities on grid...")
    
    n_valid = 0
    n_total = n_bins * n_bins
    
    for (i, σ1) in enumerate(σ1_range)
        for (j, σ2) in enumerate(σ2_range)
            try
                σs = Invariants(ms; σ1, σ2)
                
                # Check if point is physical (inside Dalitz plot)
                if Kibble(σs, ms^2) < 0
                    n_valid += 1
                    
                    # Total model intensity
                    model_intensity = unpolarized_intensity(model, σs; refζs = (1, 1, 1, 1))
                    push!(all_model_intensities, model_intensity)
                    
                    # Component intensities (using pre-created models)
                    for name in chain_names
                        component_intensity = unpolarized_intensity(component_models[name], σs; refζs = (1, 1, 1, 1))
                        push!(all_intensities[name], component_intensity)
                    end
                end
            catch e
                # Skip invalid points
                continue
            end
        end
        
        # Progress indicator every 10 rows
        if i % 10 == 0
            println("  Processed row $i / $n_bins ($(n_valid) valid points so far)")
        end
    end
    
    println()
    println("Grid calculation complete:")
    println("  Total grid points: $n_total")
    println("  Valid physical points: $n_valid")
    println("  Efficiency: $(round(n_valid/n_total*100, digits=2))%")
    println()
    
    # Print final statistics
    println("Component intensities: N = $(length(all_model_intensities))")
    
    # Calculate and print statistics for each component
    for name in chain_names
        # Mean of absolute intensities
        mean_value = mean(all_intensities[name])
        model_mean = mean(all_model_intensities)
        
        # Standard error calculation
        variance = sum((val - mean_value)^2 for val in all_intensities[name]) / length(all_intensities[name])
        std_error = sqrt(variance / length(all_intensities[name]))
        
        # Print formatted output
        println("$name Mean: $(round((mean_value/model_mean)*100, digits=3)) ± $(round((std_error/model_mean)*100, digits=3)))")
    end
    
    # Calculate final fit fractions for return value
    _int_i = map(chain_names) do name
        _intensities = all_intensities[name]
        _value = mean(_intensities)
        _err = sqrt(cov(_intensities, _intensities) / length(_intensities))
        _value ± _err
    end
    _int0 = mean(all_model_intensities)
    ff = round.(_int_i ./ _int0 .* 100; digits=2)
    
    return DataFrame("Resonance" => chain_names, "Fit Fraction [%]" => ff)
end

# ============================================================
# Test amplitude at a single point
# ============================================================
function test_simple_amplitude(model)
    ms = masses(model)
    
    # Define a specific point
    σ1 = 1.4  # GeV²
    σ2 = 3.2  # GeV²
    σs = Invariants(ms; σ1, σ2)

    # Check if the point is physical
    if Kibble(σs, ms^2) < 0
        # Calculate amplitude for this point
        amplitude_value = amplitude(model, σs; refζs = (1, 1, 1, 1))
        intensity_value = unpolarized_intensity(model, σs; refζs = (1, 1, 1, 1))

        println("Test point: σ1 = $σ1 GeV², σ2 = $σ2 GeV², σs = $σs")
        println("Amplitude: $amplitude_value")
        println("Intensity: $intensity_value")
    else
        println("Point is not in physical phase space")
    end
end

function BreitWignerTest()
    m0 = 0.8955
    Γ = 0.047299
    s_test = 1.5
    
    bw = BreitWigner(m0, Γ)
    bw_value = bw(s_test)
    println("Breit-Wigner test:")
    println("  m0 = $m0 GeV, Γ = $Γ GeV")
    println("  BW($s_test GeV²) = $bw_value")
end

# Run the tests
function main()
    println("="^70)
    println("XiC → pKπ Amplitude Model Analysis")
    println("="^70)
    println()
    
    

    # list of models
    list = [
        #"Default amplitude model",
        "Alternative model 1 - Delta resonances with free mass and width",
        #"Alternative model 2 - K(700) Relativistic BW",
        #"Alternative model 3 - K(700) with free mass and width",
        #"Alternative model 4 - K(1430) Relativistic BW",
        #"Alternative model 5 - K(1430) with free mass and width",
        #"Alternative model 6 - Both K(700) and K(1430) Relativistic BW",
        #"Alternative model 7 - Both K(700) and K(1430) with free mass and width",
        #"Alternative model 8 - K(700) Relativistic BW with free mass and width",
        #"Alternative model 9 - Both K(700) and K(1430) Relativistic BW with free mass and width",
        #"Alternative model 10 - L(1800) resonance removed",
        #"Alternative model 11 - L(1890) resonance removed",
        #"Alternative model 12 - K(1430) m=1370 MeV, Γ=180 MeV",
        #"Alternative model 13 - K(1430) m=1370 MeV, Γ=360 MeV",
        #"Alternative model 14 - K(1430) m=1430 MeV, Γ=180 MeV",
        #"Alternative model 12 - K(1430) m=1370 MeV, Γ=180 MeV",
        #"Alternative model 15 - K(1430) m=1430 MeV, Γ=360 MeV",
        #"Alternative model 16 - Multiple K mass variations 1",
        #"Alternative model 17 - Multiple K mass variations 2",
        #"Alternative model 18 - Multiple K mass variations 3",
        #"Alternative model 19 - Multiple K mass variations 4",
        #"Alternative model 20 - L(1405) free Flatte widths",
        #"Alternative model 21 - L(1600) with free mass and width",
        #"Alternative model 22 - L(1710) with free mass and width",
        #"Alternative model 23 - L(1800) contribution added with free mass and width",
        #"Alternative model 24 - L(1830) contribution added with free width",
        #"Alternative model 25 - L(1890) with free mass and width",
        #"Alternative model 26 - L(2000) with free mass and width",
        #"Alternative model 27 - L(2100) contribution added with PDG values",
        #"Alternative model 28 - L(2110) contribution added with PDG values",
        #"Alternative model 29 - S(1670) contribution added with PDG values",
        #"Alternative model 30 - S(1775) contribution added with PDG values",
        #"Alternative model 31 - Free radial parameter rXic"
    ]

    # Load the default model
    println("📂 Loading default model from YAML...")
    #model = load_model_from_yaml("Default amplitude model")
    #model = load_model_from_yaml("Alternative model 1 - Delta resonances with free mass and width")

    for model_name in list
        model = load_model_from_yaml(model_name)
        println("✓ Model loaded successfully!", model_name)
        println("  - Number of decay chains: $(length(model.chains))")
        println("  - Number of unique resonances: $(length(unique(model.names)))")
        println("  - Resonances: $(sort(unique(model.names)))")
        println()
        #=
        # Calculate fit fractions
        println("="^70)
        println("CALCULATING FIT FRACTIONS")
        println("="^70)
        println()
        
        fit_fractions_df = calculate_fit_fractions(model; num_points=10000)
        
        println()
        println("="^70)
        println("FIT FRACTIONS RESULTS (Monte Carlo)")
        println("="^70)
        println(fit_fractions_df)
        println()

        
        # Calculate fit fractions using grid method
        println("="^70)
        println("CALCULATING FIT FRACTIONS (GRID METHOD)")
        println("="^70)
        println()
        
        
        
        =#

        fit_fractions_grid_df = calculate_fit_fractions_grid(model; n_bins=200)

        println()
        println("="^70)
        println("FIT FRACTIONS RESULTS (Grid 200x200)")
        println("="^70)
        println(fit_fractions_grid_df)
        println()

        # Test at a specific point
        println("="^70)
        println("SINGLE POINT TEST")
        println("="^70)
        println()
        
        test_simple_amplitude(model)
        println()
        
        # Individual contributions
        println("="^70)
        println("INDIVIDUAL RESONANCE CONTRIBUTIONS AT TEST POINT")
        println("="^70)
        println()
        
        test_individual_resonances_short(model)
    end

   
    
end

# Test 1: Simple amplitude and intensity calculation

# Shortened version of test_individual_resonances
function test_individual_resonances_short(model)
    ms = masses(model)
    σ1, σ2 = 1.4, 3.2
    σs = Invariants(ms; σ1, σ2)
    
    if Kibble(σs, ms^2) < 0
        chain_names = unique(model.names) |> collect |> sort
        
        # First calculate total intensity
        total_intensity = unpolarized_intensity(model, σs; refζs = (1, 1, 1, 1))
        
        # Print header
        @printf("%-25s | %20s | %20s\n", "Resonance", "Intensity", "Percentage")
        println("-"^70)
        
        # Calculate and print each component
        for name in chain_names
            _model = model[model.names.==name]
            intensity_val = unpolarized_intensity(_model, σs; refζs = (1, 1, 1, 1))
            percentage = (intensity_val / total_intensity) * 100
            @printf("%-25s | %20.6f | %19.3f%%\n", name, intensity_val, percentage)
        end
        
        println("-"^70)
        @printf("%-25s | %20.6f | %19.3f%%\n", "Total", total_intensity, 100.0)
        println()
        println("Test point: σ1 = $σ1 GeV², σ2 = $σ2 GeV²")
    end
end

# Run if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

main()