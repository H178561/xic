# -------------------------------------------------------------
# Compare Complete YAML vs JSON Models at Single Dalitz Point
# Creates full models from both descriptions and compares total amplitude
# -------------------------------------------------------------

# Activate project environment
using Pkg
Pkg.activate(".")

using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.JSON
using YAML
using Test

# Import specific types for JSON reconstruction
using ThreeBodyDecaysIO.ThreeBodyDecays: BreitWignerMinL, BreitWigner, BlattWeisskopfFormFactor, BlattWeisskopf, ConstantAmplitude

println("XiC → pKπ Complete Model Comparison")
println("="^40)

# -------------------------------------------------------------
# Create complete YAML model
# -------------------------------------------------------------
println("Creating complete YAML model...")

# Load YAML files
particledict = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-particle-definitions.yaml"))
modelparameters = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-model-definitions.yaml"))
defaultparameters = modelparameters["Default amplitude model"]

# Parse YAML model
(; chains, couplings, isobarnames) = parse_model_dictionaries(defaultparameters; particledict)

# Create complete YAML model using the standard constructor
yaml_model = Lc2ppiKModel(; chains, couplings, isobarnames)

println("✓ Complete YAML model created")
println("  - Chains: $(length(chains))")
println("  - Couplings: $(length(couplings))")
println("  - Isobars: $(length(isobarnames))")

# -------------------------------------------------------------
# Load and reconstruct complete JSON model
# -------------------------------------------------------------
println("\nLoading complete JSON model...")

# Find JSON file
json_file = joinpath(@__DIR__, "..", "data", "xic2pKpi-model_new.json")
if !isfile(json_file)
    error("JSON model file not found: $json_file")
end

println("Using JSON file: $(basename(json_file))")

# Read JSON content
json_content = open(json_file) do io
    JSON.parse(io)
end

# Extract model components
decay_description = json_content["distributions"][1]["decay_description"]
functions = json_content["functions"]

println("  - Functions in JSON: $(length(functions))")

# Create workspace from functions with proper type mapping
workspace = Dict{String,Any}()
successful_functions = 0

for fn in functions
    name = fn["name"]
    type_str = fn["type"]
    
    try
        # Map JSON type names to actual workspace types based on parameters
        actual_type = if type_str == "BreitWigner"
            # Determine specific BreitWigner type based on angular momentum
            l = get(fn, "l", 0)
            if l == 0
                BreitWignerMinL
            else
                BreitWignerMinL  # Use MinL for all for now
            end
        elseif type_str == "BlattWeisskopf"
            BlattWeisskopfFormFactor
        elseif type_str == "Constant"
            ConstantAmplitude
        else
            # Try to find the type in available namespaces
            try
                eval(Symbol(type_str))
            catch
                println("  ⚠️  Unknown type '$type_str' for function '$name', skipping...")
                continue
            end
        end
        
        # Create the function instance
        workspace[name] = dict2instance(actual_type, fn)
        successful_functions += 1
        
    catch e
        println("  ⚠️  Warning: Could not create function '$name' of type '$type_str': $e")
        
        # Try alternative approaches
        if type_str == "BreitWigner"
            try
                # Try different BreitWigner variants
                for bw_type in [BreitWignerMinL, BreitWigner]
                    try
                        workspace[name] = dict2instance(bw_type, fn)
                        successful_functions += 1
                        println("  ✓ Successfully created '$name' using $bw_type")
                        break
                    catch
                        continue
                    end
                end
            catch e2
                println("  ❌ All BreitWigner attempts failed for '$name': $e2")
            end
        elseif type_str == "BlattWeisskopf"
            try
                # Try BlattWeisskopf variants
                for bw_type in [BlattWeisskopfFormFactor, BlattWeisskopf]
                    try
                        workspace[name] = dict2instance(bw_type, fn)
                        successful_functions += 1
                        println("  ✓ Successfully created '$name' using $bw_type")
                        break
                    catch
                        continue
                    end
                end
            catch e2
                println("  ❌ All BlattWeisskopf attempts failed for '$name': $e2")
            end
        end
    end
end

println("  - Successfully created $successful_functions out of $(length(functions)) functions")

# Reconstruct complete JSON model with error handling
println("\nReconstructing JSON model...")
try
    json_model = dict2instance(ThreeBodyDecay, decay_description; workspace)
    println("✓ Complete JSON model reconstructed successfully")
    json_model_created = true
catch e
    println("❌ Failed to reconstruct JSON model: $e")
    println("  Available functions in workspace: $(collect(keys(workspace)))")
    
    # Try to identify missing functions
    if haskey(decay_description, "chains")
        chains_data = decay_description["chains"]
        for chain_data in chains_data
            if haskey(chain_data, "lineshape") && !haskey(workspace, chain_data["lineshape"])
                println("  ❌ Missing function: $(chain_data["lineshape"])")
            end
        end
    end
    
    json_model = nothing
    json_model_created = false
end

# -------------------------------------------------------------
# Set up single test point on Dalitz plot
# -------------------------------------------------------------
println("\nSetting up test point...")

# Use masses from the model
ms = yaml_model.tbs.ms

# Choose a single test point in the physical region
test_point = Invariants(ms, σ1 = 1.9^2, σ2 = 1.6^2)

println("Test point:")
println("  σ₁ = $(test_point.σ1)")
println("  σ₂ = $(test_point.σ2)")

# Verify point is in physical region
if !isvalid(test_point)
    println("  ⚠️  Warning: Test point may be outside physical region")
else
    println("  ✓ Test point is in physical region")
end

# -------------------------------------------------------------
# Calculate total amplitude for both models
# -------------------------------------------------------------
println("\nCalculating total amplitudes...")

# Helicity configuration (unpolarized)
helicities = [1, 0, 0, 1]  # [λ₁, λ₂, λ₃, λ₀]

println("Using helicities: $helicities")

# Calculate YAML model amplitude
println("\nYAML model amplitude calculation...")
try
    yaml_amplitude = amplitude(yaml_model, test_point, helicities)
    println("  ✓ YAML amplitude: $yaml_amplitude")
    yaml_intensity = abs2(yaml_amplitude)
    println("  ✓ YAML intensity: $yaml_intensity")
catch e
    println("  ❌ Error calculating YAML amplitude: $e")
    yaml_amplitude = nothing
    yaml_intensity = nothing
end

# Calculate JSON model amplitude  
println("\nJSON model amplitude calculation...")
if json_model_created && json_model !== nothing
    try
        json_amplitude = amplitude(json_model, test_point, helicities)
        println("  ✓ JSON amplitude: $json_amplitude")
        json_intensity = abs2(json_amplitude)
        println("  ✓ JSON intensity: $json_intensity")
    catch e
        println("  ❌ Error calculating JSON amplitude: $e")
        json_amplitude = nothing
        json_intensity = nothing
    end
else
    println("  ❌ Cannot calculate JSON amplitude - model reconstruction failed")
    json_amplitude = nothing
    json_intensity = nothing
end

# -------------------------------------------------------------
# Compare results
# -------------------------------------------------------------
println("\n" * "="^40)
println("COMPARISON RESULTS")
println("="^40)

if yaml_amplitude !== nothing && json_amplitude !== nothing
    # Amplitude comparison
    amp_diff = abs(yaml_amplitude - json_amplitude)
    amp_rel_error = amp_diff / abs(yaml_amplitude)
    
    println("Amplitude Comparison:")
    println("  YAML:  $yaml_amplitude")
    println("  JSON:  $json_amplitude")
    println("  Difference: $amp_diff")
    println("  Relative error: $amp_rel_error")
    
    # Intensity comparison
    int_diff = abs(yaml_intensity - json_intensity)
    int_rel_error = int_diff / yaml_intensity
    
    println("\nIntensity Comparison:")
    println("  YAML:  $yaml_intensity")
    println("  JSON:  $json_intensity")
    println("  Difference: $int_diff")
    println("  Relative error: $int_rel_error")
    
    # Overall assessment
    println("\nOverall Assessment:")
    if amp_rel_error < 1e-12
        println("  ✅ EXCELLENT agreement (< 1e-12)")
        status = "EXCELLENT"
    elseif amp_rel_error < 1e-8
        println("  ✅ Very good agreement (< 1e-8)")
        status = "VERY_GOOD"
    elseif amp_rel_error < 1e-4
        println("  ⚠️  Moderate agreement (< 1e-4)")
        status = "MODERATE"
    else
        println("  ❌ Poor agreement (> 1e-4)")
        status = "POOR"
    end
    
    # Numerical test
    @testset "YAML vs JSON Model Equivalence" begin
        @test amp_rel_error < 1e-8 "Amplitude relative error should be < 1e-8"
        @test int_rel_error < 1e-8 "Intensity relative error should be < 1e-8"
    end
    
else
    println("❌ Could not compare - one or both amplitude calculations failed")
    status = "FAILED"
    amp_rel_error = NaN
    int_rel_error = NaN
end

# -------------------------------------------------------------
# Save comparison results
# -------------------------------------------------------------
results = Dict(
    "timestamp" => string(now()),
    "test_point" => Dict(
        "sigma1" => test_point.σ1,
        "sigma2" => test_point.σ2
    ),
    "helicities" => helicities,
    "yaml_model" => Dict(
        "amplitude" => yaml_amplitude === nothing ? "FAILED" : string(yaml_amplitude),
        "intensity" => yaml_intensity === nothing ? "FAILED" : yaml_intensity,
        "chains" => length(chains),
        "couplings" => length(couplings)
    ),
    "json_model" => Dict(
        "amplitude" => json_amplitude === nothing ? "FAILED" : string(json_amplitude),
        "intensity" => json_intensity === nothing ? "FAILED" : json_intensity,
        "functions" => length(functions)
    ),
    "comparison" => Dict(
        "amplitude_relative_error" => amp_rel_error,
        "intensity_relative_error" => int_rel_error,
        "status" => status
    )
)

results_file = joinpath(@__DIR__, "..", "complete_model_comparison.json")
open(results_file, "w") do io
    JSON.print(io, results, 4)
end

println("\nResults saved to: $(basename(results_file))")
println("\nComparison complete! 🎉")
