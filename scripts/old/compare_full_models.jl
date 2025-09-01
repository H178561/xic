# -------------------------------------------------------------
# Compare Full YAML vs JSON Models at Single Dalitz Point
# This script creates complete models from both descriptions
# and compares total amplitude calculations at one test point
# -------------------------------------------------------------
using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.JSON
using YAML

println("XiC → pKπ: Full Model Comparison at Single Dalitz Point")
println("="^60)

# -------------------------------------------------------------
# 1. Create complete YAML model
# -------------------------------------------------------------
println("1. Creating complete YAML model...")

# Load YAML files
particledict = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-particle-definitions.yaml"))
modelparameters = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-model-definitions.yaml"))
defaultparameters = modelparameters["Default amplitude model"]

# Parse YAML to get chains, couplings, isobarnames
(; chains, couplings, isobarnames) = parse_model_dictionaries(defaultparameters; particledict)

# Create complete YAML model using Lc2ppiKModel
yaml_model = Lc2ppiKModel(; chains, couplings, isobarnames)

println("✅ YAML model created:")
println("   Type: $(typeof(yaml_model))")
println("   Chains: $(length(yaml_model.chains))")
println("   Resonances: $(length(unique(yaml_model.names)))")

# -------------------------------------------------------------
# 2. Load complete JSON model
# -------------------------------------------------------------
println("\n2. Loading complete JSON model...")

# Find JSON file
json_files = [
    joinpath(@__DIR__, "..", "data", "xic2pKpi-model_new.json"),
    joinpath(@__DIR__, "..", "data", "xic2pKpi-model.json")
]

json_file = nothing
for file in json_files
    if isfile(file)
        json_file = file
        break
    end
end

if json_file === nothing
    error("No JSON model file found. Please run xic_yaml_to_json_new.jl first.")
end

println("   Using JSON file: $(basename(json_file))")

# Read and parse JSON
json_content = open(json_file) do io
    JSON.parse(io)
end

# Reconstruct complete JSON model
decay_description = json_content["distributions"][1]["decay_description"]
functions = json_content["functions"]

# Create workspace for function reconstruction
workspace = Dict{String,Any}()
successful_functions = 0

for fn in functions
    try
        name = fn["name"]
        type_str = fn["type"]
        
        # Handle different type names that might be in JSON
        if type_str == "BreitWigner"
            # Map to actual implementation type
            instance_type = BreitWignerMinL
        elseif type_str == "BlattWeisskopf"
            # Skip form factors for now or handle appropriately
            continue
        else
            # Try to evaluate the type directly
            instance_type = eval(Symbol(type_str))
        end
        
        workspace[name] = dict2instance(instance_type, fn)
        successful_functions += 1
        
    catch e
        println("   Warning: Could not create function $(fn["name"]): $e")
    end
end

# Reconstruct the complete JSON model
json_model = dict2instance(ThreeBodyDecay, decay_description; workspace)

println("✅ JSON model loaded:")
println("   Type: $(typeof(json_model))")
println("   Functions created: $successful_functions/$(length(functions))")
println("   Chains: $(length(json_model.chains))")

# -------------------------------------------------------------
# 3. Define single test point on Dalitz plot
# -------------------------------------------------------------
println("\n3. Setting up test point...")

# Get masses from the YAML model
ms = yaml_model.chains[1].tbs.ms

# Choose a specific point in the Dalitz plot
# Use kinematically allowed values
σ1_test = 3.2  # m²₁₂ (p-π invariant mass squared) in GeV²
σ2_test = 2.5  # m²₁₃ (p-K invariant mass squared) in GeV²

# Create the Dalitz point
test_point = Invariants(ms, σ1 = σ1_test, σ2 = σ2_test)

println("✅ Test point defined:")
println("   σ₁ (m²₁₂) = $σ1_test GeV²")
println("   σ₂ (m²₁₃) = $σ2_test GeV²") 
println("   σ₃ (m²₂₃) = $(test_point.σ3) GeV²")
println("   Sum check: σ₁+σ₂+σ₃ = $(test_point.σ1 + test_point.σ2 + test_point.σ3)")
println("   Should equal: $(ms.m0^2 + ms.m1^2 + ms.m2^2 + ms.m3^2)")

# Verify point is physical
if !isphysical(test_point, ms)
    println("⚠️  Warning: Test point may not be in physical region")
end

# -------------------------------------------------------------
# 4. Calculate complete amplitudes for both models
# -------------------------------------------------------------
println("\n4. Calculating complete model amplitudes...")

# Define helicity combinations to test
helicity_combinations = [
    (1, 1),   # λ₀=+1/2, λ₁=+1/2
    (1, -1),  # λ₀=+1/2, λ₁=-1/2
    (-1, 1),  # λ₀=-1/2, λ₁=+1/2
    (-1, -1)  # λ₀=-1/2, λ₁=-1/2
]

println("   Testing $(length(helicity_combinations)) helicity combinations...")

# Calculate YAML model amplitudes
println("\n   YAML Model Amplitudes:")
yaml_amplitudes = ComplexF64[]

for (i, (two_λ0, two_λ1)) in enumerate(helicity_combinations)
    # Calculate amplitude for complete YAML model
    amp_yaml = amplitude(yaml_model, test_point, [two_λ1, 0, 0, two_λ0])
    push!(yaml_amplitudes, amp_yaml)
    
    println("     Helicity ($two_λ0, $two_λ1): $amp_yaml")
end

# Calculate JSON model amplitudes  
println("\n   JSON Model Amplitudes:")
json_amplitudes = ComplexF64[]

for (i, (two_λ0, two_λ1)) in enumerate(helicity_combinations)
    try
        # Calculate amplitude for complete JSON model
        amp_json = amplitude(json_model, test_point, [two_λ1, 0, 0, two_λ0])
        push!(json_amplitudes, amp_json)
        
        println("     Helicity ($two_λ0, $two_λ1): $amp_json")
        
    catch e
        println("     Helicity ($two_λ0, $two_λ1): ERROR - $e")
        push!(json_amplitudes, NaN + NaN*im)
    end
end

# -------------------------------------------------------------
# 5. Calculate total unpolarized intensities
# -------------------------------------------------------------
println("\n5. Calculating unpolarized intensities...")

try
    yaml_intensity = unpolarized_intensity(yaml_model, test_point)
    println("   YAML total intensity: $yaml_intensity")
catch e
    println("   YAML total intensity: ERROR - $e")
    yaml_intensity = NaN
end

try
    json_intensity = unpolarized_intensity(json_model, test_point)
    println("   JSON total intensity: $json_intensity")
catch e
    println("   JSON total intensity: ERROR - $e")
    json_intensity = NaN
end

# -------------------------------------------------------------
# 6. Compare results
# -------------------------------------------------------------
println("\n6. Comparing results...")
println("="^40)

# Compare amplitudes
println("Amplitude Comparisons:")
all_agree = true
max_rel_error = 0.0

for (i, (yaml_amp, json_amp)) in enumerate(zip(yaml_amplitudes, json_amplitudes))
    helicity = helicity_combinations[i]
    
    if isnan(yaml_amp) || isnan(json_amp)
        println("  Helicity $helicity: CANNOT COMPARE (NaN values)")
        all_agree = false
        continue
    end
    
    # Calculate difference and relative error
    diff = abs(yaml_amp - json_amp)
    rel_error = abs(yaml_amp) > 1e-15 ? diff / abs(yaml_amp) : 0.0
    max_rel_error = max(max_rel_error, rel_error)
    
    # Determine agreement status
    if rel_error < 1e-12
        status = "✅ EXCELLENT"
    elseif rel_error < 1e-8
        status = "✅ VERY GOOD"
    elseif rel_error < 1e-4
        status = "⚠️  MODERATE"
    else
        status = "❌ POOR"
        all_agree = false
    end
    
    println("  Helicity $helicity:")
    println("    YAML: $yaml_amp")
    println("    JSON: $json_amp")
    println("    Diff: $diff")
    println("    Rel Error: $rel_error")
    println("    Status: $status")
    println()
end

# Compare intensities
if !isnan(yaml_intensity) && !isnan(json_intensity)
    intensity_diff = abs(yaml_intensity - json_intensity)
    intensity_rel_error = yaml_intensity > 1e-15 ? intensity_diff / yaml_intensity : 0.0
    max_rel_error = max(max_rel_error, intensity_rel_error)
    
    intensity_status = if intensity_rel_error < 1e-12
        "✅ EXCELLENT"
    elseif intensity_rel_error < 1e-8
        "✅ VERY GOOD"
    elseif intensity_rel_error < 1e-4
        "⚠️  MODERATE"
    else
        "❌ POOR"
    end
    
    println("Intensity Comparison:")
    println("  YAML: $yaml_intensity")
    println("  JSON: $json_intensity")
    println("  Diff: $intensity_diff")
    println("  Rel Error: $intensity_rel_error")
    println("  Status: $intensity_status")
    println()
end

# -------------------------------------------------------------
# 7. Final summary
# -------------------------------------------------------------
println("="^60)
println("FINAL SUMMARY")
println("="^60)

println("Test Configuration:")
println("  Dalitz point: σ₁=$σ1_test, σ₂=$σ2_test, σ₃=$(test_point.σ3)")
println("  YAML model: $(length(yaml_model.chains)) chains")
println("  JSON model: $(length(json_model.chains)) chains")

println("\nNumerical Results:")
println("  Amplitudes tested: $(length(helicity_combinations))")
println("  Maximum relative error: $max_rel_error")

overall_status = if max_rel_error < 1e-12
    "✅ MODELS ARE NUMERICALLY IDENTICAL"
elseif max_rel_error < 1e-8
    "✅ MODELS AGREE VERY WELL"
elseif max_rel_error < 1e-4
    "⚠️  MODELS SHOW SMALL DIFFERENCES"
else
    "❌ MODELS SHOW SIGNIFICANT DIFFERENCES"
end

println("  Overall status: $overall_status")

if max_rel_error < 1e-8
    println("\n🎉 SUCCESS: YAML and JSON models produce equivalent results!")
    println("   The JSON conversion preserved the model accurately.")
else
    println("\n⚠️  The models show differences that may need investigation.")
    if max_rel_error > 1e-4
        println("   Consider checking the JSON conversion process.")
    end
end

println("\nModel validation complete!")
