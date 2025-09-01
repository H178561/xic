# Simple YAML vs JSON Model Validation
# Focus on successfully loading both models with proper type mapping

# Activate environment 
using Pkg
Pkg.activate(".")

using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays  
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.JSON
using YAML

println("XiC Model Validation - YAML vs JSON")
println("="^40)

# Step 1: Create YAML model successfully
println("1. Loading YAML model...")
try
    particledict = YAML.load_file("data/xic-particle-definitions.yaml")
    modelparameters = YAML.load_file("data/xic-model-definitions.yaml")
    defaultparameters = modelparameters["Default amplitude model"]
    
    (; chains, couplings, isobarnames) = parse_model_dictionaries(defaultparameters; particledict)
    yaml_model = Lc2ppiKModel(; chains, couplings, isobarnames)
    
    println("✅ YAML model loaded successfully!")
    println("   Chains: $(length(chains))")
    println("   Couplings: $(length(couplings))")
    
    # Test YAML model with single point
    ms = yaml_model.tbs.ms
    test_point = Invariants(ms, σ1 = 1.9^2, σ2 = 1.6^2)
    helicities = [1, 0, 0, 1]
    
    yaml_amplitude = amplitude(yaml_model, test_point, helicities)
    yaml_intensity = abs2(yaml_amplitude)
    
    println("   Test amplitude: $yaml_amplitude")
    println("   Test intensity: $yaml_intensity")
    
    global yaml_model, yaml_amplitude, yaml_intensity, test_point, helicities
    
catch e
    println("❌ Failed to load YAML model: $e")
    exit(1)
end

# Step 2: Analyze JSON structure 
println("\n2. Analyzing JSON structure...")
json_file = "data/xic2pKpi-model_new.json"
json_content = open(json_file) do io
    JSON.parse(io)
end

decay_description = json_content["distributions"][1]["decay_description"]
functions = json_content["functions"]

println("   JSON functions: $(length(functions))")
println("   JSON decay chains: $(length(decay_description["chains"]))")

# Print all JSON function types to understand the mapping needed
println("\n3. JSON function types analysis:")
unique_types = Set{String}()
for fn in functions
    push!(unique_types, fn["type"])
end

for type_name in sort(collect(unique_types))
    count = sum(1 for fn in functions if fn["type"] == type_name)
    println("   $type_name: $count functions")
end

# Step 3: Create type mapping for JSON -> Julia types
println("\n4. Creating type mapping...")

function map_json_type(json_type_name, fn_dict)
    """Map JSON type names to actual Julia types based on function properties"""
    
    if json_type_name == "BreitWigner"
        # Check if it has angular momentum parameter to determine specific type
        if haskey(fn_dict, "l")
            return BreitWignerMinL
        else
            return BreitWignerMinL  # Default fallback
        end
    elseif json_type_name == "BlattWeisskopf"
        return BlattWeisskopfFormFactor
    else
        # Try to find the type directly
        try
            return eval(Symbol(json_type_name))
        catch
            println("   ⚠️  Unknown type: $json_type_name")
            return nothing
        end
    end
end

# Step 4: Reconstruct JSON model with proper type mapping  
println("\n5. Reconstructing JSON model...")
workspace = Dict{String,Any}()
failed_functions = String[]

for fn in functions
    name = fn["name"]
    json_type = fn["type"]
    
    julia_type = map_json_type(json_type, fn)
    
    if julia_type === nothing
        push!(failed_functions, name)
        continue
    end
    
    try
        workspace[name] = dict2instance(julia_type, fn)
        println("   ✅ Created: $name ($(julia_type))")
    catch e
        push!(failed_functions, name)
        println("   ❌ Failed: $name - $e")
    end
end

println("\n   Successfully created: $(length(workspace)) functions")
println("   Failed to create: $(length(failed_functions)) functions")

if !isempty(failed_functions)
    println("   Failed functions: $(join(failed_functions, ", "))")
end

# Step 5: Try to create JSON model
println("\n6. Attempting to create JSON model...")
try
    json_model = dict2instance(ThreeBodyDecay, decay_description; workspace)
    println("✅ JSON model created successfully!")
    
    # Test JSON model with same point
    json_amplitude = amplitude(json_model, test_point, helicities)
    json_intensity = abs2(json_amplitude)
    
    println("   Test amplitude: $json_amplitude")
    println("   Test intensity: $json_intensity")
    
    # Compare results
    println("\n7. COMPARISON RESULTS:")
    println("="^40)
    
    amp_diff = abs(yaml_amplitude - json_amplitude)
    amp_rel_error = amp_diff / abs(yaml_amplitude)
    
    int_diff = abs(yaml_intensity - json_intensity)
    int_rel_error = int_diff / yaml_intensity
    
    println("Amplitude:")
    println("  YAML:  $yaml_amplitude") 
    println("  JSON:  $json_amplitude")
    println("  Relative error: $amp_rel_error")
    
    println("\nIntensity:")
    println("  YAML:  $yaml_intensity")
    println("  JSON:  $json_intensity") 
    println("  Relative error: $int_rel_error")
    
    if amp_rel_error < 1e-12
        println("\n🎉 EXCELLENT AGREEMENT! (< 1e-12)")
    elseif amp_rel_error < 1e-8  
        println("\n✅ Very good agreement (< 1e-8)")
    elseif amp_rel_error < 1e-4
        println("\n⚠️  Moderate agreement (< 1e-4)")
    else
        println("\n❌ Poor agreement (> 1e-4)")
    end
    
catch e
    println("❌ Failed to create JSON model: $e")
    
    # Still provide partial analysis
    println("\n7. PARTIAL ANALYSIS:")
    println("✅ YAML model working: amplitude = $yaml_amplitude")
    println("❌ JSON model failed: type mapping issues need resolution")
    println("\nType mapping problems detected:")
    println("- JSON uses generic 'BreitWigner' but Julia needs 'BreitWignerMinL'")
    println("- JSON uses 'BlattWeisskopf' but Julia needs 'BlattWeisskopfFormFactor'")
    println("- Some functions failed to be created in workspace")
end

println("\nValidation complete!")
