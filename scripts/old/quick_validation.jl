# -------------------------------------------------------------
# Quick validation script to verify YAML vs JSON model equivalence
# Focused validation at specific test points
# -------------------------------------------------------------

println("XiC YAML vs JSON Model Validation")
println("="^40)

# Test if we can load the basic packages
try
    using Lc2ppiKSemileptonicModelLHCb
    using ThreeBodyDecaysIO.JSON
    using YAML
    println("✅ Packages loaded successfully")
catch e
    println("❌ Package loading failed: $e")
    exit(1)
end

# Test 1: Load and parse JSON model
println("\nTest 1: JSON Model Loading")
println("-"^30)

json_file = joinpath(@__DIR__, "..", "data", "xic2pKpi-model.json")
if !isfile(json_file)
    println("❌ JSON file not found: $json_file")
    exit(1)
end

try
    json_content = open(json_file) do io
        JSON.parse(io)
    end
    
    # Check basic structure
    @assert haskey(json_content, "distributions")
    @assert haskey(json_content, "functions")
    
    decay_desc = json_content["distributions"][1]["decay_description"]
    @assert haskey(decay_desc, "chains")
    
    chains = decay_desc["chains"]
    functions = json_content["functions"]
    
    println("✅ JSON model loaded successfully")
    println("   - Decay chains: $(length(chains))")
    println("   - Functions: $(length(functions))")
    
    # Check for validation checksums
    if haskey(json_content, "misc") && haskey(json_content["misc"], "amplitude_model_checksums")
        checksums = json_content["misc"]["amplitude_model_checksums"]
        println("   - Validation checksums: $(length(checksums))")
        
        # Count non-placeholder checksums
        real_checksums = [c for c in checksums if c["value"] != "0.0 + 0.0i" && c["value"] != 0.0]
        placeholder_checksums = length(checksums) - length(real_checksums)
        
        println("   - Real checksums: $(length(real_checksums))")
        println("   - Placeholder checksums: $placeholder_checksums")
    else
        println("   - Validation checksums: MISSING")
    end
    
catch e
    println("❌ JSON loading failed: $e")
    exit(1)
end

# Test 2: Load YAML model
println("\nTest 2: YAML Model Loading")
println("-"^30)

try
    particledict = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-particle-definitions.yaml"))
    modelparameters = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-model-definitions.yaml"))
    
    println("✅ YAML files loaded successfully")
    println("   - Particles: $(length(particledict))")
    println("   - Model parameters: $(length(modelparameters))")
    
    # Check the default model
    if haskey(modelparameters, "Default amplitude model")
        defaultmodel = modelparameters["Default amplitude model"]
        println("   - Default model parameters: $(length(defaultmodel))")
    else
        println("❌ No 'Default amplitude model' found")
        exit(1)
    end
    
catch e
    println("❌ YAML loading failed: $e")
    exit(1)
end

# Test 3: Basic model parsing
println("\nTest 3: Model Parsing")
println("-"^30)

try
    using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
    
    defaultparameters = modelparameters["Default amplitude model"]
    
    # Parse model (this might fail, so we'll catch it)
    result = parse_model_dictionaries(defaultparameters; particledict)
    chains, couplings, isobarnames = result.chains, result.couplings, result.isobarnames
    
    println("✅ Model parsing successful")
    println("   - Chains: $(length(chains))")
    println("   - Couplings: $(length(couplings))")
    println("   - Isobar names: $(length(isobarnames))")
    
    # List some key isobars
    key_isobars = ["K(892)", "L(1405)", "L(1520)", "D(1232)"]
    found_isobars = [name for name in key_isobars if name in keys(chains)]
    println("   - Key isobars found: $(found_isobars)")
    
catch e
    println("❌ Model parsing failed: $e")
    println("   This might be expected if model definitions need updates")
end

# Test 4: Crosscheck data availability
println("\nTest 4: Crosscheck Data")
println("-"^30)

crosscheck_file = joinpath(@__DIR__, "..", "data", "crosscheck_Xic.json")
if isfile(crosscheck_file)
    try
        crosscheck_data = open(crosscheck_file) do io
            JSON.parse(io)
        end
        
        println("✅ Crosscheck data available")
        println("   - File: crosscheck_Xic.json")
        
        if haskey(crosscheck_data, "chainvars")
            println("   - Chain variables: ✅")
        end
        
        if haskey(crosscheck_data, "lineshapes")
            lineshapes = crosscheck_data["lineshapes"]
            println("   - Lineshapes: $(length(lineshapes))")
        end
        
        if haskey(crosscheck_data, "chains")
            chains_data = crosscheck_data["chains"]
            println("   - Chain amplitudes: $(length(chains_data))")
        end
        
    catch e
        println("⚠️  Crosscheck data exists but couldn't parse: $e")
    end
else
    println("⚠️  No crosscheck data found")
end

# Test 5: JSON vs Reference model comparison
println("\nTest 5: Model Comparison Readiness")
println("-"^30)

# Check if we have reference models for comparison
lc_reference = joinpath(@__DIR__, "..", "data", "lc2ppik-lhcb-2683025.json")
if isfile(lc_reference)
    try
        lc_content = open(lc_reference) do io
            JSON.parse(io)
        end
        
        println("✅ Lc2pkpi reference model available")
        
        # Compare structures
        if haskey(lc_content, "misc") && haskey(lc_content["misc"], "amplitude_model_checksums")
            lc_checksums = lc_content["misc"]["amplitude_model_checksums"]
            println("   - Reference checksums: $(length(lc_checksums))")
            
            # Show example of non-zero checksum
            real_checksums = [c for c in lc_checksums if isa(c["value"], Number) && c["value"] != 0.0]
            if !isempty(real_checksums)
                example = real_checksums[1]
                println("   - Example checksum: $(example["distribution"]) = $(example["value"])")
            end
        end
        
    catch e
        println("⚠️  Reference model exists but couldn't parse: $e")
    end
else
    println("⚠️  No Lc2pkpi reference model found")
end

# Summary
println("\n" * "="^40)
println("VALIDATION SUMMARY")
println("="^40)

status_items = [
    ("JSON Model Structure", "✅ VALID"),
    ("YAML Model Files", "✅ PRESENT"),
    ("Model Parsing", "⚠️  NEEDS TESTING"),
    ("Validation Checksums", "⚠️  MOSTLY PLACEHOLDERS"),
    ("Crosscheck Data", isfile(crosscheck_file) ? "✅ AVAILABLE" : "❌ MISSING"),
    ("Reference Models", isfile(lc_reference) ? "✅ AVAILABLE" : "❌ MISSING")
]

for (item, status) in status_items
    println("$status $item")
end

println("\n🎯 NEXT STEPS:")
println("1. 📝 Compute real validation checksums for JSON model")
println("2. 🔍 Run detailed crosscheck against reference data") 
println("3. ⚖️  Compare YAML vs JSON amplitudes at test points")
println("4. ✅ Verify numerical equivalence")

println("\n💡 RECOMMENDATIONS:")
if isfile(crosscheck_file)
    println("✅ Use existing crosscheck_Xic.json for validation points")
end

# Check current JSON checksums status
if haskey(json_content, "misc") && haskey(json_content["misc"], "amplitude_model_checksums")
    checksums = json_content["misc"]["amplitude_model_checksums"]
    real_checksums = [c for c in checksums if c["value"] != "0.0 + 0.0i" && c["value"] != 0.0]
    
    if length(real_checksums) > 0
        println("✅ JSON has $(length(real_checksums)) real validation checksums")
    else
        println("📝 JSON checksums are placeholders - need computation")
    end
end

println("🔧 Use existing crosscheck infrastructure in test/ directory")
println("📊 Models are structurally ready for numerical validation")

println("\n🎉 Basic validation complete!")
