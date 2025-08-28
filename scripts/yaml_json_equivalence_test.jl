# YAML vs JSON Model Equivalence Test
# This script loads both YAML and JSON model descriptions and validates they produce identical results

import Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.JSON
using YAML
using Test

# Import required types explicitly
using ThreeBodyDecaysIO.ThreeBodyDecays: BreitWignerMinL, BlattWeisskopfFormFactor, ConstantAmplitude

println("XiC → pKπ YAML vs JSON Model Equivalence Test")
println("="^50)

# ============================================================================
# Load YAML Model
# ============================================================================
println("Loading YAML model...")

particledict = YAML.load_file(joinpath(dirname(@__DIR__), "data", "xic-particle-definitions.yaml"))
modelparameters = YAML.load_file(joinpath(dirname(@__DIR__), "data", "xic-model-definitions.yaml"))
defaultparameters = modelparameters["Default amplitude model"]

(; chains, couplings, isobarnames) = parse_model_dictionaries(defaultparameters; particledict)
yaml_model = Lc2ppiKModel(; chains, couplings, isobarnames)

println("✓ YAML model loaded successfully")
println("  - Chains: $(length(chains))")
println("  - Couplings: $(length(couplings))")
println("  - Isobars: $(length(isobarnames))")

# ============================================================================
# Load JSON Model  
# ============================================================================
println("\nLoading JSON model...")

json_file = joinpath(dirname(@__DIR__), "data", "xic2pKpi-model_new.json")
if !isfile(json_file)
    error("JSON model file not found: $json_file")
end

json_content = open(json_file) do io
    JSON.parse(io)
end

decay_description = json_content["distributions"][1]["decay_description"]
functions = json_content["functions"]

println("Using JSON file: $(basename(json_file))")
println("  - Functions in JSON: $(length(functions))")

# Create workspace from functions
workspace = Dict{String,Any}()
successful_functions = 0

for fn in functions
    name = fn["name"]
    type_str = fn["type"]
    
    try
        # Map JSON type names to actual workspace types
        actual_type = if type_str == "BreitWigner"
            BreitWignerMinL
        elseif type_str == "BlattWeisskopf"
            BlattWeisskopfFormFactor
        elseif type_str == "Constant"
            ConstantAmplitude
        else
            eval(Symbol(type_str))
        end
        
        workspace[name] = dict2instance(actual_type, fn)
        successful_functions += 1
        
    catch e
        println("  ⚠️  Warning: Could not create function '$name' of type '$type_str'")
        # Continue without this function
    end
end

println("  - Successfully created $successful_functions out of $(length(functions)) functions")

# Reconstruct JSON model
json_model = nothing
try
    json_model = dict2instance(ThreeBodyDecay, decay_description; workspace)
    println("✓ JSON model reconstructed successfully")
catch e
    println("❌ Failed to reconstruct JSON model: $e")
    println("This comparison cannot proceed without a valid JSON model.")
    exit(1)
end

# ============================================================================
# Model Comparison Tests
# ============================================================================
println("\n" * "="^50)
println("EQUIVALENCE TESTING")
println("="^50)

# Test points for comparison
test_points = [
    # Reference point from run-Xic2pKpi.jl
    (σ1 = 0.7980703453578917, σ2 = 3.6486261122281745),
    # Additional test points
    (σ1 = 1.5, σ2 = 3.2),
    (σ1 = 2.0, σ2 = 2.5),
    (σ1 = 1.0, σ2 = 4.0)
]

helicity_configs = [
    [1, 0, 0, 1],    # Standard configuration
    [0, 1, 0, 1],    # Alternative helicities
    [-1, 0, 0, 1]    # Another configuration
]

println("Testing $(length(test_points)) points with $(length(helicity_configs)) helicity configurations...")

all_tests_passed = true
total_tests = 0
passed_tests = 0

# Use YAML masses for creating test points (both models should have same masses)
ms_yaml = masses(yaml_model)

# Check if masses match
println("\nMass comparison:")
ms_yaml = masses(yaml_model)
ms_json = masses(json_model)
println("  YAML masses: $ms_yaml")
println("  JSON masses: $ms_json")

mass_tolerance = 1e-10
# Compare masses as tuples/arrays
yaml_mass_array = [ms_yaml.m1, ms_yaml.m2, ms_yaml.m3, ms_yaml.m0]
json_mass_array = [ms_json.m1, ms_json.m2, ms_json.m3, ms_json.m0]
masses_match = all(abs(yaml_mass_array[i] - json_mass_array[i]) < mass_tolerance for i in 1:4)
if masses_match
    println("  ✓ Masses match within tolerance")
else
    println("  ❌ Masses do not match")
    all_tests_passed = false
end

println("\nAmplitude comparison:")

for (i, test_point) in enumerate(test_points)
    σs = Invariants(ms_yaml; test_point...)
    
    # Check if point is physical
    if !isphysical(σs, ms_yaml)
        println("Test point $i: SKIPPED (not in physical region)")
        continue
    end
    
    println("Test point $i: σ₁=$(test_point.σ1), σ₂=$(test_point.σ2)")
    
    point_passed = true
    
    for (j, helicities) in enumerate(helicity_configs)
        total_tests += 1
        
        try
            # Calculate amplitudes
            yaml_amplitude = amplitude(yaml_model, σs, helicities)
            json_amplitude = amplitude(json_model, σs, helicities)
            
            # Calculate intensities
            yaml_intensity = abs2(yaml_amplitude)
            json_intensity = abs2(json_amplitude)
            
            # Compare amplitudes
            amp_diff = abs(yaml_amplitude - json_amplitude)
            rel_error = amp_diff / abs(yaml_amplitude)
            
            # Compare intensities
            int_diff = abs(yaml_intensity - json_intensity)
            int_rel_error = int_diff / yaml_intensity
            
            tolerance = 1e-8
            
            if rel_error < tolerance && int_rel_error < tolerance
                println("  ✓ Helicity $j: PASS (rel_error: $(rel_error))")
                passed_tests += 1
            else
                println("  ❌ Helicity $j: FAIL")
                println("    YAML amplitude: $yaml_amplitude")
                println("    JSON amplitude: $json_amplitude")
                println("    Relative error: $rel_error")
                println("    Intensity rel error: $int_rel_error")
                point_passed = false
                all_tests_passed = false
            end
            
        catch e
            println("  ❌ Helicity $j: ERROR - $e")
            point_passed = false
            all_tests_passed = false
        end
    end
    
    if !point_passed
        println("  → Point $i overall: FAILED")
    else
        println("  → Point $i overall: PASSED")
    end
end

# ============================================================================
# Individual Resonance Comparison
# ============================================================================
println("\nIndividual resonance comparison:")

# Get common resonances (simplified test)
yaml_names = Set(yaml_model.names)
json_names = Set(json_model.names) 

common_names = intersect(yaml_names, json_names)
println("Common resonances: $(length(common_names))")

if length(common_names) > 0
    # Test a few common resonances at the reference point
    test_σs = Invariants(ms_yaml; σ1 = 0.7980703453578917, σ2 = 3.6486261122281745)
    
    resonance_tests = 0
    resonance_passed = 0
    
    for name in collect(common_names)[1:min(5, length(common_names))]  # Test first 5
        try
            yaml_res_model = yaml_model[yaml_model.names.==name]
            json_res_model = json_model[json_model.names.==name]
            
            yaml_res_amp = amplitude(yaml_res_model, test_σs, [1, 0, 0, 1])
            json_res_amp = amplitude(json_res_model, test_σs, [1, 0, 0, 1])
            
            res_diff = abs(yaml_res_amp - json_res_amp)
            res_rel_error = res_diff / abs(yaml_res_amp)
            
            resonance_tests += 1
            
            if res_rel_error < 1e-8
                println("  ✓ $name: PASS")
                resonance_passed += 1
            else
                println("  ❌ $name: FAIL (rel_error: $res_rel_error)")
            end
            
        catch e
            println("  ⚠️  $name: Cannot compare ($e)")
        end
    end
    
    println("Individual resonance tests: $resonance_passed/$resonance_tests passed")
else
    println("  ⚠️  No common resonance names found between models")
end

# ============================================================================
# Final Results
# ============================================================================
println("\n" * "="^50)
println("FINAL RESULTS")
println("="^50)

println("Overall tests: $passed_tests/$total_tests passed")

if all_tests_passed && passed_tests == total_tests
    println("🎉 SUCCESS: YAML and JSON models are numerically equivalent!")
    println("   ✓ All amplitude comparisons passed")
    println("   ✓ All intensity comparisons passed") 
    println("   ✓ Models produce identical results within tolerance")
    exit_code = 0
else
    println("❌ FAILURE: YAML and JSON models are NOT equivalent")
    println("   Some tests failed - models produce different results")
    exit_code = 1
end

println("\nTolerance used: 1e-8")
println("JSON file: $(basename(json_file))")
println("YAML files: xic-particle-definitions.yaml, xic-model-definitions.yaml")

# ============================================================================
# Test Suite Integration
# ============================================================================
@testset "YAML vs JSON Model Equivalence" begin
    @test masses_match "Model masses should match"
    @test passed_tests == total_tests "All amplitude tests should pass"
    @test all_tests_passed "Overall equivalence test should pass"
end

println("\n" * "="^50)
if exit_code == 0
    println("✅ YAML ↔ JSON validation COMPLETED SUCCESSFULLY")
else
    println("❌ YAML ↔ JSON validation FAILED")
end
println("="^50)

exit(exit_code)
