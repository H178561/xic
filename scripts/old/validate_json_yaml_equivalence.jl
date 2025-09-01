#!/usr/bin/env julia

"""
Compute and validate checksums for the XiC JSON model to ensure 
it produces identical results to the YAML model.
"""

using Pkg
Pkg.activate(".")

using JSON3
using ThreeBodyDecaysIO
using ThreeBodyDecays
using Parameters

# Load the JSON model
json_file = "data/xic2pKpi-model.json"
println("Loading JSON model from: $json_file")

if !isfile(json_file)
    error("JSON model file not found: $json_file")
end

json_data = JSON3.read(read(json_file, String))
dpd_json = read_model(json_file)

println("Successfully loaded JSON model")
println("Number of decay chains: ", length(dpd_json.chains))

# Load the YAML model for comparison
particle_file = "data/xic-particle-definitions.yaml"
model_file = "data/xic-model-definitions.yaml"

println("\nLoading YAML model for comparison...")
dpd_yaml = read_model(particle_file, model_file)
println("Successfully loaded YAML model")
println("Number of decay chains: ", length(dpd_yaml.chains))

# Define test kinematics (same as in crosscheck)
σs = [1.9^2, 2.2^2, 2.5^2]  # Three kinematic points
m1sq, m2sq, m3sq = (0.938)^2, (0.494)^2, (0.140)^2  # p, K, π masses
m0sq = (2.468)^2  # Ξc mass

println("\nComputing amplitude matrices at test points...")

# Compute amplitudes at test points
amplitude_checksums = Complex{Float64}[]
all_results = []

for (i, σ1) in enumerate(σs)
    println("Computing at kinematic point $i: σ₁ = $σ1")
    
    # Compute σ2 and σ3 from σ1 using Dalitz constraints
    σ2 = m0sq + m1sq + m2sq + m3sq - σ1 - (sqrt(σ1) - sqrt(m3sq))^2
    σ3 = m0sq + m1sq + m2sq + m3sq - σ1 - σ2
    
    # Amplitude calculations
    amp_json = amplitude(dpd_json, (σ1, σ2))
    amp_yaml = amplitude(dpd_yaml, (σ1, σ2))
    
    # Store results
    push!(all_results, (
        σ1 = σ1,
        σ2 = σ2, 
        σ3 = σ3,
        amp_json = amp_json,
        amp_yaml = amp_yaml,
        ratio = abs(amp_yaml) > 1e-10 ? amp_json / amp_yaml : NaN,
        abs_diff = abs(amp_json - amp_yaml),
        rel_diff = abs(amp_yaml) > 1e-10 ? abs(amp_json - amp_yaml) / abs(amp_yaml) : NaN
    ))
    
    # Accumulate for checksum
    push!(amplitude_checksums, amp_json)
    
    # Print comparison
    @printf("  JSON amplitude:  %12.6f %+12.6fi\n", real(amp_json), imag(amp_json))
    @printf("  YAML amplitude:  %12.6f %+12.6fi\n", real(amp_yaml), imag(amp_yaml))
    @printf("  Ratio:           %12.6f %+12.6fi\n", real(all_results[end].ratio), imag(all_results[end].ratio))
    @printf("  Abs difference:  %12.6e\n", all_results[end].abs_diff)
    @printf("  Rel difference:  %12.6e\n", all_results[end].rel_diff)
    println()
end

# Compute overall checksum
overall_checksum = sum(amplitude_checksums)
println("Overall amplitude checksum: $(real(overall_checksum)) + $(imag(overall_checksum))i")

# Analysis of differences
println("\n" * "="^60)
println("VALIDATION ANALYSIS")
println("="^60)

max_abs_diff = maximum(r.abs_diff for r in all_results if !isnan(r.abs_diff))
max_rel_diff = maximum(r.rel_diff for r in all_results if !isnan(r.rel_diff))

println("Maximum absolute difference: $max_abs_diff")
println("Maximum relative difference: $max_rel_diff")

# Check if models are equivalent
tolerance_abs = 1e-12
tolerance_rel = 1e-10

if max_abs_diff < tolerance_abs && max_rel_diff < tolerance_rel
    println("✅ MODELS ARE EQUIVALENT within tolerances")
    println("   Absolute tolerance: $tolerance_abs")
    println("   Relative tolerance: $tolerance_rel")
    
    # Update JSON file with computed checksum
    println("\nUpdating JSON model with computed checksums...")
    
    # Read the raw JSON to preserve formatting
    json_text = read(json_file, String)
    json_parsed = JSON3.read(json_text)
    
    # Update the amplitude_model_checksums section
    if haskey(json_parsed, "amplitude_model_checksums")
        checksums = json_parsed["amplitude_model_checksums"]
        
        # Update with computed values
        checksums["kinematic_point_1"]["amplitude"] = Dict(
            "real" => real(amplitude_checksums[1]),
            "imag" => imag(amplitude_checksums[1])
        )
        checksums["kinematic_point_2"]["amplitude"] = Dict(
            "real" => real(amplitude_checksums[2]),
            "imag" => imag(amplitude_checksums[2])
        )
        checksums["kinematic_point_3"]["amplitude"] = Dict(
            "real" => real(amplitude_checksums[3]),
            "imag" => imag(amplitude_checksums[3])
        )
        checksums["overall_checksum"] = Dict(
            "real" => real(overall_checksum),
            "imag" => imag(overall_checksum)
        )
        
        # Write updated JSON
        backup_file = json_file * ".backup"
        cp(json_file, backup_file)
        println("Created backup: $backup_file")
        
        open(json_file, "w") do io
            JSON3.pretty(io, json_parsed, JSON3.AlignmentContext(indent=2))
        end
        
        println("✅ Updated $json_file with computed checksums")
    else
        println("⚠️  No amplitude_model_checksums section found in JSON")
    end
    
else
    println("❌ MODELS ARE NOT EQUIVALENT")
    println("   Differences exceed tolerances")
    
    # Detailed analysis of specific differences
    println("\nDetailed difference analysis:")
    for (i, r) in enumerate(all_results)
        if !isnan(r.rel_diff) && r.rel_diff > tolerance_rel
            println("  Point $i: rel_diff = $(r.rel_diff)")
            println("    JSON: $(r.amp_json)")
            println("    YAML: $(r.amp_yaml)")
        end
    end
end

println("\n" * "="^60)
println("VALIDATION COMPLETE")
println("="^60)
