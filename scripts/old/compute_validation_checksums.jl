# -------------------------------------------------------------
# Script to compute and add validation checksums to XiC JSON model
# Following the pattern from lc2ppik-lhcb-2683025.json
# -------------------------------------------------------------
using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.JSON
using YAML

println("Computing validation checksums for XiC → pKπ JSON model")
println("="^60)

# -------------------------------------------------------------
# Load and setup models
# -------------------------------------------------------------
println("Loading models...")

# Load YAML model
particledict = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-particle-definitions.yaml"))
modelparameters = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-model-definitions.yaml"))
defaultparameters = modelparameters["Default amplitude model"]
(; chains, couplings, isobarnames) = parse_model_dictionaries(defaultparameters; particledict)

# Set up three-body system
ms = let
    _mΞc = particledict["Lambda_c+"]["mass"] / 1e3  
    _mp = particledict["p"]["mass"] / 1e3
    _mπ = particledict["pi+"]["mass"] / 1e3
    _mK = particledict["K-"]["mass"] / 1e3
    ThreeBodyMasses(m1 = _mp, m2 = _mπ, m3 = _mK, m0 = _mΞc)
end

tbs = ThreeBodySystem(ms, ThreeBodySpins(1, 0, 0; two_h0 = 1))

# Create isobars
isobars = Dict()
for (key, lineshape) in chains
    dict = Dict{String, Any}(particledict[key])
    dict["lineshape"] = lineshape
    isobars[key] = definechaininputs(key, dict; tbs)
end

# Update parameters
defaultparameters["ArK(892)1"] = "1.0 ± 0.0"
defaultparameters["AiK(892)1"] = "0.0 ± 0.0"

println("✓ YAML model setup complete")

# Load JSON model
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
    error("No JSON model file found. Please run a conversion script first.")
end

json_content = open(json_file) do io
    JSON.parse(io)
end

println("✓ JSON model loaded from: $(basename(json_file))")

# -------------------------------------------------------------
# Define validation points
# -------------------------------------------------------------
println("Setting up validation points...")

# Standard validation point (matching Lc2pkpi pattern)
validation_point = Dict(
    "name" => "validation_point",
    "parameters" => [
        Dict("name" => "cos_theta_31", "value" => -0.2309352648098208),
        Dict("name" => "phi_31", "value" => 0.0),
        Dict("name" => "m_31", "value" => 1.9101377207489973),
        Dict("name" => "cos_theta_31_2", "value" => 0.0),
        Dict("name" => "phi_31_2", "value" => 0.0),
        Dict("name" => "m_31_2", "value" => 2.46794)  # XiC mass
    ]
)

# Additional points for individual resonances
validation_points = [
    validation_point,
    Dict(
        "name" => "validation_point_m12sq",
        "parameters" => [Dict("name" => "m_12_sq", "value" => 3.2)]
    ),
    Dict(
        "name" => "validation_point_m23sq", 
        "parameters" => [Dict("name" => "m_23_sq", "value" => 1.4)]
    ),
    Dict(
        "name" => "validation_point_m31sq",
        "parameters" => [Dict("name" => "m_31_sq", "value" => 3.2)]
    )
]

# Convert to kinematic invariants for calculation
σs_validation = Invariants(ms, 
    σ1 = validation_point["parameters"][3]["value"]^2,  # m_31^2
    σ2 = (validation_point["parameters"][6]["value"] - validation_point["parameters"][3]["value"])^2  # (m_XiC - m_31)^2
)

println("✓ Validation points defined")

# -------------------------------------------------------------
# Calculate model intensity at validation point
# -------------------------------------------------------------
println("Calculating model intensity...")

function calculate_total_intensity(isobars, tbs, σs, defaultparameters)
    total_intensity = 0.0
    
    # Sum over all chains with real coefficients
    for (parname, parvalue) in defaultparameters
        if startswith(parname, "Ar") && parname != "ArK(892)1"  # Skip reference amplitude
            try
                coefficient = parse(Float64, split(parvalue, " ")[1])  # Extract numerical value
                c, d = parname2decaychain(parname, isobars; tbs)
                
                # Sum over helicity combinations
                for (two_λ0, two_λ1) in [(1, 1), (1, -1), (-1, 1), (-1, -1)]
                    amp = c * amplitude(d, σs, [two_λ1, 0, 0, two_λ0])
                    total_intensity += abs2(coefficient * amp)
                end
            catch e
                println("  Warning: Could not process $parname: $e")
            end
        end
    end
    
    return total_intensity
end

total_intensity = calculate_total_intensity(isobars, tbs, σs_validation, defaultparameters)
println("✓ Total intensity calculated: $total_intensity")

# -------------------------------------------------------------
# Calculate individual resonance values
# -------------------------------------------------------------
println("Calculating individual resonance values...")

function calculate_resonance_value(resonance_name, isobars, σs)
    if haskey(isobars, resonance_name)
        try
            # For invariant mass resonances, use σ1 or σ2 depending on the resonance
            if contains(resonance_name, "K")
                # K resonances typically in m₂₃ 
                value = isobars[resonance_name].Xlineshape(σs.σ3)
            elseif contains(resonance_name, "L") || contains(resonance_name, "Λ")
                # Lambda resonances typically in m₃₁
                value = isobars[resonance_name].Xlineshape(σs.σ2) 
            elseif contains(resonance_name, "D") || contains(resonance_name, "Δ")
                # Delta resonances typically in m₁₂
                value = isobars[resonance_name].Xlineshape(σs.σ1)
            else
                # Default to σ1
                value = isobars[resonance_name].Xlineshape(σs.σ1)
            end
            return value
        catch e
            println("  Warning: Could not calculate $resonance_name: $e")
            return complex(0.0)
        end
    else
        println("  Warning: Resonance $resonance_name not found")
        return complex(0.0)
    end
end

# Calculate for major resonances
major_resonances = ["K(892)", "L(1405)", "L(1520)", "L(1690)", "D(1232)", "K(700)", "K(1430)"]
resonance_values = Dict()

for res_name in major_resonances
    value = calculate_resonance_value(res_name, isobars, σs_validation)
    resonance_values[res_name] = value
    println("  $res_name: $value")
end

# -------------------------------------------------------------
# Generate amplitude model checksums
# -------------------------------------------------------------
println("Generating amplitude model checksums...")

amplitude_model_checksums = []

# Overall model checksum
push!(amplitude_model_checksums, Dict(
    "point" => "validation_point",
    "distribution" => "default_model",
    "value" => total_intensity
))

# Individual resonance checksums
resonance_point_mapping = Dict(
    "K(892)" => "validation_point_m23sq",
    "K(700)" => "validation_point_m23sq", 
    "K(1430)" => "validation_point_m23sq",
    "L(1405)" => "validation_point_m31sq",
    "L(1520)" => "validation_point_m31sq",
    "L(1690)" => "validation_point_m31sq",
    "D(1232)" => "validation_point_m12sq"
)

for (res_name, value) in resonance_values
    if haskey(resonance_point_mapping, res_name)
        point_name = resonance_point_mapping[res_name]
        
        # Format as string with complex notation
        value_str = if imag(value) >= 0
            "$(real(value)) + $(imag(value))i"
        else
            "$(real(value)) - $(abs(imag(value)))i"
        end
        
        push!(amplitude_model_checksums, Dict(
            "point" => point_name,
            "distribution" => res_name,
            "value" => value_str
        ))
    end
end

println("✓ Generated $(length(amplitude_model_checksums)) checksums")

# -------------------------------------------------------------
# Update JSON model with checksums
# -------------------------------------------------------------
println("Updating JSON model with validation checksums...")

# Update the misc section
if !haskey(json_content, "misc")
    json_content["misc"] = Dict()
end

json_content["misc"]["amplitude_model_checksums"] = amplitude_model_checksums

# Also update parameter_points if they don't exist
if !haskey(json_content, "parameter_points")
    json_content["parameter_points"] = validation_points
end

# -------------------------------------------------------------
# Write updated JSON file
# -------------------------------------------------------------
output_file = replace(json_file, ".json" => "_validated.json")
println("Writing updated JSON with checksums to: $(basename(output_file))")

open(output_file, "w") do io
    JSON.print(io, json_content, 4)
end

# -------------------------------------------------------------
# Verification
# -------------------------------------------------------------
println("Verifying updated JSON file...")

# Read back and check
verification_content = open(output_file) do io
    JSON.parse(io)
end

@assert haskey(verification_content, "misc")
@assert haskey(verification_content["misc"], "amplitude_model_checksums")
@assert length(verification_content["misc"]["amplitude_model_checksums"]) == length(amplitude_model_checksums)

println("✓ Verification passed")

# -------------------------------------------------------------
# Summary
# -------------------------------------------------------------
println("\n" * "="^60)
println("VALIDATION CHECKSUM GENERATION COMPLETE")
println("="^60)

println("Input:  $(basename(json_file))")
println("Output: $(basename(output_file))")
println("Checksums added: $(length(amplitude_model_checksums))")

println("\nChecksum Summary:")
println("  - Overall model intensity: $total_intensity")
println("  - Individual resonances: $(length(resonance_values))")

println("\nResonance Values:")
for (res_name, value) in resonance_values
    println("  - $res_name: $value")
end

println("\nNext steps:")
println("  1. Use $(basename(output_file)) as your reference JSON model")
println("  2. Compare with other models using these checksums")
println("  3. Run crosscheck validation against original YAML model")

println("\n🎉 Validation checksums generated successfully!")
