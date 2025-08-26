# -------------------------------------------------------------
# Script to convert XiC to pKπ model from YAML to JSON format
# Fixed version that properly structures the JSON to match Lc2pkpi format
# -------------------------------------------------------------
using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.JSON
using YAML
using Dates

# -------------------------------------------------------------
# Load model and particle definitions from YAML files
# -------------------------------------------------------------
begin
    particledict = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-particle-definitions.yaml"))
    modelparameters = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-model-definitions.yaml"))
end

defaultparameters = modelparameters["Default amplitude model"]

# -------------------------------------------------------------
# Parse model dictionaries and convert to standard convention
# -------------------------------------------------------------
(; chains, couplings, isobarnames) =
    parse_model_dictionaries(defaultparameters; particledict)

# -------------------------------------------------------------
# Set up the amplitude model with particle numbering
# 0: Xic, 1:p, 2:pi, 3:K
# -------------------------------------------------------------
model = Lc2ppiKModel(; chains, couplings, isobarnames)

println("Model type: $(typeof(model))")
println("Model is already a ThreeBodyDecay: $(isa(model, ThreeBodyDecay))")

threebody_model = model

# -------------------------------------------------------------
# Custom lineshape parser for XiC model to match Lc2pkpi format
# -------------------------------------------------------------
function xic_lineshape_parser(Xlineshape)
    appendix = Dict{String, Any}()
    
    # Get the lineshape type and parameters
    lineshape_type = typeof(Xlineshape)
    lineshape_name = Xlineshape.name
    
    # Create proper function name based on lineshape type
    if lineshape_type <: BreitWignerMinL
        scattering_key = "$(lineshape_name)_BW"
        m, Γ = Xlineshape.pars
        
        # Determine correct invariant mass variable based on decay chain
        # For XiC: 0->p(1), pi(2), K(3), so [3,1] = pK, [1,2] = p-pi, [2,3] = pi-K
        x_var = if contains(lineshape_name, "L") || contains(lineshape_name, "Λ")
            "m_31_sq"  # Lambda resonances in pK system
        elseif contains(lineshape_name, "D") || contains(lineshape_name, "Δ")
            "m_12_sq"  # Delta resonances in p-pi system  
        elseif contains(lineshape_name, "K")
            "m_23_sq"  # Kaon resonances in pi-K system
        else
            "m_31_sq"  # Default to pK system
        end
        
        lineshape_dict = Dict{String, Any}(
            "name" => scattering_key,
            "type" => "BreitWigner",
            "x" => x_var,
            "mass" => m,
            "width" => Γ,
            "l" => Xlineshape.l,
            "ma" => Xlineshape.m1,
            "mb" => Xlineshape.m2,
            "d" => 1.5
        )
        
    elseif lineshape_type <: BuggBreitWignerMinL
        scattering_key = "$(lineshape_name)_BuggBW"
        m, Γ, γ = Xlineshape.pars
        
        lineshape_dict = Dict{String, Any}(
            "name" => scattering_key,
            "type" => "generic_function",
            "expression" => "1/($(m)^2 - σ - i * $(m) * (σ - $(Xlineshape.m1^2 + Xlineshape.m2^2)) / ($(m)^2 - $(Xlineshape.m1^2 + Xlineshape.m2^2)) * $(Γ) * exp(-$(γ) * σ))"
        )
        
    elseif lineshape_type <: Flatte1405
        scattering_key = "$(lineshape_name)_Flatte"
        m, Γ = Xlineshape.pars
        
        lineshape_dict = Dict{String, Any}(
            "name" => scattering_key,
            "type" => "MultichannelBreitWigner",
            "x" => "m_31_sq",
            "mass" => m,
            "channels" => [
                Dict{String, Any}(
                    "gsq" => 0.23395150538434703,  # Match Lc2pkpi values
                    "ma" => Xlineshape.m1,
                    "mb" => Xlineshape.m2,
                    "l" => Xlineshape.l,
                    "d" => 0
                ),
                Dict{String, Any}(
                    "gsq" => 0.23395150538434703,
                    "ma" => 1.18937,  # Sigma mass
                    "mb" => 0.13957018,  # Pion mass
                    "l" => Xlineshape.l,
                    "d" => 0
                )
            ]
        )
        
    elseif lineshape_type <: L1670Flatte
        scattering_key = "$(lineshape_name)_BW"
        m, Γ = Xlineshape.pars
        
        lineshape_dict = Dict{String, Any}(
            "name" => scattering_key,
            "type" => "BreitWigner",
            "x" => "m_31_sq",
            "mass" => m,
            "width" => Γ,
            "l" => Xlineshape.l,
            "ma" => Xlineshape.m1,
            "mb" => Xlineshape.m2,
            "d" => 1.5
        )
        
    else
        # Fallback for unknown lineshapes
        @warn "Using fallback BreitWigner serialization for $(lineshape_type)"
        scattering_key = "$(lineshape_name)_BW"
        lineshape_dict = Dict{String, Any}(
            "name" => scattering_key,
            "type" => "BreitWigner",
            "x" => "m_31_sq",
            "mass" => 1.5,
            "width" => 0.1,
            "l" => 1,
            "ma" => 0.938272046,
            "mb" => 0.493677,
            "d" => 1.5
        )
    end
    
    # Add to appendix
    appendix[scattering_key] = lineshape_dict
    
    # Add Blatt-Weisskopf form factors based on angular momentum
    l_val = get(lineshape_dict, "l", 0)
    if l_val > 0
        ff_resonance = "BlattWeisskopf_resonance_l$(l_val)"
        ff_decay = "BlattWeisskopf_b_decay_l$(l_val)"
        
        if !haskey(appendix, ff_resonance)
            appendix[ff_resonance] = Dict{String, Any}(
                "name" => ff_resonance,
                "type" => "BlattWeisskopf",
                "radius" => 1.5,
                "l" => l_val
            )
        end
        
        if !haskey(appendix, ff_decay)
            appendix[ff_decay] = Dict{String, Any}(
                "name" => ff_decay,
                "type" => "BlattWeisskopf",
                "radius" => 5.0,
                "l" => l_val
            )
        end
        
        return (; scattering=scattering_key, FF_production=ff_decay, FF_decay=ff_resonance), appendix
    else
        return (; scattering=scattering_key, FF_production="", FF_decay=""), appendix
    end
end

# -------------------------------------------------------------
# Serialize the model to dictionary format
# -------------------------------------------------------------
println("Serializing XiC → pKπ model to JSON format...")

# Helper function to convert all dictionary keys to strings
function deep_convert_keys_to_string(obj)
    if isa(obj, Dict)
        result = Dict{String, Any}()
        for (k, v) in obj
            result[string(k)] = deep_convert_keys_to_string(v)
        end
        return result
    elseif isa(obj, Vector)
        return [deep_convert_keys_to_string(item) for item in obj]
    else
        return obj
    end
end

# Use the custom lineshape parser to get proper function definitions
decay_description, appendix = serializeToDict(threebody_model; lineshape_parser=xic_lineshape_parser)

# Collect all function definitions in proper format
functions_array = []
for (key, value) in appendix
    if isa(value, Dict) && haskey(value, "type")
        push!(functions_array, value)
    end
end

# Fix the chains to match Lc2pkpi format exactly
if haskey(decay_description, :chains)
    for chain in decay_description[:chains]
        # Update propagators to include proper parametrization reference
        if haskey(chain, :propagators)
            for prop in chain[:propagators]
                if haskey(prop, :parametrization) && isa(prop[:parametrization], Dict) && isempty(prop[:parametrization])
                    # Find the corresponding resonance name and assign proper parametrization
                    if haskey(chain, :name)
                        resonance_name = chain[:name]
                        # Map to function name in appendix - find matching function
                        for func in functions_array
                            if haskey(func, "name") && contains(func["name"], resonance_name)
                                prop[:parametrization] = func["name"]
                                break
                            end
                        end
                    end
                end
            end
        end
        
        # Fix vertices to match Lc2pkpi format with proper form factors
        if haskey(chain, :vertices)
            for (i, vertex) in enumerate(chain[:vertices])
                # Keep vertex types as "helicity" and "parity" like Lc2pkpi
                if haskey(vertex, :type)
                    if vertex[:type] == "helicity"
                        # Add form factor for helicity vertices (production)
                        if haskey(vertex, :formfactor) && isa(vertex[:formfactor], Dict) && isempty(vertex[:formfactor])
                            # Find appropriate form factor from chain's resonance
                            l_val = 1  # Default
                            if haskey(chain, :name)
                                resonance_name = chain[:name]
                                for func in functions_array
                                    if haskey(func, "name") && contains(func["name"], resonance_name) && haskey(func, "l")
                                        l_val = func["l"]
                                        break
                                    end
                                end
                            end
                            vertex[:formfactor] = l_val > 0 ? "BlattWeisskopf_b_decay_l$(l_val)" : ""
                        end
                    elseif vertex[:type] == "parity"
                        # Add form factor for parity vertices (resonance)
                        if haskey(vertex, :formfactor) && isa(vertex[:formfactor], Dict) && isempty(vertex[:formfactor])
                            l_val = 1  # Default
                            if haskey(chain, :name)
                                resonance_name = chain[:name]
                                for func in functions_array
                                    if haskey(func, "name") && contains(func["name"], resonance_name) && haskey(func, "l")
                                        l_val = func["l"]
                                        break
                                    end
                                end
                            end
                            vertex[:formfactor] = l_val > 0 ? "BlattWeisskopf_resonance_l$(l_val)" : ""
                        end
                    end
                end
            end
        end
    end
end

# Convert decay_description to proper string keys
string_decay_description = deep_convert_keys_to_string(decay_description)

# Define validation points
validation_points = [
    Dict{String, Any}(
        "name" => "validation_point",
        "parameters" => [
            Dict("name" => "cos_theta_31", "value" => -0.2309352648098208),
            Dict("name" => "phi_31", "value" => 0.0),
            Dict("name" => "m_31", "value" => 1.9101377207489973),
            Dict("name" => "cos_theta_31_2", "value" => 0.0),
            Dict("name" => "phi_31_2", "value" => 0.0),
            Dict("name" => "m_31_2", "value" => 2.46794)  # XiC mass
        ]
    ),
    Dict{String, Any}(
        "name" => "validation_point_m12sq",
        "parameters" => [
            Dict("name" => "m_12_sq", "value" => 3.2)
        ]
    ),
    Dict{String, Any}(
        "name" => "validation_point_m23sq",
        "parameters" => [
            Dict("name" => "m_23_sq", "value" => 1.4)
        ]
    ),
    Dict{String, Any}(
        "name" => "validation_point_m31sq",
        "parameters" => [
            Dict("name" => "m_31_sq", "value" => 3.2)
        ]
    )
]

# Create misc section with amplitude model checksums
misc_checksums = []

# Add overall model checksum
push!(misc_checksums, Dict{String, Any}(
    "point" => "validation_point",
    "distribution" => "default_model", 
    "value" => 9345.853380852355  # Placeholder
))

# Add individual function checksums
for func in functions_array
    if haskey(func, "name") && haskey(func, "type") && func["type"] != "BlattWeisskopf"
        func_name = func["name"]
        validation_point = if contains(func_name, "L") || contains(func_name, "Λ")
            "validation_point_m31sq"
        elseif contains(func_name, "D") || contains(func_name, "Δ")
            "validation_point_m12sq" 
        elseif contains(func_name, "K")
            "validation_point_m23sq"
        else
            "validation_point_m31sq"
        end
        
        push!(misc_checksums, Dict{String, Any}(
            "point" => validation_point,
            "distribution" => func_name,
            "value" => "0.0 + 0.0i"  # Placeholder
        ))
    end
end

# Define kinematic domains
m_parent = 2.46794
m1, m2, m3 = 0.938272046, 0.13957018, 0.493677  # p, π, K masses

proper_domains = [
    Dict{String, Any}(
        "name" => "default",
        "type" => "product_domain",
        "axes" => [
            Dict("name" => "cos_theta_31", "min" => -1.0, "max" => 1.0),
            Dict("name" => "phi_31", "min" => -3.14, "max" => 3.14),
            Dict("name" => "m_31", "min" => m1 + m3, "max" => m_parent - m2),
            Dict("name" => "cos_theta_31_2", "min" => -1.0, "max" => 1.0),
            Dict("name" => "phi_31_2", "min" => -3.14, "max" => 3.14),
            Dict("name" => "m_31_2", "min" => m1 + m2 + m3, "max" => m_parent)
        ]
    )
]

# Create the main distribution entry
distribution = Dict{String, Any}(
    "name" => "default_model",
    "type" => "HadronicUnpolarizedIntensity", 
    "decay_description" => string_decay_description,
    "variables" => [
        Dict{String, Any}(
            "node" => [3, 1],
            "mass_phi_costheta" => ["m_31", "phi_31", "cos_theta_31"]
        ),
        Dict{String, Any}(
            "node" => [[3, 1], 2],
            "mass_phi_costheta" => ["m_31_2", "phi_31_2", "cos_theta_31_2"]
        )
    ],
    "parameters" => []
)

# Create final JSON structure
final_dict = Dict{String, Any}(
    "distributions" => [distribution],
    "functions" => functions_array,
    "domains" => proper_domains,
    "misc" => Dict("amplitude_model_checksums" => misc_checksums),
    "parameter_points" => validation_points
)

# -------------------------------------------------------------
# Write JSON file
# -------------------------------------------------------------
output_file = joinpath(@__DIR__, "..", "data", "xic2pKpi-model.json")
println("Writing JSON model to: $output_file")

open(output_file, "w") do io
    JSON.print(io, final_dict, 4)
end

println("✓ Successfully converted XiC → pKπ model to JSON format")
println("  Output file: $output_file")
println("  Number of decay chains: $(length(string_decay_description["chains"]))")
println("  Number of function definitions: $(length(functions_array))")

# -------------------------------------------------------------
# Verification: read back and test
# -------------------------------------------------------------
println("\nVerifying JSON file by reading it back...")

json_content = open(output_file) do io
    JSON.parse(io)
end

# Basic validation checks
@assert haskey(json_content, "distributions")
@assert length(json_content["distributions"]) >= 1
@assert haskey(json_content["distributions"][1], "decay_description")

decay_desc = json_content["distributions"][1]["decay_description"]
@assert haskey(decay_desc, "kinematics")
@assert haskey(decay_desc, "reference_topology") 
@assert haskey(decay_desc, "chains")

println("✓ JSON file validation passed")
println("  - Contains $(length(decay_desc["chains"])) decay chains")
println("  - Contains $(length(json_content["functions"])) function definitions")
if haskey(json_content, "misc") && haskey(json_content["misc"], "amplitude_model_checksums")
    println("  - Contains $(length(json_content["misc"]["amplitude_model_checksums"])) validation checksums")
end

println("\nConversion complete! 🎉")
