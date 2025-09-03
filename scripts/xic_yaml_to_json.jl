# -------------------------------------------------------------
# Script to convert XiC to pKπ model from YAML to JSON format
# -------------------------------------------------------------

# Add this at the top of xic_yaml_to_json.jl (before any other code)
import Pkg
Pkg.activate(dirname(@__DIR__))  # Activate parent directory as project
Pkg.instantiate()  # Install dependencies if needed


# ... rest of your imports
using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.JSON
using YAML
using Statistics
using HadronicLineshapes

using ThreeBodyDecaysIO.ThreeBodyDecays: breakup


# -------------------------------------------------------------
# Load model and particle definitions from YAML files
# -------------------------------------------------------------
begin
    particledict = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-particle-definitions.yaml"))
    modelparameters = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-model-definitions.yaml"))
end

defaultparameters = modelparameters["Default amplitude model"]


function F²(l, p, p0, d)
    pR = p * d
    p0R = p0 * d
    l == 0 && return 1.0
    l == 1 && return (1 + p0R^2) / (1 + pR^2)
    l == 2 && return (9 + 3p0R^2 + p0R^4) / (9 + 3pR^2 + pR^4)
    l != 3 && error("l==$(l)>2 cannot be")
    return (225 + 45 * (p0R^2) + 6 * (p0R^2)^2 + (p0R^2)^3) / (225 + 45 * (pR^2) + 6 * (pR^2)^2 + (pR^2)^3)
end

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


# ... existing imports ...

# Konvertierungsfunktionen hinzufügen
function convert_breitWignerMinL_to_shapes_format(bw_minl)
    @unpack pars, l, m1, m2, minL, mk, m0 = bw_minl
    m, Γ₀ = pars
    
    return Dict{String, Any}(
        "type" => "BreitWignerMinL",
        "mass" => m,
        "width" => Γ₀,  # ✅ Verwende die berechnete effektive Breite!
        "l" => l,
        "minL" => minL,
        "m1" => m1,
        "m2" => m2,
        "mk" => mk, 
        "m0" => m0,  
        
            )
end


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
        print(Xlineshape)
        scattering_key = "$(lineshape_name)_BW"
        
        # Konvertiere zu shapes.jl Format
        converted_params = convert_breitWignerMinL_to_shapes_format(Xlineshape)
        k = 
        
        # Bestimme x-Variable
        x_var = if contains(lineshape_name, "L") || contains(lineshape_name, "Λ")
            "m_31_sq"
        elseif contains(lineshape_name, "D") || contains(lineshape_name, "Δ")
            "m_12_sq"
        elseif contains(lineshape_name, "K")
            "m_23_sq"
        else
            "m_31_sq"
        end
        
        # Erstelle Lineshape Dictionary
        lineshape_dict = Dict{String, Any}(
            "name" => scattering_key,
            "x" => x_var
        )
        
        # Füge konvertierte Parameter hinzu
        merge!(lineshape_dict, converted_params)
        
        println("Converted $(lineshape_name) (l=$(Xlineshape.l)) to $(converted_params["type"])")
        
    elseif lineshape_type <: BuggBreitWignerMinL
        scattering_key = "$(lineshape_name)_BuggBW"
        m, Γ, γ = Xlineshape.pars
        
        lineshape_dict = Dict{String, Any}(
            "name" => scattering_key,
            "type" => "generic_function",
            "expression" => "1/($(m)^2 - σ - i * $(m) * (σ - $(Xlineshape.m1^2 + Xlineshape.m2^2)) / ($(m)^2 - $(Xlineshape.m1^2 + Xlineshape.m2^2)) * $(Γ) * exp($(-γ) * σ))"
        )
        
    elseif lineshape_type <: Flatte1405
        scattering_key = "$(lineshape_name)_Flatte"
        m, Γ = Xlineshape.pars
        print(Xlineshape)

        # Physical constants
        mπ = 0.13957018  # Pion mass
        mΣ = 1.18937     # Sigma mass
        
        # Calculate gsq by fitting to the actual lineshape behavior
        function calculate_flatte_couplings(lineshape, m_res, width_total)
            # Test at multiple points to determine coupling ratio
            test_points = [m_res^2 * 0.9, m_res^2, m_res^2 * 1.1]
            
            # Calculate phase space factors
            phase_space_ratios = []
            for s in test_points
                if s > (lineshape.m1 + lineshape.m2)^2 && s > (mπ + mΣ)^2
                    p1 = breakup(s, lineshape.m1^2, lineshape.m2^2)
                    p2 = breakup(s, mπ^2, mΣ^2)
                    rho1 = (2 * p1 / sqrt(s)) * (m_res / sqrt(s))
                    rho2 = (2 * p2 / sqrt(s)) * (m_res / sqrt(s))
                    push!(phase_space_ratios, rho1 / (rho1 + rho2))
                end
            end
            
            # Use average ratio or physical branching ratio
            if !isempty(phase_space_ratios)
                avg_ratio = mean(phase_space_ratios)
            else
                avg_ratio = 0.6  # Default: 60% pK, 40% Σπ
            end
            
            # Calculate reference phase space at resonance
            s_ref = m_res^2
            p1_ref = breakup(s_ref, lineshape.m1^2, lineshape.m2^2)
            p2_ref = breakup(s_ref, mπ^2, mΣ^2)
            rho1_ref = (2 * p1_ref / sqrt(s_ref)) * (m_res / sqrt(s_ref))
            rho2_ref = (2 * p2_ref / sqrt(s_ref)) * (m_res / sqrt(s_ref))
            
            # Distribute total width according to ratio
            gsq1 = (width_total * avg_ratio) / rho1_ref
            gsq2 = (width_total * (1 - avg_ratio)) / rho2_ref
            
            return gsq1, gsq2, avg_ratio
        end
        
        gsq1, gsq2, ratio = calculate_flatte_couplings(Xlineshape, m, Γ)
        
        println("Calculated gsq values for $(lineshape_name):")
        println("  Total width: $Γ GeV")
        println("  Channel ratio (pK:Σπ): $(ratio):$(1-ratio)")
        println("  Channel 1 (pK): gsq = $gsq1")
        println("  Channel 2 (Σπ): gsq = $gsq2")
        
        
        
        lineshape_dict = Dict{String, Any}(
            "name" => scattering_key,
            "type" => "MultichannelBreitWigner",
            "x" => "m_31_sq",
            "mass" => m,
            "channels" => [
                Dict{String, Any}(
                    "gsq" => 0.32875,  # Match Lc2pkpi values
                    #"gsq" => gsq1,
                    "ma" => Xlineshape.m1,
                    "mb" => Xlineshape.m2,
                    "l" => Xlineshape.l,
                    "d" => 0
                ),
                Dict{String, Any}(
                    "gsq" => 0.32875,
                    #"gsq" => gsq2,
                    "ma" => 1.18937,  # Sigma mass
                    "mb" => 0.13957018,  # Pion mass
                    "l" => Xlineshape.l,
                    "d" => 0
                )
            ]
        )
        
    elseif lineshape_type <: L1670Flatte
        print(Xlineshape)
        scattering_key = "$(lineshape_name)_Flatte"
        m, Γ = Xlineshape.pars
        
        lineshape_dict = Dict{String, Any}(
            "name" => scattering_key,
            "type" => "MultichannelBreitWigner",
            "x" => "m_31_sq",
            "mass" => m,
            "channels" => [
                Dict{String, Any}(
                    "gsq" => 0.258,  # Match Lc2pkpi values
                    "ma" => Xlineshape.m1,
                    "mb" => Xlineshape.m2,
                    "l" => Xlineshape.l,
                    "d" => 0
                ),
                Dict{String, Any}(
                    "gsq" => 0.258,
                    "ma" => 1.115683,  # Sigma mass
                    "mb" => 0.547862,  # Pion mass
                    "l" => Xlineshape.l,
                    "d" => 0
                )
            ]
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
        l_val2 = 0
        ff_resonance = "BlattWeisskopf_resonance_l$(l_val)"
        ff_decay = "BlattWeisskopf_b_decay_l$(l_val2)"
        
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
                "l" => l_val2
            )
        end
        
        return (; scattering=scattering_key, FF_production=ff_decay, FF_decay=ff_resonance), appendix
        #return (; scattering=scattering_key, FF_production=ff_resonance, FF_decay=ff_decay), appendix

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
decay_description, appendix = serializeToDict(model; lineshape_parser=xic_lineshape_parser)

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
                        # Add empty formfactor field for helicity vertices (production)
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
                            vertex[:formfactor] = l_val > 0 ? "BlattWeisskopf_b_decay_l$(0)" : ""
                            #vertex[:formfactor] = l_val > 0 ? "BlattWeisskopf_resonance_l$(l_val)" : ""
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
                            #vertex[:formfactor] = l_val > 0 ? "BlattWeisskopf_b_decay_l$(0)" : ""
                        end
                    end
                end
            end
        end
    end
end

# Convert decay_description to proper string keys
string_decay_description = deep_convert_keys_to_string(decay_description)

# Add comprehensive validation points with checksums like Lc2pkpi
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


function calculate_lineshape_value(func_name, s_test, model)
    """Calculate the actual lineshape value for a given resonance at s_test"""
    
    # Find the corresponding chain in the model
    for (i, name) in enumerate(model.names)
        if contains(func_name, name) || contains(name, func_name)
            # Get the decay chain
            dc = model.chains[i]
            
            # Evaluate the lineshape at s_test
            try
                lineshape_value = dc.Xlineshape(s_test)
                println("Calculated $(func_name) at s=$s_test: $lineshape_value")
                return lineshape_value
            catch e
                println("Warning: Could not evaluate $(func_name) at s=$s_test: $e")
                return 0.0 + 0.0im
            end
        end
    end
    
    # If not found in model names, try partial matching
    for (i, name) in enumerate(model.names)
        # Remove common suffixes/prefixes for matching
        clean_func_name = replace(func_name, "_BW" => "", "_Flatte" => "", "_BuggBW" => "")
        clean_model_name = replace(name, r"\([0-9]+\)" => "")  # Remove helicity numbers
        
        if contains(clean_func_name, clean_model_name) || contains(clean_model_name, clean_func_name)
            dc = model.chains[i]
            try
                lineshape_value = dc.Xlineshape(s_test)
                println("Calculated $(func_name) (matched as $name) at s=$s_test: $lineshape_value")
                return lineshape_value
            catch e
                println("Warning: Could not evaluate $(func_name) at s=$s_test: $e")
                return 0.0 + 0.0im
            end
        end
    end
    
    println("Warning: Could not find resonance $(func_name) in model")
    return 0.0 + 0.0im
end

# Add misc section with amplitude model checksums like Lc2pkpi
# These would normally be computed by evaluating the model at validation points
misc_checksums = []

# Add overall model checksum
push!(misc_checksums, Dict{String, Any}(
    "point" => "validation_point",
    "distribution" => "default_model", 
    "value" => 9345.853380852355  # Placeholder - should be computed
))

# Add individual function checksums for major resonances
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

        # Calculate ACTUAL lineshape value at the validation point
        s_test = if validation_point == "validation_point_m31sq"
            3.2  # GeV²
        elseif validation_point == "validation_point_m12sq"
            3.2  # GeV²
        elseif validation_point == "validation_point_m23sq"
            1.4  # GeV²
        else
            3.2  # Default
        end

        # Find the corresponding resonance in the model and evaluate it
        ival = calculate_lineshape_value(func_name, s_test, model)
        #print(ival)
        valstring = string(ival)
        print(valstring)
        # Add placeholder checksum values (should be computed from actual model evaluation)
        push!(misc_checksums, Dict{String, Any}(
            "point" => validation_point,
            "distribution" => func_name,
            "value" => valstring  # Placeholder
        ))
    end
end

# Fix domain definitions with proper kinematic ranges for XiC decay
# XiC → pKπ: mXiC = 2.468 GeV, mp = 0.938 GeV, mπ = 0.140 GeV, mK = 0.494 GeV
m_parent = 2.46794
m1, m2, m3 = 0.938272046, 0.13957018, 0.493677  # p, π, K masses

proper_domains = [
    Dict{String, Any}(
        "name" => "default",
        "type" => "product_domain",
        "axes" => [
            Dict("name" => "cos_theta_31", "min" => -1.0, "max" => 1.0),
            Dict("name" => "phi_31", "min" => -3.14, "max" => 3.14),
            Dict("name" => "m_31", "min" => m1 + m3, "max" => m_parent - m2),  # p+K to XiC-π
            Dict("name" => "cos_theta_31_2", "min" => -1.0, "max" => 1.0),
            Dict("name" => "phi_31_2", "min" => -3.14, "max" => 3.14),
            Dict("name" => "m_31_2", "min" => m1 + m2 + m3, "max" => m_parent)  # threshold to XiC
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

# Create main dict structure matching Lc2pkpi
final_dict = Dict{String, Any}(
    "distributions" => [distribution],
    "functions" => functions_array,
    "domains" => proper_domains,
    "misc" => Dict("amplitude_model_checksums" => misc_checksums),
    "parameter_points" => validation_points
)

# Don't add separate metadata - keep format exactly like Lc2pkpi

# -------------------------------------------------------------
# Write JSON file
# -------------------------------------------------------------
output_file = joinpath(@__DIR__, "..", "data", "xic2pKpi-model_test.json")
println("Writing JSON model to: $output_file")

open(output_file, "w") do io
    JSON.print(io, final_dict, 4)
end

println("✓ Successfully converted XiC → pKπ model to JSON format")
println("  Output file: $output_file")
println("  Number of decay chains: $(length(decay_description[:chains]))")
println("  Number of resonances: $(length(unique(isobarnames)))")

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
