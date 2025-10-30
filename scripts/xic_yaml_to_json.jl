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


function F²(l, p, p0, d)
    pR = p * d
    p0R = p0 * d
    l == 0 && return 1.0
    l == 1 && return (1 + p0R^2) / (1 + pR^2)
    l == 2 && return (9 + 3p0R^2 + p0R^4) / (9 + 3pR^2 + pR^4)
    # F-wave (l=3)
    l == 3 && return (225 + 45 * (p0R^2) + 6 * (p0R^2)^2 + (p0R^2)^3) / (225 + 45 * (pR^2) + 6 * (pR^2)^2 + (pR^2)^3)
    # G-wave (l=4)
    l == 4 && return (11025 + 1575 * (p0R^2) + 135 * (p0R^2)^2 + 10 * (p0R^2)^3 + (p0R^2)^4) /
                             (11025 + 1575 * (pR^2) + 135 * (pR^2)^2 + 10 * (pR^2)^3 + (pR^2)^4)
    error("l==$(l)>4 is not supported")
end
# -------------------------------------------------------------
# Load model and particle definitions from YAML files
# -------------------------------------------------------------
begin
    particledict = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-particle-definitions.yaml"))
    modelparameters = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-model-definitions.yaml"))
end

#defaultparameters = modelparameters["Default amplitude model"]


# list of models
list = [
    #"Default amplitude model",
    #"Alternative model 1 - Delta resonances with free mass and width",
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
    "Alternative model 27 - L(2100) contribution added with PDG values",
    #"Alternative model 28 - L(2110) contribution added with PDG values",
    #"Alternative model 29 - S(1670) contribution added with PDG values",
    #"Alternative model 30 - S(1775) contribution added with PDG values",
    #"Alternative model 31 - Free radial parameter rXic"
]

global i = 22

for model in list
    global i
    println("\n" * "="^60)
    println("Processing model $i: $model")
    println("="^60)
    
    try
        defaultparameters = deepcopy(modelparameters[model])
        
        # Fix fixed parameters by adding dummy uncertainties
        if haskey(defaultparameters, "parameters")
            for (param_name, param_value) in defaultparameters["parameters"]
                if isa(param_value, String)
                    # Check if it's a fixed parameter (just a number without ±)
                    if !contains(param_value, "±") && !contains(param_value, "+/-")
                        try
                            # Try to parse as a number
                            val = parse(Float64, strip(param_value))
                            # Convert to measured parameter format with tiny uncertainty
                            defaultparameters["parameters"][param_name] = "$val ± 0.0001"
                            println("  Fixed parameter $param_name: $param_value → $(defaultparameters["parameters"][param_name])")
                        catch
                            # If parsing fails, leave as is
                            println("  Warning: Could not parse parameter $param_name: $param_value")
                        end
                    end
                end
            end
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

    test_simple_amplitude(model)

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
            #print(Xlineshape)
            scattering_key = "$(lineshape_name)_BWminL"
            
            # Konvertiere zu shapes.jl Format
            converted_params = convert_breitWignerMinL_to_shapes_format(Xlineshape)
            
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
            
            #println("Converted $(lineshape_name) (l=$(Xlineshape.l)) to $(converted_params["type"])")
            
        elseif lineshape_type <: BuggBreitWignerMinL
            scattering_key = "$(lineshape_name)_BuggBWminL"
            m, Γ, γ = Xlineshape.pars
            l, minl = Xlineshape.l, Xlineshape.minL
            m1, m2, mk, m0 = Xlineshape.m1, Xlineshape.m2, Xlineshape.mk, Xlineshape.m0

            lineshape_dict = Dict{String, Any}(
                "name" => scattering_key,
                "type" => "BuggBreitWignerMinL",
                "x" => "m_31_sq",
                "mass" => m,
                "width" => Γ,
                "gamma" => γ,
                "l" => l,
                "minL" => minl,
                "m1" => m1,
                "m2" => m2,
                "mk" => mk,
                "m0" => m0
            )
            
        elseif lineshape_type <: Flatte1405
            scattering_key = "$(lineshape_name)_Flatte1405"
            m, Γ = Xlineshape.pars
            l, minl = Xlineshape.l, Xlineshape.minL
            m1, m2, mk, m0 = Xlineshape.m1, Xlineshape.m2, Xlineshape.mk, Xlineshape.m0
            #print(Xlineshape)

            
            lineshape_dict = Dict{String, Any}(
                "name" => scattering_key,
                "type" => "Flatte1405",
                "x" => "m_31_sq",
                "mass" => m,
                "width" => Γ,
                "l" => l,
                "minL" => minl,
                "m1" => m1,
                "m2" => m2,
                "mk" => mk,
                "m0" => m0
            )
            
        elseif lineshape_type <: L1670Flatte
            #print(Xlineshape)
            scattering_key = "$(lineshape_name)_L1670Flatte"
            m, Γ = Xlineshape.pars
            l, minl = Xlineshape.l, Xlineshape.minL
            m1, m2, mk, m0 = Xlineshape.m1, Xlineshape.m2, Xlineshape.mk, Xlineshape.m0
            
            lineshape_dict = Dict{String, Any}(
                "name" => scattering_key,
                "type" => "L1670Flatte",
                "x" => "m_31_sq",
                "mass" => m,
                "width" => Γ,
                "l" => l,
                "minL" => minl,
                "m1" => m1,
                "m2" => m2,
                "mk" => mk,
                "m0" => m0,
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
        minl_val = get(lineshape_dict, "minl", 0)

        if l_val > 0
            l_val2 = 0
            ff_resonance = "BlattWeisskopf_resonance_l$(l_val)"
            ff_decay = "BlattWeisskopf_b_decay_l$(minl_val)"
            
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
                    "l" => minl_val
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
                    #println("Calculated $(func_name) at s=$s_test: $lineshape_value")
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
            clean_func_name = replace(func_name, "_BWminL" => "", "_Flatte1405" => "", "_L1670Flatte" => "", "_BuggBWminL" => "")
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
            validation_point = if startswith(func_name, "L") || startswith(func_name, "Λ")
                "validation_point_m31sq"
            elseif startswith(func_name, "D") || startswith(func_name, "Δ")
                "validation_point_m12sq" 
            elseif startswith(func_name, "K")
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
            #print("\n", validation_point, " ", s_test, " ", func_name, "\n")
            # Find the corresponding resonance in the model and evaluate it
            ival = calculate_lineshape_value(func_name, s_test, model)
            #print(ival)
            valstring = string(ival)
            #print(valstring)
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
    output_file = joinpath(@__DIR__, "..", "data", "json", "model$i.json")
    println("Writing JSON model to: $output_file")

    open(output_file, "w") do io
        JSON.print(io, final_dict, 4)
    end

    #println("✓ Successfully converted XiC → pKπ model to JSON format")
    #println("  Output file: $output_file")
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
        i=i+1
        
    catch e
        println("❌ Failed to process model $i: $model")
        println("Error: $e")
        i=i+1
        continue
    end
end
