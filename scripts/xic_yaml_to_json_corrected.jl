# -------------------------------------------------------------
# Script to convert XiC to pKπ model from YAML to JSON format
# Following the pattern from test_write_read_model.jl
# -------------------------------------------------------------
using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.JSON
using YAML

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
println("Number of decay chains: $(length(model.chains))")

# -------------------------------------------------------------
# Custom lineshape parser for XiC model following test pattern
# -------------------------------------------------------------
function xic_lineshape_parser(Xlineshape)
    appendix = Dict()
    
    # Get the lineshape name for the scattering key
    lineshape_name = string(typeof(Xlineshape).name.name)
    if hasfield(typeof(Xlineshape), :name)
        lineshape_name = Xlineshape.name
    end
    
    # Create scattering function name
    scattering_key = "$(lineshape_name)_func"
    
    # Serialize the lineshape to dictionary format
    serialized_lineshape, lineshape_appendix = serializeToDict(NamedArgFunc(Xlineshape, ["m12"]))
    
    # Update the scattering key in the serialized lineshape
    if haskey(lineshape_appendix, "m12")
        lineshape_appendix[scattering_key] = lineshape_appendix["m12"]
        delete!(lineshape_appendix, "m12")
        
        # Set correct invariant mass variable based on resonance type
        if haskey(lineshape_appendix[scattering_key], "x")
            x_var = if contains(lineshape_name, "L") || contains(lineshape_name, "Λ")
                "m_31_sq"  # Lambda resonances in pK system
            elseif contains(lineshape_name, "D") || contains(lineshape_name, "Δ")
                "m_12_sq"  # Delta resonances in p-pi system  
            elseif contains(lineshape_name, "K")
                "m_23_sq"  # Kaon resonances in pi-K system
            else
                "m_31_sq"  # Default to pK system
            end
            lineshape_appendix[scattering_key]["x"] = x_var
        end
    end
    
    merge!(appendix, lineshape_appendix)
    
    # Add form factors based on angular momentum
    l_val = 1  # Default
    if hasfield(typeof(Xlineshape), :l)
        l_val = Xlineshape.l
    end
    
    FF_decay = "BlattWeisskopf(resonance)"
    FF_production = "BlattWeisskopf(b-decay)"
    
    if l_val > 0
        ff_dict = Dict(
            "BlattWeisskopf(resonance)" =>
                Dict(:type => "BlattWeisskopf", :l => l_val, :radius => 1.5),
            "BlattWeisskopf(b-decay)" =>
                Dict(:type => "BlattWeisskopf", :l => l_val, :radius => 5.0),
        )
        merge!(appendix, ff_dict)
    else
        FF_decay = ""
        FF_production = ""
    end
    
    return (; scattering=scattering_key, FF_production, FF_decay), appendix
end

# -------------------------------------------------------------
# Serialize the model to dictionary format
# -------------------------------------------------------------
println("Serializing XiC → pKπ model to JSON format...")

# Use the custom lineshape parser following the test pattern
decay_description, appendix = serializeToDict(model; lineshape_parser=xic_lineshape_parser)

# Update appendix with correct variable names for each function
for (key, value) in appendix
    if isa(value, Dict) && haskey(value, :x) && value[:x] == "m12"
        # Determine correct variable based on function name
        if contains(key, "L") || contains(key, "Λ")
            value[:x] = "m_31_sq"  # Lambda resonances
        elseif contains(key, "D") || contains(key, "Δ")
            value[:x] = "m_12_sq"  # Delta resonances
        elseif contains(key, "K")
            value[:x] = "m_23_sq"  # Kaon resonances
        else
            value[:x] = "m_31_sq"  # Default
        end
    end
end

# Add HS3 fields to create complete JSON structure
dict = add_hs3_fields(decay_description, appendix, "xic_default_model")

# -------------------------------------------------------------
# Write JSON file
# -------------------------------------------------------------
output_file = joinpath(@__DIR__, "..", "data", "xic2pKpi-model.json")
println("Writing JSON model to: $output_file")

open(output_file, "w") do io
    JSON.print(io, dict, 4)
end

println("✓ Successfully converted XiC → pKπ model to JSON format")
println("  Output file: $output_file")
println("  Number of decay chains: $(length(decay_description[:chains]))")
println("  Number of functions: $(length(dict["functions"]))")

# -------------------------------------------------------------
# Verification: read back and test like in test_write_read_model.jl
# -------------------------------------------------------------
println("\nVerifying JSON file by reading it back...")

json_content = open(output_file) do io
    JSON.parse(io)
end

# Basic validation checks from test pattern
decay_desc = json_content["distributions"][1]["decay_description"]
@assert haskey(decay_desc, "kinematics")
@assert haskey(decay_desc, "reference_topology") 
@assert haskey(decay_desc, "chains")

chains = decay_desc["chains"]
@assert length(chains) >= 1

for chain in chains
    @assert haskey(chain, "vertices")
    @assert haskey(chain, "propagators")
    @assert length(chain["vertices"]) == 2
    @assert length(chain["propagators"]) == 1
end

println("✓ JSON file validation passed")
println("  - Contains $(length(chains)) decay chains")
println("  - Contains $(length(json_content["functions"])) function definitions")

# Test parsing back to model (following test pattern)
println("\nTesting model reconstruction...")

input = copy(json_content)
decay_description_read = input["distributions"][1]["decay_description"]
functions = input["functions"]

workspace = Dict{String,Any}()
for fn in functions
    name = fn["name"]
    type_str = fn["type"]
    instance_type = eval(Symbol(type_str))
    workspace[name] = dict2instance(instance_type, fn)
end

reconstructed_model = dict2instance(ThreeBodyDecay, decay_description_read; workspace)
@assert reconstructed_model isa ThreeBodyDecay

println("✓ Model reconstruction successful")
println("\nConversion complete! 🎉")
