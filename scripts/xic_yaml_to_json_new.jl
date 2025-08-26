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
function xic_lineshape_parser(Xlineshape)
    appendix = Dict()
    
    # Get the lineshape name
    lineshape_name = string(typeof(Xlineshape).name.name)
    if hasfield(typeof(Xlineshape), :name)
        lineshape_name = Xlineshape.name
    end
    
    # Create scattering function name
    scattering_key = "$(lineshape_name)_BW"
    
    # For custom lineshapes that don't have serializeToDict methods,
    # we need to manually create the function definition
    try
        # Try to serialize the actual lineshape (works for standard ones)
        scattering_dict = serializeToDict(NamedArgFunc(Xlineshape, ["m12"]))[1]
        appendix[scattering_key] = scattering_dict
    catch
        # Fall back to manual creation for custom lineshapes
        lineshape_dict = Dict{Symbol, Any}(
            :type => "BreitWigner",  # Default type
            :x => "m12"
        )
        
        # Extract parameters if available
        if hasfield(typeof(Xlineshape), :pars)
            pars = Xlineshape.pars
            if length(pars) >= 2
                lineshape_dict[:mass] = pars[1]
                lineshape_dict[:width] = pars[2]
            end
        end
        
        # Add mass info if available
        if hasfield(typeof(Xlineshape), :m1) && hasfield(typeof(Xlineshape), :m2)
            lineshape_dict[:ma] = Xlineshape.m1
            lineshape_dict[:mb] = Xlineshape.m2
        end
        
        # Add angular momentum if available
        if hasfield(typeof(Xlineshape), :l)
            lineshape_dict[:l] = Xlineshape.l
        else
            lineshape_dict[:l] = 1  # Default
        end
        
        lineshape_dict[:d] = 1.5  # Default radius
        
        appendix[scattering_key] = lineshape_dict
    end
    
    # Set correct invariant mass variable
    x_var = if contains(lineshape_name, "L") || contains(lineshape_name, "Λ")
        "m_31_sq"  # Lambda resonances in pK system
    elseif contains(lineshape_name, "D") || contains(lineshape_name, "Δ")
        "m_12_sq"  # Delta resonances in p-pi system  
    elseif contains(lineshape_name, "K")
        "m_23_sq"  # Kaon resonances in pi-K system
    else
        "sigma"  # Use generic sigma like in test
    end
    
    # Update the x variable in the function definition
    if haskey(appendix[scattering_key], :x)
        appendix[scattering_key][:x] = x_var
    end
    
    # Add form factors like in test_write_read_model.jl
    FF_decay = "BlattWeisskopf(resonance)"
    FF_production = "BlattWeisskopf(b-decay)"
    ff_dict = Dict(
        "BlattWeisskopf(resonance)" =>
            Dict(:type => "BlattWeisskopf", :l => 1, :radius => 1.5),
        "BlattWeisskopf(b-decay)" =>
            Dict(:type => "BlattWeisskopf", :l => 1, :radius => 5.0),
    )
    merge!(appendix, ff_dict)
    
    return (; scattering=scattering_key, FF_production, FF_decay), appendix
end

# -------------------------------------------------------------
# Serialize model following the test pattern exactly
# -------------------------------------------------------------
println("Serializing XiC → pKπ model to JSON format...")

# Use our custom lineshape parser that handles custom lineshapes
decay_description, appendix = serializeToDict(model; lineshape_parser=xic_lineshape_parser)

# Update all x variables to use sigma like in test
for (key, value) in appendix
    if isa(value, Dict) && haskey(value, :x)
        value[:x] = "sigma"
    end
end

# Create the JSON structure using add_hs3_fields like in test
dict = add_hs3_fields(decay_description, appendix, "xic_default_model")

# -------------------------------------------------------------
# Write JSON file
# -------------------------------------------------------------
output_file = joinpath(@__DIR__, "..", "data", "xic2pKpi-model_new.json")
println("Writing JSON model to: $output_file")

open(output_file, "w") do io
    JSON.print(io, dict, 4)
end

println("✓ Successfully converted XiC → pKπ model to JSON format")
println("  Output file: $output_file")
if haskey(decay_description, :chains)
    println("  Number of decay chains: $(length(decay_description[:chains]))")
end
println("  Available keys in dict: $(keys(dict))")
if haskey(dict, :functions)
    println("  Number of functions: $(length(dict[:functions]))")
else
    println("  No 'functions' key found")
end

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
