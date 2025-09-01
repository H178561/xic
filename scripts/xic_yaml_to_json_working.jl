# -------------------------------------------------------------
# Working script to convert XiC to pKπ model from YAML to JSON format
# Fixed to handle all actual lineshape types
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
println("Loading YAML files...")
begin
    particledict = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-particle-definitions.yaml"))
    modelparameters = YAML.load_file(joinpath(@__DIR__, "..", "data", "xic-model-definitions.yaml"))
end

defaultparameters = modelparameters["Default amplitude model"]

# -------------------------------------------------------------
# Parse model dictionaries and convert to standard convention
# -------------------------------------------------------------
println("Parsing model dictionaries...")
(; chains, couplings, isobarnames) =
    parse_model_dictionaries(defaultparameters; particledict)

println("Found $(length(chains)) decay chains")
println("Resonances: $(join(isobarnames, ", "))")

# -------------------------------------------------------------
# Set up the amplitude model with particle numbering
# 0: Xic, 1:p, 2:pi, 3:K
# -------------------------------------------------------------
model = Lc2ppiKModel(; chains, couplings, isobarnames)

println("Model created successfully")
println("Model type: $(typeof(model))")
println("Number of decay chains: $(length(model.chains))")

# -------------------------------------------------------------
# Custom lineshape parser that handles all available lineshape types
# -------------------------------------------------------------
function xic_lineshape_parser(Xlineshape)
    appendix = Dict()
    
    # Get the lineshape type and name
    lineshape_type = typeof(Xlineshape)
    lineshape_name = Xlineshape.name
    
    println("Processing lineshape: $lineshape_name of type $lineshape_type")
    
    # Handle different lineshape types specifically
    if lineshape_type <: BreitWignerMinL
        scattering_key = "$(lineshape_name)_BW"
        m, Γ = Xlineshape.pars
        
        # Determine correct invariant mass variable based on resonance type
        x_var = if contains(lineshape_name, "L") || contains(lineshape_name, "Λ")
            "m_31_sq"  # Lambda resonances in pK system
        elseif contains(lineshape_name, "D") || contains(lineshape_name, "Δ")
            "m_12_sq"  # Delta resonances in p-pi system  
        elseif contains(lineshape_name, "K")
            "m_23_sq"  # Kaon resonances in pi-K system
        else
            "m_31_sq"  # Default to pK system
        end
        
        lineshape_dict = Dict{Symbol, Any}(
            :name => scattering_key,
            :type => "BreitWigner",
            :x => x_var,
            :mass => m,
            :width => Γ,
            :l => Xlineshape.l,
            :ma => Xlineshape.m1,
            :mb => Xlineshape.m2,
            :d => 1.5
        )
        
    elseif lineshape_type <: BuggBreitWignerMinL
        scattering_key = "$(lineshape_name)_BuggBW"
        m, Γ, γ = Xlineshape.pars
        
        # Create expression for Bugg lineshape
        threshold_sq = (Xlineshape.m1 + Xlineshape.m2)^2
        lineshape_dict = Dict{Symbol, Any}(
            :name => scattering_key,
            :type => "generic_function",
            :expression => "1/($(m)^2 - σ - i * $(m) * (σ - $(threshold_sq)) / ($(m)^2 - $(threshold_sq)) * $(Γ) * exp(-$(γ) * σ))"
        )
        
    elseif lineshape_type <: Flatte1405
        scattering_key = "$(lineshape_name)_Flatte"
        m, Γ = Xlineshape.pars
        
        lineshape_dict = Dict{Symbol, Any}(
            :name => scattering_key,
            :type => "MultichannelBreitWigner",
            :x => "m_31_sq",
            :mass => m,
            :channels => [
                Dict{Symbol, Any}(
                    :gsq => 0.23395150538434703,
                    :ma => Xlineshape.m1,
                    :mb => Xlineshape.m2,
                    :l => Xlineshape.l,
                    :d => 0
                ),
                Dict{Symbol, Any}(
                    :gsq => 0.23395150538434703,
                    :ma => 1.18937,  # Sigma mass
                    :mb => 0.13957018,  # Pion mass
                    :l => Xlineshape.l,
                    :d => 0
                )
            ]
        )
        
    elseif lineshape_type <: L1670Flatte
        scattering_key = "$(lineshape_name)_L1670Flatte"
        m, Γ = Xlineshape.pars
        
        # L1670 special case - treat as modified BreitWigner for now
        # TODO: Implement proper L1670Flatte in ThreeBodyDecaysIO
        lineshape_dict = Dict{Symbol, Any}(
            :name => scattering_key,
            :type => "BreitWigner",
            :x => "m_31_sq",
            :mass => m,
            :width => Γ,
            :l => Xlineshape.l,
            :ma => Xlineshape.m1,
            :mb => Xlineshape.m2,
            :d => 1.5
        )
        
    else
        @warn "Unknown lineshape type $(lineshape_type) for $(lineshape_name)"
        scattering_key = "$(lineshape_name)_Unknown"
        lineshape_dict = Dict{Symbol, Any}(
            :name => scattering_key,
            :type => "generic_function",
            :expression => "1.0"  # Placeholder
        )
    end
    
    # Add to appendix
    appendix[scattering_key] = lineshape_dict
    
    # Add Blatt-Weisskopf form factors based on angular momentum
    l_val = get(lineshape_dict, :l, 0)
    if l_val > 0
        ff_resonance = "BlattWeisskopf_resonance_l$(l_val)"
        ff_decay = "BlattWeisskopf_b_decay_l$(l_val)"
        
        if !haskey(appendix, ff_resonance)
            appendix[ff_resonance] = Dict{Symbol, Any}(
                :name => ff_resonance,
                :type => "BlattWeisskopf",
                :radius => 1.5,
                :l => l_val
            )
        end
        
        if !haskey(appendix, ff_decay)
            appendix[ff_decay] = Dict{Symbol, Any}(
                :name => ff_decay,
                :type => "BlattWeisskopf",
                :radius => 5.0,
                :l => l_val
            )
        end
        
        return (; scattering=scattering_key, FF_production=ff_decay, FF_decay=ff_resonance), appendix
    else
        return (; scattering=scattering_key, FF_production="", FF_decay=""), appendix
    end
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
output_file = joinpath(@__DIR__, "..", "data", "xic2pKpi-model_working.json")
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
    
    # Count different types of functions
    function_types = Dict{String, Int}()
    for func in dict[:functions]
        if haskey(func, "type")
            func_type = func["type"]
            function_types[func_type] = get(function_types, func_type, 0) + 1
        end
    end
    
    println("  Function types:")
    for (ftype, count) in function_types
        println("    $ftype: $count")
    end
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

# List all the resonances found
println("\nResonances found in conversion:")
resonance_names = Set{String}()
for func in json_content["functions"]
    if haskey(func, "name") && !contains(func["name"], "BlattWeisskopf")
        # Extract resonance name from function name
        func_name = func["name"]
        if contains(func_name, "_")
            resonance_name = split(func_name, "_")[1]
            push!(resonance_names, resonance_name)
        end
    end
end

println("  Unique resonances: $(length(resonance_names))")
for res_name in sort(collect(resonance_names))
    println("    - $res_name")
end

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
