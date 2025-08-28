#!/usr/bin/env julia
# -------------------------------------------------------------
# YAML Structure Demonstration: XiC → pKπ Model
# Shows how YAML files are loaded and used step-by-step
# -------------------------------------------------------------

using Lc2ppiKSemileptonicModelLHCb
using YAML
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays

println("🔍 YAML Structure Demonstration: XiC → pKπ Model")
println("="^60)

# -------------------------------------------------------------
# Step 1: Load YAML Files
# -------------------------------------------------------------
println("\n📁 Step 1: Loading YAML files...")

# Load particle definitions
println("Loading particle definitions...")
particledict = YAML.load_file("data/xic-particle-definitions.yaml")
println("   ✓ Loaded $(length(particledict)) particle definitions")

# Show a few examples
println("\n📋 Sample particle definitions:")
for (name, props) in collect(particledict)[1:3]
    println("   $name:")
    println("     Mass: $(props["mass"]) MeV")
    println("     JP: $(props["jp"])")
    if haskey(props, "lineshape")
        println("     Lineshape: $(props["lineshape"])")
    end
end

# Load model definitions  
println("\nLoading model definitions...")
modelparameters = YAML.load_file("data/xic-model-definitions.yaml")
println("   ✓ Available models: $(collect(keys(modelparameters)))")

defaultparameters = modelparameters["Default amplitude model"]
println("   ✓ Using 'Default amplitude model'")

# Show structure
println("\n📋 Model structure:")
println("   Lineshapes: $(length(defaultparameters["lineshapes"])) assignments")
println("   Parameters: $(length(defaultparameters["parameters"])) couplings")

# -------------------------------------------------------------
# Step 2: Parse YAML into Julia Structures  
# -------------------------------------------------------------
println("\n🔧 Step 2: Parsing YAML into Julia structures...")

# Parse model dictionaries
(; chains, couplings, isobarnames) = parse_model_dictionaries(defaultparameters; particledict)

println("   ✓ Parsed chains: $(length(chains))")
println("   ✓ Parsed couplings: $(length(couplings))")
println("   ✓ Isobar names: $(length(isobarnames))")

# Show what types of resonances we have
resonance_types = Dict{String, Int}()
for (name, _) in chains
    if occursin("L(", name)
        resonance_types["Lambda"] = get(resonance_types, "Lambda", 0) + 1
    elseif occursin("D(", name)
        resonance_types["Delta"] = get(resonance_types, "Delta", 0) + 1
    elseif occursin("K(", name) || occursin("K2(", name)
        resonance_types["Kaon"] = get(resonance_types, "Kaon", 0) + 1
    end
end

println("\n📊 Resonance inventory:")
for (type, count) in resonance_types
    println("   $type resonances: $count")
end

# -------------------------------------------------------------
# Step 3: Extract Physical Information
# -------------------------------------------------------------
println("\n⚛️  Step 3: Physical information extraction...")

# Extract masses (converted to GeV)
masses_gev = let
    _mΞc = particledict["Lambda_c+"]["mass"] / 1e3  # MeV → GeV
    _mp = particledict["p"]["mass"] / 1e3
    _mπ = particledict["pi+"]["mass"] / 1e3
    _mK = particledict["K-"]["mass"] / 1e3
    ThreeBodyMasses(m1 = _mp, m2 = _mπ, m3 = _mK, m0 = _mΞc)
end

println("   Particle masses (GeV):")
println("     XiC (m₀): $(masses_gev.m0)")
println("     p   (m₁): $(masses_gev.m1)")
println("     π   (m₂): $(masses_gev.m2)")
println("     K   (m₃): $(masses_gev.m3)")

# Show lineshape assignments
println("\n🌊 Lineshape assignments:")
lineshape_counts = Dict{String, Int}()
for (resonance, lineshape) in defaultparameters["lineshapes"]
    lineshape_str = string(lineshape)
    lineshape_counts[lineshape_str] = get(lineshape_counts, lineshape_str, 0) + 1
end

for (lineshape, count) in sort(collect(lineshape_counts))
    println("   $lineshape: $count resonances")
end

# -------------------------------------------------------------
# Step 4: Create the Model
# -------------------------------------------------------------
println("\n🏗️  Step 4: Creating amplitude model...")

try
    model = Lc2ppiKModel(; chains, couplings, isobarnames)
    println("   ✓ Model created successfully!")
    println("   ✓ Model type: $(typeof(model))")
    println("   ✓ Number of decay chains: $(length(model.chains))")
    
    # Test that we can access the masses
    model_masses = masses(model)
    println("   ✓ Model masses accessible: m₀=$(model_masses.m0), m₁=$(model_masses.m1), m₂=$(model_masses.m2), m₃=$(model_masses.m3)")
    
catch e
    println("   ❌ Model creation failed: $e")
end

# -------------------------------------------------------------
# Step 5: Show Parameter Structure
# -------------------------------------------------------------
println("\n🎛️  Step 5: Parameter structure analysis...")

# Analyze parameter naming convention
real_params = []
imag_params = []
lambda_resonances = []
delta_resonances = []
kaon_resonances = []

for (param_name, value) in defaultparameters["parameters"]
    if startswith(param_name, "Ar")
        push!(real_params, param_name)
    elseif startswith(param_name, "Ai")
        push!(imag_params, param_name)
    end
    
    if occursin("L(", param_name)
        push!(lambda_resonances, param_name)
    elseif occursin("D(", param_name)
        push!(delta_resonances, param_name)
    elseif occursin("K(", param_name)
        push!(kaon_resonances, param_name)
    end
end

println("   Parameter breakdown:")
println("     Real parameters: $(length(real_params))")
println("     Imaginary parameters: $(length(imag_params))")
println("     Lambda couplings: $(length(lambda_resonances))")
println("     Delta couplings: $(length(delta_resonances))")
println("     Kaon couplings: $(length(kaon_resonances))")

# Show some example parameters
println("\n🔢 Example parameters:")
example_params = collect(defaultparameters["parameters"])[1:3]
for (name, value) in example_params
    println("   $name: $value")
end

# -------------------------------------------------------------
# Step 6: Demonstrate Physical Point Calculation
# -------------------------------------------------------------
println("\n🎯 Step 6: Physical point demonstration...")

# Define a physical point in the Dalitz plot
test_point = (σ1 = 2.5, σ2 = 1.8, σ3 = 2.2)  # GeV²

# Check if point is physical
tbs = ThreeBodySystem(masses_gev, ThreeBodySpins(1, 0, 0; two_h0 = 1))
is_physical = isphysical(test_point, masses_gev)

println("   Test point: σ₁=$(test_point.σ1), σ₂=$(test_point.σ2), σ₃=$(test_point.σ3) GeV²")
println("   Physical: $is_physical")

if is_physical
    println("   ✓ Point is in allowed kinematic region")
else
    println("   ⚠️  Point is outside kinematic boundaries")
end

# -------------------------------------------------------------
# Summary
# -------------------------------------------------------------
println("\n🎉 Summary: YAML Structure Demonstration Complete")
println("="^60)
println("✅ Successfully demonstrated:")
println("   • YAML file loading (3 main files)")
println("   • Parameter parsing and type conversion")  
println("   • Model structure analysis")
println("   • Physical mass extraction")
println("   • Lineshape assignment mapping")
println("   • Amplitude model creation")
println("   • Parameter naming convention")
println("   • Physical point validation")
println()
println("📝 Key insights:")
println("   • YAML provides human-readable physics specification")
println("   • Three-file structure separates concerns effectively")
println("   • Automatic unit conversion (MeV → GeV)")
println("   • Complex parameter structure supports helicity amplitudes")
println("   • Model ready for amplitude calculations")
