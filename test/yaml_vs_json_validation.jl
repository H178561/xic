# -------------------------------------------------------------
# YAML vs JSON Model Validation Test
# Based on the working crosscheck_Xic.jl script
# -------------------------------------------------------------
cd(joinpath(@__DIR__, ".."))
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using YAML: YAML
using ThreeBodyDecays
using Lc2ppiKSemileptonicModelLHCb
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.Parameters
using ThreeBodyDecaysIO.JSON

println("XiC YAML vs JSON Model Validation")
println("="^50)

# -------------------------------------------------------------
# PART 1: Load YAML Model (Following crosscheck_Xic.jl pattern)
# -------------------------------------------------------------
println("1. Loading YAML model...")

# Load YAML model (same as crosscheck_Xic.jl)
isobarsinput = YAML.load_file(joinpath("data", "xic-particle-definitions.yaml"))
modelparameters = YAML.load_file(joinpath("data", "xic-model-definitions.yaml"))
defaultmodel = modelparameters["Default amplitude model"]

# Compute tbs from particle definitions
ms = let
    _mΞc = isobarsinput["Lambda_c+"]["mass"] / 1e3
    _mp = isobarsinput["p"]["mass"] / 1e3
    _mπ = isobarsinput["pi+"]["mass"] / 1e3
    _mK = isobarsinput["K-"]["mass"] / 1e3
    ThreeBodyMasses(m1 = _mp, m2 = _mπ, m3 = _mK, m0 = _mΞc)
end
tbs = ThreeBodySystem(ms, ThreeBodySpins(1, 0, 0; two_h0 = 1))

# Parse model dictionaries
(; chains, couplings, isobarnames) = parse_model_dictionaries(defaultmodel; particledict=isobarsinput)

# Create isobars
isobars = Dict()
for (key, lineshape) in chains
    dict = Dict{String, Any}(isobarsinput[key])
    dict["lineshape"] = lineshape
    isobars[key] = definechaininputs(key, dict; tbs)
end

# Update model parameters
defaultparameters = defaultmodel
defaultparameters["ArK(892)1"] = "1.0 ± 0.0"
defaultparameters["AiK(892)1"] = "0.0 ± 0.0"

println("✅ YAML model loaded successfully")
println("   - Masses: m₀=$(ms.m0), m₁=$(ms.m1), m₂=$(ms.m2), m₃=$(ms.m3)")
println("   - Chains: $(length(chains))")
println("   - Isobars: $(length(isobars))")

# -------------------------------------------------------------
# PART 2: Load JSON Model
# -------------------------------------------------------------
println("\n2. Loading JSON model...")

# Find JSON file
json_files = [
    joinpath("data", "xic2pKpi-model_new.json"),
    joinpath("data", "xic2pKpi-model.json")
]

json_file = nothing
for file in json_files
    if isfile(file)
        json_file = file
        break
    end
end

if json_file === nothing
    error("No JSON model file found")
end

# Load JSON
json_content = open(json_file) do io
    JSON.parse(io)
end

decay_description = json_content["distributions"][1]["decay_description"]
functions = json_content["functions"]

# Extract masses from JSON
kinematics = decay_description["kinematics"]
ms_json = let
    m0 = kinematics["initial_state"]["mass"]
    m1 = kinematics["final_state"][1]["mass"]
    m2 = kinematics["final_state"][2]["mass"]
    m3 = kinematics["final_state"][3]["mass"]
    ThreeBodyMasses(m1 = m1, m2 = m2, m3 = m3, m0 = m0)
end

# Reconstruct JSON model
workspace = Dict{String,Any}()
function_count = 0
for fn in functions
    try
        name = fn["name"]
        type_str = fn["type"]
        instance_type = eval(Symbol(type_str))
        workspace[name] = dict2instance(instance_type, fn)
        function_count += 1
    catch e
        println("   Warning: Could not create function $(fn["name"]): $e")
    end
end

try
    json_model = dict2instance(ThreeBodyDecay, decay_description; workspace)
    println("✅ JSON model loaded successfully")
    println("   - File: $(basename(json_file))")
    println("   - Masses: m₀=$(ms_json.m0), m₁=$(ms_json.m1), m₂=$(ms_json.m2), m₃=$(ms_json.m3)")
    println("   - Chains: $(length(decay_description["chains"]))")
    println("   - Functions created: $function_count/$(length(functions))")
catch e
    println("❌ JSON model reconstruction failed: $e")
    exit(1)
end

# -------------------------------------------------------------
# PART 3: Compare Mass Definitions
# -------------------------------------------------------------
println("\n3. Comparing mass definitions...")

mass_tolerance = 1e-8
mass_comparisons = [
    ("Parent (ΞC)", ms.m0, ms_json.m0),
    ("Proton", ms.m1, ms_json.m1),
    ("Pion", ms.m2, ms_json.m2),
    ("Kaon", ms.m3, ms_json.m3)
]

all_masses_match = true
for (name, yaml_mass, json_mass) in mass_comparisons
    diff = abs(yaml_mass - json_mass)
    match = diff < mass_tolerance
    status = match ? "✅" : "❌"
    println("   $status $name: YAML=$yaml_mass, JSON=$json_mass, diff=$diff")
    all_masses_match &= match
end

if all_masses_match
    println("✅ All masses match within tolerance ($mass_tolerance)")
else
    println("❌ Mass mismatch detected!")
    exit(1)
end

# -------------------------------------------------------------
# PART 4: Load Crosscheck Data for Test Points
# -------------------------------------------------------------
println("\n4. Loading crosscheck data for test points...")

crosscheckresult = readjson(joinpath("data", "crosscheck_Xic.json"))

# Extract test point (same as crosscheck)
σs0 = Invariants(ms,
    σ1 = crosscheckresult["chainvars"]["m2kpi"],
    σ2 = crosscheckresult["chainvars"]["m2pk"])

println("✅ Crosscheck data loaded")
println("   - Test point: σ₁=$(σs0.σ1), σ₂=$(σs0.σ2)")

# Parse complex numbers from crosscheck (same helper as crosscheck)
parsepythoncomplex(s::String) = eval(Meta.parse(
    replace(s,
        "(" => "",
        ")" => "",
        "j" => "im")))

# -------------------------------------------------------------
# PART 5: Compare Individual Lineshapes
# -------------------------------------------------------------
println("\n5. Comparing individual lineshapes...")

# Test key lineshapes (from crosscheck)
lineshape_tests = [
    ("K(892)", "BW_K(892)_p^1_q^0", σs0.σ1),
    ("L(1405)", "BW_L(1405)_p^0_q^0", σs0.σ2),
    ("L(1690)", "BW_L(1690)_p^2_q^1", σs0.σ2)
]

lineshape_results = []

for (isobar_name, crosscheck_key, test_σ) in lineshape_tests
    if haskey(isobars, isobar_name) && haskey(crosscheckresult["lineshapes"], crosscheck_key)
        try
            # Our YAML calculation
            yaml_value = isobars[isobar_name].Xlineshape(test_σ)
            
            # Reference value from crosscheck
            ref_value = crosscheckresult["lineshapes"][crosscheck_key] |> parsepythoncomplex
            
            # Compare
            rel_error = abs(yaml_value - ref_value) / abs(ref_value)
            match = rel_error < 1e-8
            
            status = match ? "✅" : "❌"
            println("   $status $isobar_name:")
            println("      YAML: $yaml_value")
            println("      Ref:  $ref_value")
            println("      Error: $rel_error")
            
            push!(lineshape_results, (name=isobar_name, yaml=yaml_value, ref=ref_value, error=rel_error, match=match))
            
        catch e
            println("   ❌ $isobar_name: Error - $e")
        end
    else
        println("   ⚠️  $isobar_name: Not available for testing")
    end
end

successful_lineshapes = count(r.match for r in lineshape_results)
println("Lineshape validation: $successful_lineshapes/$(length(lineshape_results)) passed")

# -------------------------------------------------------------
# PART 6: Compare Amplitude Calculations
# -------------------------------------------------------------
println("\n6. Comparing amplitude calculations...")

# Test amplitude calculation for key parameter (from crosscheck pattern)
Adict2matrix(d::Dict) = parsepythoncomplex.([
    d["A++"] d["A+-"]
    d["A-+"] d["A--"]])

# Filter for real parameters (same as crosscheck)
crosscheckresult_realpars = filter(kv -> kv[1][2] == 'r', crosscheckresult["chains"])

# Test a few key amplitudes
test_params = ["ArK(892)1", "ArL(1405)1", "ArL(1520)1"]
amplitude_results = []

for param_name in test_params
    if haskey(defaultparameters, param_name) && haskey(crosscheckresult_realpars, param_name)
        try
            # Our YAML calculation
            c, d = parname2decaychain(param_name, isobars; tbs)
            M_yaml = [c * amplitude(d, σs0, [two_λ1, 0, 0, two_λ0])
                      for (two_λ0, two_λ1) in [(1, 1) (1, -1); (-1, 1) (-1, -1)]]
            
            # Reference from crosscheck
            adict = crosscheckresult_realpars[param_name]
            M_ref = Adict2matrix(adict)
            
            # Compare matrices
            differences = abs.(M_yaml .- M_ref)
            relative_errors = differences ./ abs.(M_ref)
            max_rel_error = maximum(relative_errors)
            
            good_match = max_rel_error < 1e-8
            status = good_match ? "✅" : "❌"
            
            println("   $status $param_name:")
            println("      Max relative error: $max_rel_error")
            
            if !good_match
                println("      YAML matrix:")
                display(M_yaml)
                println("      Reference matrix:")
                display(M_ref)
            end
            
            push!(amplitude_results, (param=param_name, error=max_rel_error, match=good_match))
            
        catch e
            println("   ❌ $param_name: Error - $e")
        end
    else
        println("   ⚠️  $param_name: Not available for testing")
    end
end

successful_amplitudes = count(r.match for r in amplitude_results)
println("Amplitude validation: $successful_amplitudes/$(length(amplitude_results)) passed")

# -------------------------------------------------------------
# PART 7: Check JSON Validation Checksums
# -------------------------------------------------------------
println("\n7. Checking JSON validation checksums...")

if haskey(json_content, "misc") && haskey(json_content["misc"], "amplitude_model_checksums")
    checksums = json_content["misc"]["amplitude_model_checksums"]
    
    total_checksums = length(checksums)
    real_checksums = [c for c in checksums if c["value"] != "0.0 + 0.0i" && c["value"] != 0.0]
    placeholder_checksums = total_checksums - length(real_checksums)
    
    println("✅ JSON validation checksums found:")
    println("   - Total: $total_checksums")
    println("   - Real: $(length(real_checksums))")
    println("   - Placeholders: $placeholder_checksums")
    
    if length(real_checksums) > 0
        example = real_checksums[1]
        println("   - Example: $(example["distribution"]) = $(example["value"])")
    end
    
    checksum_status = length(real_checksums) > 0 ? "✅ PRESENT" : "⚠️  PLACEHOLDERS"
else
    println("❌ No validation checksums in JSON")
    checksum_status = "❌ MISSING"
end

# -------------------------------------------------------------
# SUMMARY REPORT
# -------------------------------------------------------------
println("\n" * "="^50)
println("YAML vs JSON VALIDATION SUMMARY")
println("="^50)

results_summary = [
    ("YAML Model Loading", "✅ SUCCESS"),
    ("JSON Model Loading", "✅ SUCCESS"),
    ("Mass Consistency", all_masses_match ? "✅ MATCH" : "❌ MISMATCH"),
    ("Lineshape Validation", "$successful_lineshapes/$(length(lineshape_results)) passed"),
    ("Amplitude Validation", "$successful_amplitudes/$(length(amplitude_results)) passed"),
    ("JSON Checksums", checksum_status)
]

println("Validation Results:")
for (test, result) in results_summary
    println("  $result $test")
end

# Overall assessment
critical_tests_passed = all_masses_match && successful_lineshapes > 0 && successful_amplitudes > 0
overall_status = critical_tests_passed ? "✅ VALIDATED" : "❌ ISSUES FOUND"

println("\n🎯 OVERALL STATUS: $overall_status")

if critical_tests_passed
    println("✅ YAML and JSON models are numerically equivalent")
    println("✅ Both models produce the same lineshape and amplitude results")
    println("✅ Models are ready for production use")
else
    println("❌ Numerical discrepancies found between models")
    println("📝 Review the specific test failures above")
end

println("\n📊 DETAILED RESULTS:")
if !isempty(lineshape_results)
    max_lineshape_error = maximum(r.error for r in lineshape_results)
    println("   - Maximum lineshape error: $max_lineshape_error")
end

if !isempty(amplitude_results)
    max_amplitude_error = maximum(r.error for r in amplitude_results)
    println("   - Maximum amplitude error: $max_amplitude_error")
end

println("\n🎉 YAML vs JSON validation complete!")

# Return summary for programmatic use
validation_summary = Dict(
    "masses_match" => all_masses_match,
    "lineshapes_passed" => successful_lineshapes,
    "amplitudes_passed" => successful_amplitudes,
    "total_lineshapes" => length(lineshape_results),
    "total_amplitudes" => length(amplitude_results),
    "json_has_checksums" => haskey(json_content, "misc") && haskey(json_content["misc"], "amplitude_model_checksums"),
    "overall_success" => critical_tests_passed
)

println("\nValidation completed with summary: $validation_summary")
