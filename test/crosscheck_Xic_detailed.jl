cd(joinpath(@__DIR__, ".."))
using Pkg
Pkg.activate(".")
Pkg.instantiate()
#
using YAML: YAML
using Plots
using LaTeXStrings
import Plots.PlotMeasures.mm
#
using Parameters
using Measurements
using DataFrames
#
using ThreeBodyDecays

using Lc2ppiKSemileptonicModelLHCb

# Load the same setup as crosscheck_Xic.jl
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

# Define isobars
isobars = Dict()
for (key, lineshape) in defaultmodel["lineshapes"]
    dict = Dict{String, Any}(isobarsinput[key])
    dict["lineshape"] = lineshape
    isobars[key] = definechaininputs(key, dict; tbs)
end

# Update model parameters
defaultparameters = defaultmodel["parameters"]
defaultparameters["ArK(892)1"] = "1.0 ± 0.0"
defaultparameters["AiK(892)1"] = "0.0 ± 0.0"

# Note: Xic model doesn't have the same shape parameters as Lc model
# so we skip the parameter updates for now

# Load crosscheck data
crosscheckresult = readjson(joinpath("data", "crosscheck_Xic.json"))

σs0 = Invariants(ms,
    σ1 = crosscheckresult["chainvars"]["m2kpi"],
    σ2 = crosscheckresult["chainvars"]["m2pk"])

parsepythoncomplex(s::String) = eval(Meta.parse(
    replace(s,
        "(" => "",
        ")" => "",
        "j" => "im")))

Adict2matrix(d::Dict) = parsepythoncomplex.(
    [d["A++"] d["A+-"]
                 d["A-+"] d["A--"]])

# Detailed analysis function
function analyze_wave(parname, adict, isobars, tbs, σs0)
    println("="^80)
    println("ANALYZING WAVE: $parname")
    println("="^80)

    # Get our calculation
    c, d = parname2decaychain(parname, isobars; tbs)
    println("Coupling factor: $c")
    println("Decay chain: $(typeof(d))")

    # Calculate our matrix elements
    M_DPD = [c * amplitude(d, σs0, [two_λ1, 0, 0, two_λ0])
             for (two_λ0, two_λ1) in
             [(1, 1) (1, -1)
        (-1, 1) (-1, -1)]]

    # Get reference matrix
    M_LHCb′ = amplitudeLHCb2DPD(Adict2matrix(adict))

    println("\nMATRIX ELEMENTS COMPARISON:")
    println("Our DPD calculation:")
    display(M_DPD)
    println("\nReference LHCb' calculation:")
    display(M_LHCb′)

    # Improved comparison logic
    tolerance = 1e-3
    zerotol = 1e-10
    ratios = similar(M_DPD)
    matches = falses(2, 2)
    redflags = falses(2, 2)
    for i in 1:2, j in 1:2
        our_val = M_DPD[i, j]
        ref_val = M_LHCb′[i, j]
        if abs(our_val) < zerotol && abs(ref_val) < zerotol
            matches[i, j] = true
            ratios[i, j] = 1.0
        elseif abs(our_val) < zerotol || abs(ref_val) < zerotol
            matches[i, j] = false
            redflags[i, j] = true
            ratios[i, j] = NaN
        else
            ratios[i, j] = our_val / ref_val
            matches[i, j] = abs(ratios[i, j] - 1) < tolerance
        end
    end

    println("\nRATIOS (Our/Reference):")
    display(ratios)
    println("\nELEMENTS MATCHING (within $tolerance):")
    display(matches)
    println("\nRED FLAG ELEMENTS (one zero, one nonzero):")
    display(redflags)

    # Detailed element-by-element analysis
    println("\nDETAILED ELEMENT ANALYSIS:")
    for i in 1:2, j in 1:2
        λ0 = i == 1 ? 1 : -1
        λ1 = j == 1 ? 1 : -1
        our_val = M_DPD[i, j]
        ref_val = M_LHCb′[i, j]
        ratio = ratios[i, j]
        matches_el = matches[i, j]
        redflag_el = redflags[i, j]
        println("Element [λ0=$λ0, λ1=$λ1]:")
        println("  Our value:     $our_val")
        println("  Reference:     $ref_val")
        println("  Ratio:         $ratio")
        println("  Matches:       $matches_el")
        println("  Red flag:      $redflag_el")
        println("  Difference:    $(abs(our_val - ref_val))")
        println()
    end

    # Helicity flip investigation (permutations)
    perms = [M_DPD, M_DPD', reverse(M_DPD, dims = 1), reverse(M_DPD, dims = 2)]
    permnames = ["original", "transpose", "flip rows", "flip cols"]
    best_match = 0
    best_perm = 1
    for (k, perm) in enumerate(perms)
        perm_matches = falses(2, 2)
        for i in 1:2, j in 1:2
            our_val = perm[i, j]
            ref_val = M_LHCb′[i, j]
            if abs(our_val) < zerotol && abs(ref_val) < zerotol
                perm_matches[i, j] = true
            elseif abs(our_val) < zerotol || abs(ref_val) < zerotol
                perm_matches[i, j] = false
            else
                perm_matches[i, j] = abs(our_val / ref_val - 1) < tolerance
            end
        end
        nmatch = count(perm_matches)
        if nmatch > best_match
            best_match = nmatch
            best_perm = k
        end
    end
    println("\nHELICITY FLIP INVESTIGATION:")
    println("Best permutation: $(permnames[best_perm]) with $best_match/4 matches")
    println("Permutation matrix:")
    display(perms[best_perm])

    return (parname = parname, M_DPD = M_DPD, M_LHCb′ = M_LHCb′, ratios = ratios, matches = matches, redflags = redflags, best_perm = permnames[best_perm], best_match = best_match)
end

# Analyze all waves
crosscheckresult_realpars = filter(kv -> kv[1][2] == 'r', crosscheckresult["chains"])

println("DETAILED XIC CROSSCHECK ANALYSIS")
println("="^80)
println("Total waves to analyze: $(length(crosscheckresult_realpars))")
println()

# Store results for summary
analysis_results = []

for (parname, adict) in crosscheckresult_realpars
    result = analyze_wave(parname, adict, isobars, tbs, σs0)
    push!(analysis_results, result)

    # Check if any elements don't match
    if !all(result.matches)
        println("⚠️  WAVE HAS MISMATCHES: $parname")
        println("   Non-matching elements: $(findall(.!result.matches))")
    else
        println("✅ WAVE MATCHES PERFECTLY: $parname")
    end
    println()
end

# Summary statistics
println("SUMMARY STATISTICS")
println("="^80)

total_elements = length(analysis_results) * 4
matching_elements = sum(sum(r.matches) for r in analysis_results)
mismatching_elements = total_elements - matching_elements

println("Total matrix elements: $total_elements")
println("Matching elements: $matching_elements")
println("Mismatching elements: $mismatching_elements")
println("Match rate: $(round(matching_elements/total_elements*100, digits=2))%")

# Waves with issues
problematic_waves = [r.parname for r in analysis_results if !all(r.matches)]
println("\nWaves with mismatches: $(length(problematic_waves))")
for wave in problematic_waves
    println("  - $wave")
end

# Save detailed results
using JSON3
detailed_results = Dict(
    "summary" => Dict(
        "total_elements" => total_elements,
        "matching_elements" => matching_elements,
        "mismatching_elements" => mismatching_elements,
        "match_rate" => matching_elements / total_elements,
        "problematic_waves" => problematic_waves,
    ),
    "wave_analysis" => Dict(
        r.parname => Dict(
            "M_DPD" => [string(x) for x in r.M_DPD],
            "M_LHCb" => [string(x) for x in r.M_LHCb′],
            "ratios" => [string(x) for x in r.ratios],
            "matches" => r.matches,
            "redflags" => r.redflags,
            "all_match" => all(r.matches),
            "best_perm" => r.best_perm,
            "best_match" => r.best_match,
        ) for r in analysis_results
    ),
)

open(joinpath("research", "detailed_xic_analysis.json"), "w") do f
    write(f, JSON3.write(detailed_results, pretty = true))
end

println("\nDetailed results saved to research/detailed_xic_analysis.json")
