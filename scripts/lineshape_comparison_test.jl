# =====================================================
# Direct Comparison Test: BreitWignerMinL vs MultichannelBreitWigner
# =====================================================

import Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

using Lc2ppiKSemileptonicModelLHCb
using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.ThreeBodyDecays
using ThreeBodyDecaysIO.HadronicLineshapes
using Statistics
using Printf

# Include local files
include(joinpath(@__DIR__, "..", "shapes.jl"))
include(joinpath(@__DIR__, "..", "src", "lineshapes.jl"))

function test_lineshape_comparison()
    println("="^60)
    println("LINESHAPE COMPARISON: BreitWignerMinL vs MultichannelBreitWigner")
    println("="^60)
    
    # Define test parameters for L(1600)
    m = 1.6        # GeV - Resonance mass
    Γ₀ = 0.2       # GeV - Nominal width
    l = 1          # Orbital angular momentum
    minL = 0       # Minimal angular momentum for decay
    
    # Particle masses (GeV)
    mp = 0.938272046    # Proton
    mK = 0.49367700000000003  # Kaon
    mπ = 0.13957018     # Pion
    mΛc = 2.28646       # Lambda_c
    
    # Form factor radii
    dR = 1.5    # Resonance radius
    dΛc = 5.0   # Lambda_c radius
    
    println("Test Parameters:")
    println("  Resonance: L(1600)")
    println("  Mass: $m GeV")
    println("  Width: $Γ₀ GeV") 
    println("  l: $l, minL: $minL")
    println("  Masses: p=$mp, K=$mK, π=$mπ, Λc=$mΛc GeV")
    println("  Radii: dR=$dR, dΛc=$dΛc GeV⁻¹")
    println()
    
    # ================================================
    # 1. Create YAML BreitWignerMinL
    # ================================================
    yaml_bw = BreitWignerMinL(
        pars = (m, Γ₀),
        l = l,
        minL = minL,
        name = "L(1600)",
        m1 = mp,      # Proton
        m2 = mK,      # Kaon
        mk = mπ,      # Pion (spectator)
        m0 = mΛc      # Lambda_c
    )
    
    # ================================================
    # 2. Create JSON MultichannelBreitWigner (Method 1: Direct)
    # ================================================
    p0 = breakup(m, mp, mK)
    gsq_direct = m * Γ₀ / (2 * p0)  # Direct conversion
    
    json_bw_direct = MultichannelBreitWigner(
        m,
        SVector([(gsq = gsq_direct, ma = mp, mb = mK, l = l, d = dR)])
    )
    
    # ================================================
    # 3. Create JSON MultichannelBreitWigner (Method 2: Form-factor corrected)
    # ================================================
    FF = HadronicLineshapes.BlattWeisskopf{l}(dR)
    FF_at_pole = FF(p0)
    gsq_corrected = m * Γ₀ / (2 * p0) / FF_at_pole^2
    
    json_bw_corrected = MultichannelBreitWigner(
        m,
        SVector([(gsq = gsq_corrected, ma = mp, mb = mK, l = l, d = dR)])
    )
    
    # ================================================
    # 4. Create JSON MultichannelBreitWigner (Method 3: Empirical)
    # ================================================
    gsq_empirical = Γ₀ * p0 / m  # Your empirically found formula
    
    json_bw_empirical = MultichannelBreitWigner(
        m,
        SVector([(gsq = gsq_empirical, ma = mp, mb = mK, l = l, d = dR)])
    )
    
    println("Conversion Parameters:")
    println("  p0 = $p0 GeV")
    println("  FF(p0) = $FF_at_pole")
    println("  gsq_direct = $gsq_direct")
    println("  gsq_corrected = $gsq_corrected") 
    println("  gsq_empirical = $gsq_empirical")
    println()
    
    # ================================================
    # 5. Test at different energy points
    # ================================================
    test_energies = [1.4, 2.0, 2.5, 3.0, 3.2, 4.0]  # GeV²
    
    println("Energy Point Tests:")
    println("-"^80)
    @printf "%-8s %-20s %-20s %-20s %-20s %-20s\n" "s [GeV²]" "YAML" "JSON_direct" "JSON_corrected" "JSON_empirical" "YAML/JSON_ratio"
    println("-"^80)
    
    for s in test_energies
        # Calculate lineshape values
        yaml_value = yaml_bw(s)
        json_direct_value = json_bw_direct(s)
        json_corrected_value = json_bw_corrected(s) 
        json_empirical_value = json_bw_empirical(s)
        
        # Calculate ratio of magnitudes
        ratio_direct = abs(yaml_value) / abs(json_direct_value)
        ratio_corrected = abs(yaml_value) / abs(json_corrected_value)
        ratio_empirical = abs(yaml_value) / abs(json_empirical_value)
        
        @printf "%-8.1f %-20s %-20s %-20s %-20s %-8.3f/%-8.3f/%-8.3f\n" s "$(round(abs(yaml_value), digits=4))" "$(round(abs(json_direct_value), digits=4))" "$(round(abs(json_corrected_value), digits=4))" "$(round(abs(json_empirical_value), digits=4))" ratio_direct ratio_corrected ratio_empirical
    end
    
    println("-"^80)
    println()
    
    # ================================================
    # 6. Detailed analysis at s = 3.2 GeV²
    # ================================================
    s_test = 3.2
    println("DETAILED ANALYSIS at s = $s_test GeV²:")
    println("-"^50)
    
    yaml_detailed = yaml_bw(s_test)
    json_direct_detailed = json_bw_direct(s_test)
    json_corrected_detailed = json_bw_corrected(s_test)
    json_empirical_detailed = json_bw_empirical(s_test)
    
    println("Complex Values:")
    println("  YAML:           $yaml_detailed")
    println("  JSON_direct:    $json_direct_detailed")
    println("  JSON_corrected: $json_corrected_detailed")
    println("  JSON_empirical: $json_empirical_detailed")
    println()
    
    println("Magnitudes:")
    println("  |YAML|:           $(abs(yaml_detailed))")
    println("  |JSON_direct|:    $(abs(json_direct_detailed))")
    println("  |JSON_corrected|: $(abs(json_corrected_detailed))")
    println("  |JSON_empirical|: $(abs(json_empirical_detailed))")
    println()
    
    println("Ratios (YAML/JSON):")
    println("  Direct:    $(yaml_detailed / json_direct_detailed)")
    println("  Corrected: $(yaml_detailed / json_corrected_detailed)")
    println("  Empirical: $(yaml_detailed / json_empirical_detailed)")
    println()
    
    # ================================================
    # 7. Form Factor Analysis
    # ================================================
    println("FORM FACTOR ANALYSIS:")
    println("-"^30)
    
    sqrt_s = sqrt(s_test)
    p_current = breakup(sqrt_s, mp, mK)
    q_current = breakup(mΛc, sqrt_s, mπ)
    q0 = breakup(mΛc, m, mπ)
    
    FF_current = FF(p_current)
    kinetic_factor = (p_current/p0)^l * (q_current/q0)^minL
    
    println("  √s = $sqrt_s GeV")
    println("  p(s) = $p_current GeV")
    println("  p0 = $p0 GeV") 
    println("  q(s) = $q_current GeV")
    println("  q0 = $q0 GeV")
    println("  FF(p) = $FF_current")
    println("  Kinetic factor = $kinetic_factor")
    println()
    
    # Test form factor correction
    yaml_ff_corrected = yaml_detailed / (FF_current^2 * kinetic_factor)
    println("  YAML / (FF² × kinetic) = $yaml_ff_corrected")
    println("  Difference from JSON_direct: $(abs(yaml_ff_corrected - json_direct_detailed))")
    println()
    
    # ================================================
    # 8. Summary and Recommendations
    # ================================================
    println("SUMMARY:")
    println("-"^20)
    
    # Find best matching method
    diffs = [
        abs(yaml_detailed - json_direct_detailed),
        abs(yaml_detailed - json_corrected_detailed), 
        abs(yaml_detailed - json_empirical_detailed)
    ]
    
    best_method = argmin(diffs)
    methods = ["Direct", "Form-factor corrected", "Empirical"]
    
    println("  Best matching method: $(methods[best_method])")
    println("  Smallest difference: $(diffs[best_method])")
    println()
    
    if diffs[best_method] < 1e-10
        println("  ✓ Excellent agreement found!")
    elseif diffs[best_method] < 1e-6
        println("  ✓ Good agreement found!")
    else
        println("  ⚠ Significant differences remain")
        println("    Recommended gsq for perfect match: $(gsq_empirical * abs(yaml_detailed) / abs(json_empirical_detailed))")
    end
    
    println()
    println("="^60)
    println("TEST COMPLETED")
    println("="^60)
end

# Run the test
if abspath(PROGRAM_FILE) == @__FILE__
    test_lineshape_comparison()
end

# Also provide a simple function call
test_lineshape_comparison()
