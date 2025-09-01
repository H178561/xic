#!/usr/bin/env julia
# Script to compare BreitWignerMinL vs MultichannelBreitWigner for L(1600)
# This helps understand the differences and fix the conversion

import Pkg
Pkg.activate(dirname(@__DIR__))
Pkg.instantiate()

using Lc2ppiKSemileptonicModelLHCb
using ThreeBodyDecaysIO
using ThreeBodyDecaysIO.HadronicLineshapes
using ThreeBodyDecays
using Plots
using YAML

# Load particle definitions to get masses
particledict = YAML.load_file(joinpath(dirname(@__DIR__), "data", "xic-particle-definitions.yaml"))

# Extract masses (convert MeV to GeV)
m_Xic = 2.46794  # 2.46794 GeV
m_p = particledict["p"]["mass"] / 1000.0      # 0.938272 GeV  
m_pi = particledict["pi+"]["mass"] / 1000.0    # 0.13957 GeV
m_K = particledict["K-"]["mass"] / 1000.0      # 0.493677 GeV

# L(1600) parameters from particle definitions
m_L1600 = particledict["L(1600)"]["mass"] / 1000.0  # 1.6 GeV
Γ_L1600 = particledict["L(1600)"]["width"] / 1000.0 # 0.2 GeV

println("Particle masses (GeV):")
println("  Ξc: $m_Xic")
println("  p:  $m_p") 
println("  π:  $m_pi")
println("  K:  $m_K")
println("  L(1600): mass = $m_L1600 GeV, width = $Γ_L1600 GeV")
println()

# =============================================================================
# 1. Create BreitWignerMinL (YAML-style)
# =============================================================================

# L(1600) decays to p+K with L=1, and is produced in Ξc decay with minL=0
include(joinpath(dirname(@__DIR__), "src", "lineshapes.jl"))

# Create YAML-style BreitWignerMinL
yaml_lineshape = BreitWignerMinL(
    pars = (m_L1600, Γ_L1600),
    l = 1,       # orbital angular momentum in p+K system
    minL = 0,    # minimum L in Ξc → L(1600) π decay
    name = "L(1600)",
    m1 = m_p,    # proton mass
    m2 = m_K,    # kaon mass  
    mk = m_pi,   # pion mass (bachelor particle)
    m0 = m_Xic   # parent mass
)

println("YAML-style BreitWignerMinL created:")
println("  pars: $(yaml_lineshape.pars)")
println("  l: $(yaml_lineshape.l)")
println("  minL: $(yaml_lineshape.minL)")
println("  m1 (p): $(yaml_lineshape.m1)")
println("  m2 (K): $(yaml_lineshape.m2)")
println("  mk (π): $(yaml_lineshape.mk)")
println("  m0 (Ξc): $(yaml_lineshape.m0)")
println()

# =============================================================================
# 2. Create MultichannelBreitWigner (JSON-style)
# =============================================================================

# Calculate effective coupling for MultichannelBreitWigner
# The gsq should be chosen so that at the resonance pole, the width matches Γ_L1600

# Reference momentum at resonance pole
p0 = breakup(m_L1600, m_p, m_K)
dR = 1.5  # Blatt-Weisskopf radius

# Form factor at pole
FF_pole = HadronicLineshapes.BlattWeisskopf{1}(dR)
FF_at_pole = FF_pole(p0)

# Calculate effective coupling
# From MultichannelBreitWigner: mΓ = gsq * 2*p/√s * (p/p0)^(2l+1) * FF(p)²/FF(p0)²
# At pole (s = m²): Γ = gsq * 2*p0/m * FF(p0)²/FF(p0)² = gsq * 2*p0/m
# Therefore: gsq = m * Γ / (2*p0)
gsq_eff = m_L1600 * Γ_L1600 / (2 * p0)

# Create JSON-style MultichannelBreitWigner
json_lineshape = MultichannelBreitWigner(
    m_L1600,
    SVector((
        gsq = gsq_eff,
        ma = m_p,
        mb = m_K, 
        l = 1,
        d = dR
    ))
)

println("JSON-style MultichannelBreitWigner created:")
println("  mass: $(json_lineshape.m)")
println("  gsq_eff: $gsq_eff")
println("  reference momentum p0: $p0")
println("  channel: ma=$m_p, mb=$m_K, l=1, d=$dR")
println()

# =============================================================================
# 3. Compare lineshapes over energy range
# =============================================================================

# Define energy range for comparison
s_min = (m_p + m_K)^2    # Threshold
s_max = (m_Xic - m_pi)^2 # Maximum allowed
s_range = range(s_min, s_max, length=1000)

# Calculate lineshape values
yaml_values = [yaml_lineshape(s) for s in s_range]
json_values = [json_lineshape(s) for s in s_range]

println("Comparison at specific points:")
test_points = [1.5, 2.0, 2.5, 3.0, 3.2, 3.5]

for s_test in test_points
    if s_test >= s_min && s_test <= s_max
        yaml_val = yaml_lineshape(s_test)
        json_val = json_lineshape(s_test)
        
        println("s = $s_test GeV²:")
        println("  YAML: $yaml_val")
        println("  JSON: $json_val")
        println("  |YAML|: $(abs(yaml_val))")
        println("  |JSON|: $(abs(json_val))")
        println("  Ratio: $(yaml_val / json_val)")
        println("  |Ratio|: $(abs(yaml_val) / abs(json_val))")
        println()
    end
end

# =============================================================================
# 4. Detailed breakdown at s = 3.2 GeV²
# =============================================================================

s_test = 3.2
sqrt_s = sqrt(s_test)

println("="^60)
println("DETAILED BREAKDOWN AT s = $s_test GeV²")
println("="^60)

# Current momenta
p_current = breakup(sqrt_s, m_p, m_K)
q_current = breakup(m_Xic, sqrt_s, m_pi)

# Reference momenta  
p0 = breakup(m_L1600, m_p, m_K)
q0 = breakup(m_Xic, m_L1600, m_pi)

println("Momenta:")
println("  p_current = $p_current GeV")
println("  p0 = $p0 GeV") 
println("  q_current = $q_current GeV")
println("  q0 = $q0 GeV")
println()

# Form factors
dR = 1.5    # resonance radius
dΛc = 5.0   # decay radius
FF_res = HadronicLineshapes.BlattWeisskopf{1}(dR)
FF_decay = HadronicLineshapes.BlattWeisskopf{0}(dΛc)  # minL = 0

FF_res_current = FF_res(p_current)
FF_res_pole = FF_res(p0)
FF_decay_current = FF_decay(q_current)
FF_decay_pole = FF_decay(q0)

println("Form factors:")
println("  FF_res(p_current) = $FF_res_current")
println("  FF_res(p0) = $FF_res_pole")
println("  FF_decay(q_current) = $FF_decay_current")  
println("  FF_decay(q0) = $FF_decay_pole")
println()

# YAML formula breakdown
println("YAML BreitWignerMinL breakdown:")

# 1. Momentum-dependent width
Γ_current = Γ_L1600 * (p_current/p0)^3 * m_L1600/sqrt_s * (FF_res_current/FF_res_pole)^2
println("  Γ_current = $Γ_current GeV")

# 2. Basic Breit-Wigner
basic_BW = 1 / (m_L1600^2 - s_test - 1im * m_L1600 * Γ_current)
println("  Basic BW = $basic_BW")

# 3. Additional factors
momentum_factor = (p_current/p0)^1 * (q_current/q0)^0  # l=1, minL=0
ff_factor = sqrt((FF_res_current/FF_res_pole)^2 * (FF_decay_current/FF_decay_pole)^2)

println("  Momentum factor = $momentum_factor")
println("  FF factor = $ff_factor")

# Complete YAML result
yaml_complete = basic_BW * momentum_factor * ff_factor
yaml_actual = yaml_lineshape(s_test)

println("  Complete calculation = $yaml_complete")
println("  Actual YAML result = $yaml_actual")
println("  Match: $(abs(yaml_complete - yaml_actual) < 1e-10)")
println()

# JSON formula breakdown
println("JSON MultichannelBreitWigner breakdown:")

# Width calculation
json_momentum_factor = (p_current/p0)^3  # (2l+1) = 3
json_ff_ratio = (FF_res_current/FF_res_pole)^2
mΓ_json = gsq_eff * 2*p_current/sqrt_s * json_momentum_factor * json_ff_ratio

println("  mΓ_json = $mΓ_json GeV")
println("  Γ_json = $(mΓ_json/m_L1600) GeV")

json_complete = 1 / (m_L1600^2 - s_test - 1im * m_L1600 * mΓ_json/m_L1600)
json_actual = json_lineshape(s_test)

println("  Complete calculation = $json_complete")
println("  Actual JSON result = $json_actual")
println("  Match: $(abs(json_complete - json_actual) < 1e-10)")
println()

# =============================================================================
# 5. Analysis and correction
# =============================================================================

println("="^60)
println("ANALYSIS")
println("="^60)

println("The key difference is in the additional kinematic factors:")
println("  YAML includes: (p/p0)^l * (q/q0)^minL * √[FF_res * FF_decay] = $momentum_factor * $ff_factor")
println("  JSON only includes the basic Breit-Wigner with momentum-dependent width")
println()

println("To make them match, the JSON lineshape should be multiplied by:")
correction_factor = momentum_factor * ff_factor
println("  Correction factor = $correction_factor")
println("  JSON corrected = $(json_actual * correction_factor)")
println("  YAML actual = $yaml_actual")
println("  Ratio after correction = $(yaml_actual / (json_actual * correction_factor))")

# =============================================================================
# 6. Create plots
# =============================================================================

println()
println("Creating comparison plots...")

# Plot 1: Magnitude comparison
p1 = plot(sqrt.(s_range), abs.(yaml_values), label="YAML BreitWignerMinL", lw=2)
plot!(p1, sqrt.(s_range), abs.(json_values), label="JSON MultichannelBW", lw=2, ls=:dash)
plot!(p1, sqrt.(s_range), abs.(json_values) .* abs(correction_factor), label="JSON corrected", lw=2, ls=:dot)
xlabel!(p1, "√s (GeV)")
ylabel!(p1, "|Lineshape|")
title!(p1, "L(1600) Lineshape Magnitude")

# Plot 2: Phase comparison  
p2 = plot(sqrt.(s_range), angle.(yaml_values), label="YAML BreitWignerMinL", lw=2)
plot!(p2, sqrt.(s_range), angle.(json_values), label="JSON MultichannelBW", lw=2, ls=:dash)
xlabel!(p2, "√s (GeV)")
ylabel!(p2, "Phase (rad)")
title!(p2, "L(1600) Lineshape Phase")

# Plot 3: Real and imaginary parts
p3 = plot(sqrt.(s_range), real.(yaml_values), label="YAML Real", lw=2)
plot!(p3, sqrt.(s_range), real.(json_values), label="JSON Real", lw=2, ls=:dash)
plot!(p3, sqrt.(s_range), imag.(yaml_values), label="YAML Imag", lw=2)
plot!(p3, sqrt.(s_range), imag.(json_values), label="JSON Imag", lw=2, ls=:dash)
xlabel!(p3, "√s (GeV)")
ylabel!(p3, "Lineshape")
title!(p3, "L(1600) Real and Imaginary Parts")

# Combine plots
final_plot = plot(p1, p2, p3, layout=(3,1), size=(800, 900))

# Save plot
plot_path = joinpath(dirname(@__DIR__), "plots", "l1600_lineshape_comparison.png")
savefig(final_plot, plot_path)
println("Plot saved to: $plot_path")

println()
println("="^60)
println("SUMMARY")
println("="^60)
println("1. YAML BreitWignerMinL includes additional kinematic factors")
println("2. JSON MultichannelBreitWigner only has the basic resonance shape")
println("3. To match them, multiply JSON result by: $correction_factor")
println("4. This factor comes from: (p/p0)^l * (q/q0)^minL * √[FF_corrections]")
println("5. The conversion needs to account for these vertex form factors")
