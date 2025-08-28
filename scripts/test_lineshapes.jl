#!/usr/bin/env julia

"""
Test script for individual lineshape functions with Mandelstam variables
"""

using Lc2ppiKSemileptonicModelLHCb
using Plots
using Printf

# Create test kinematic setup
const ms = ThreeBodyMasses(m1=0.938, m2=0.494, m3=0.140, m0=2.286)  # p, K, π, Λc
const tbs = ThreeBodySystem(ms, ThreeBodySpins(1, 0, 0; two_h0=1))

println("Testing individual lineshape functions")
println("====================================")

# Test 1: Simple BreitWigner lineshape
println("\n1. Testing BreitWigner lineshape for K(892)")

# Create a K(892) BreitWigner lineshape
k892_bw = BreitWignerMinL(
    pars=(0.8955, 0.0473),  # mass, width
    l=1,                    # orbital angular momentum 
    minL=1,                 # minimal L
    name="K(892)_test",
    m1=0.140,              # π mass
    m2=0.494,              # K mass  
    mk=0.938,              # proton mass (spectator)
    m0=2.286               # Λc mass
)

# Test over a range of σ values (Mandelstam variable s_23)
σ_range = range(0.5, 1.5, length=1000)  # GeV²
amplitudes = [k892_bw(σ) for σ in σ_range]

# Plot the lineshape
plot(sqrt.(σ_range), abs2.(amplitudes), 
     xlabel="m(πK) [GeV]", 
     ylabel="|A|²", 
     title="K(892) BreitWigner Lineshape",
     label="K(892) BW",
     lw=2)

# Test specific Mandelstam values
test_values = [0.6, 0.8, 1.0, 1.2]
println("\nK(892) BreitWigner test values:")
for σ in test_values
    amp = k892_bw(σ)
    @printf("σ = %.2f GeV²  →  A = %.4f + %.4fi  →  |A|² = %.6f\n", 
            σ, real(amp), imag(amp), abs2(amp))
end

# Test 2: Λ(1520) BreitWigner
println("\n2. Testing BreitWigner lineshape for Λ(1520)")

l1520_bw = BreitWignerMinL(
    pars=(1.5185, 0.0152),  # mass, width from your JSON
    l=2,                    # orbital angular momentum
    minL=2,                 # minimal L
    name="L(1520)_test", 
    m1=0.938,              # proton mass
    m2=0.494,              # K mass
    mk=0.140,              # π mass (spectator) 
    m0=2.286               # Λc mass
)

# Test over different σ range for Λ(1520)
σ_range_l = range(2.0, 3.0, length=1000)  # GeV² for pK system
amplitudes_l = [l1520_bw(σ) for σ in σ_range_l]

plot!(sqrt.(σ_range_l), abs2.(amplitudes_l), 
      xlabel="m(pK) [GeV]", 
      ylabel="|A|²",
      label="Λ(1520) BW",
      lw=2)

println("\nΛ(1520) BreitWigner test values:")
test_values_l = [2.2, 2.3, 2.4, 2.5]
for σ in test_values_l
    amp = l1520_bw(σ)
    @printf("σ = %.2f GeV²  →  A = %.4f + %.4fi  →  |A|² = %.6f\n", 
            σ, real(amp), imag(amp), abs2(amp))
end

# Test 3: Compare with complex Mandelstam variables (for phase space)
println("\n3. Testing with complex Mandelstam variables")

# Test near the pole with small imaginary part
σ_complex = 0.8955^2 + 0.01im  # Near K(892) pole
amp_complex = k892_bw(σ_complex)
@printf("σ = %.4f + %.4fi  →  A = %.6f + %.6fi\n", 
        real(σ_complex), imag(σ_complex), real(amp_complex), imag(amp_complex))

# Test 4: Physical phase space check
println("\n4. Testing physical phase space boundaries")

function test_physical_boundaries(lineshape, name, σ_min, σ_max)
    println("\nTesting $name boundaries:")
    
    # Test at threshold
    amp_min = lineshape(σ_min)
    @printf("At threshold σ_min = %.4f: A = %.6f + %.6fi, |A|² = %.8f\n",
            σ_min, real(amp_min), imag(amp_min), abs2(amp_min))
    
    # Test at upper limit 
    amp_max = lineshape(σ_max)
    @printf("At upper limit σ_max = %.4f: A = %.6f + %.6fi, |A|² = %.8f\n",
            σ_max, real(amp_max), imag(amp_max), abs2(amp_max))
    
    # Test at resonance
    m_res = sqrt(lineshape.pars[1]^2)  # Resonance mass
    σ_res = m_res^2
    if σ_min <= σ_res <= σ_max
        amp_res = lineshape(σ_res)
        @printf("At resonance σ_res = %.4f: A = %.6f + %.6fi, |A|² = %.8f\n",
                σ_res, real(amp_res), imag(amp_res), abs2(amp_res))
    end
end

# Test K(892) boundaries (πK system)
σ_min_k892 = (0.140 + 0.494)^2  # (mπ + mK)²
σ_max_k892 = (2.286 - 0.938)^2  # (mΛc - mp)²
test_physical_boundaries(k892_bw, "K(892)", σ_min_k892, σ_max_k892)

# Test Λ(1520) boundaries (pK system)  
σ_min_l1520 = (0.938 + 0.494)^2  # (mp + mK)²
σ_max_l1520 = (2.286 - 0.140)^2  # (mΛc - mπ)²
test_physical_boundaries(l1520_bw, "Λ(1520)", σ_min_l1520, σ_max_l1520)

savefig("lineshape_test.png")
println("\nPlot saved as 'lineshape_test.png'")

# Test 5: Function evaluation timing
println("\n5. Performance test")
using BenchmarkTools

println("Timing K(892) BreitWigner evaluation:")
σ_test = 0.9
@btime $k892_bw($σ_test)

println("\nTiming Λ(1520) BreitWigner evaluation:")
σ_test_l = 2.3
@btime $l1520_bw($σ_test_l)

println("\nLineshape testing complete!")
