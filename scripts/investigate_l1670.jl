#!/usr/bin/env julia

"""
Investigate the L(1670) discrepancies in detail.
The L(1670) resonances show systematic differences with ratios of ~-3.35+0.77i
"""

using Pkg
Pkg.activate(".")

using JSON3, Printf

println("INVESTIGATING L(1670) DISCREPANCIES")
println("="^80)

# Load the detailed analysis results
data_file = "data/crosscheck_Xic.json"
if !isfile(data_file)
    error("Reference data file not found: $data_file")
end

reference_data = JSON3.read(read(data_file, String))

# Focus on L(1670) waves
l1670_waves = ["ArL(1670)1", "ArL(1670)2"]

println("\nAnalyzing L(1670) waves:")
for wave_name in l1670_waves
    println("\n", "="^50)
    println("WAVE: $wave_name")
    println("="^50)
    
    if haskey(reference_data, wave_name)
        ref_data = reference_data[wave_name]
        
        println("Reference matrix:")
        ref_matrix = ref_data["matrix"]
        for i in 1:2
            for j in 1:2
                val = ref_matrix[i][j]
                if isa(val, Dict)
                    real_part = get(val, "real", 0.0)
                    imag_part = get(val, "imag", 0.0)
                    @printf("%8.3f%+8.3fi  ", real_part, imag_part)
                else
                    @printf("%8.3f%+8.3fi  ", real(val), imag(val))
                end
            end
            println()
        end
        
        # Check if this is using L1670Flatte
        println("\nDecay chain information:")
        if haskey(ref_data, "decay_chain_type")
            println("Decay chain type: ", ref_data["decay_chain_type"])
        end
        
        # Look for any special parameters
        if haskey(ref_data, "parameters")
            println("Parameters:")
            for (key, value) in ref_data["parameters"]
                println("  $key: $value")
            end
        end
    else
        println("❌ Wave $wave_name not found in reference data!")
    end
end

println("\n" * "="^80)
println("ANALYSIS SUMMARY")
println("="^80)

println("""
The L(1670) resonances show systematic discrepancies with all matrix elements
having ratios of approximately -3.35 + 0.77i. This suggests:

1. A fundamental difference in the L(1670) implementation
2. Possible phase difference or sign convention
3. Different parametrization (Flatte vs BreitWigner)
4. Different coupling constants or normalization

Key observations:
- L(1670)1: 0/4 matches (complete mismatch)
- L(1670)2: 0/4 matches (complete mismatch)
- Consistent ratio pattern: -3.35 + 0.77i across all elements
- This is the only resonance type showing complete disagreement

Recommended investigation:
1. Check L1670Flatte implementation vs reference
2. Verify coupling constants and parameters
3. Check phase conventions
4. Compare with other Flatte resonances (L(1405) works fine)
""")
