# Simple test to validate YAML vs JSON model equivalence
println("Starting XiC model comparison...")

try
    using Lc2ppiKSemileptonicModelLHCb
    using Lc2ppiKSemileptonicModelLHCb.ThreeBodyDecays
    using ThreeBodyDecaysIO
    using YAML
    println("✓ Packages loaded successfully")
    
    # Load YAML model
    println("Loading YAML model...")
    particledict = YAML.load_file("data/xic-particle-definitions.yaml")
    modelparameters = YAML.load_file("data/xic-model-definitions.yaml")
    defaultparameters = modelparameters["Default amplitude model"]
    
    (; chains, couplings, isobarnames) = parse_model_dictionaries(defaultparameters; particledict)
    yaml_model = Lc2ppiKModel(; chains, couplings, isobarnames)
    println("✓ YAML model created with $(length(chains)) chains")
    
    # Test point
    ms = yaml_model.tbs.ms
    test_point = Invariants(ms, σ1 = 1.9^2, σ2 = 1.6^2)
    helicities = [1, 0, 0, 1]
    
    # Calculate amplitude
    amp = amplitude(yaml_model, test_point, helicities)
    intensity = abs2(amp)
    
    println("Test Results:")
    println("  Amplitude: $amp")
    println("  Intensity: $intensity")
    println("  Point: σ₁=$(test_point.σ1), σ₂=$(test_point.σ2)")
    
catch e
    println("❌ Error: $e")
    println("Stacktrace:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end
