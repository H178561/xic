function Lc2ppiKModel(; chains, couplings, isobarnames)
    v = collect(enumerate(isobarnames))
    # parse K2(number) or K(number), or L(number), or S(number), or anyletter(number)
    regex = r"[A-Za-z\d]+\((\d+)\)"
    sort!(v, by=x -> eval(Meta.parse(match(regex, x[2])[1])))
    # Sort by resonance type: L (Lambda), S (Sigma), D (Delta), K (Kaon)
    sort!(v, by=x -> findfirst(x[2][1], "LSDK"))
    order = getindex.(v, 1)

    mypairs = isobarnames[order] .=> zip(couplings[order], chains[order])
    x = Vector{Pair{String,Tuple{Complex,AbstractDecayChain}}}(mypairs)
    ThreeBodyDecay(x)
end
