
breakup(m², m1², m2²) = sqrtKallenFact(sqrt(m²), sqrt(m1²), sqrt(m2²)) / (2 * sqrt(m²))

function F²(l, p, p0, d)
    pR = p * d
    p0R = p0 * d
    l == 0 && return 1.0
    l == 1 && return (1 + p0R^2) / (1 + pR^2)
    l == 2 && return (9 + 3p0R^2 + p0R^4) / (9 + 3pR^2 + pR^4)
    l != 3 && error("l==$(l)>2 cannot be")
    return (225 + 45 * (p0R^2) + 6 * (p0R^2)^2 + (p0R^2)^3) / (225 + 45 * (pR^2) + 6 * (pR^2)^2 + (pR^2)^3)
end


abstract type Lineshape end

@with_kw struct BreitWignerMinL{T} <: Lineshape
    pars::T
    l::Int
    minL::Int
    #
    name::String
    #
    m1::Float64
    m2::Float64
    mk::Float64
    m0::Float64
end
BreitWignerMinL(pars::T; kw...) where {T} = BreitWignerMinL(; pars, kw...)
function (BW::BreitWignerMinL)(σ)
    dR, dΛc = 1.5, 5.0 # /GeV
    m, Γ₀ = BW.pars
    @unpack l, minL = BW
    @unpack m1, m2, mk, m0 = BW
    p, p0 = breakup(σ, m1^2, m2^2), breakup(m^2, m1^2, m2^2)
    q, q0 = breakup(m0^2, σ, mk^2), breakup(m0^2, m^2, mk^2)
    Γ = Γ₀ * (p / p0)^(2l + 1) * m / sqrt(σ) * F²(l, p, p0, dR)
    1 / (m^2 - σ - 1im * m * Γ) * (p / p0)^l * (q / q0)^minL *
    sqrt(F²(l, p, p0, dR) * F²(minL, q, q0, dΛc))
end

# BuggBreitWignerMinL
@with_kw struct BuggBreitWignerMinL{T} <: Lineshape
    pars::T
    l::Int
    minL::Int
    #
    name::String
    #
    m1::Float64
    m2::Float64
    mk::Float64
    m0::Float64
end
BuggBreitWignerMinL(pars::T; kw...) where {
    T <: NamedTuple{X, Tuple{Float64, Float64}}} where {X} =
    BuggBreitWignerMinL(; pars = merge(pars, (γ = 1.1,)), kw...)
#
function (BW::BuggBreitWignerMinL)(σ)
    σA = mK^2 - mπ^2 / 2
    m, Γ₀, γ = BW.pars
    @unpack m1, m2 = BW
    Γ = (σ - σA) / (m^2 - σA) * Γ₀ * exp(-γ * σ)# * breakup(σ,m1^2,m2^2)/(2*sqrt(σ))
    1 / (m^2 - σ - 1im * m * Γ)
end

# Flatte1405
@with_kw struct Flatte1405{T} <: Lineshape
    pars::T
    l::Int
    minL::Int
    #
    name::String
    #
    m1::Float64
    m2::Float64
    mk::Float64
    m0::Float64
end
#
Flatte1405(pars::T; kw...) where {T} = Flatte1405(; pars, kw...)
function (BW::Flatte1405)(σ)
    m, Γ₀ = BW.pars
    @unpack m1, m2, m0, mk = BW
    p, p0 = breakup(σ, m1^2, m2^2), breakup(m^2, mπ^2, mΣ^2)
    p′, p0′ = breakup(σ, mπ^2, mΣ^2), breakup(m^2, mπ^2, mΣ^2)
    Γ1 = Γ₀ * (p / p0) * m / sqrt(σ)
    Γ2 = Γ₀ * (p′ / p0′) * m / sqrt(σ)
    Γ = Γ1 + Γ2
    1 / (m^2 - σ - 1im * m * Γ)
end

@with_kw struct TFCode1670{T} <: Lineshape
    pars::T
    l::Int
    minL::Int
    #
    name::String
    #
    m1::Float64
    m2::Float64
    mk::Float64
    m0::Float64
end
TFCode1670(pars::T; kw...) where {T} = TFCode1670(; pars, kw...)

k(m, ma, mb) = breakup(m^2, ma^2, mb^2)
function (d::TFCode1670)(σ::Float64)
    mf, _ = d.pars
    # 
    # the value of the coupling comes from belle width / k to match G1 (const) for nominal mass
    g1 = 0.0272
    g2, ma2, mb2 = 0.258, 1.115683, 0.547862

    iϵ = 1e-8im
    m = sqrt(σ)
    # the line does not look natural, but it's intentional
    D = mf^2 - m^2 - 1im * (g1 + g2^2 * 2k(m + iϵ, ma2, mb2) / m)
    return 1 / D
end


@with_kw struct L1670Flatte{T} <: Lineshape
    pars::T
    l::Int
    minL::Int
    #
    name::String
    #
    m1::Float64
    m2::Float64
    mk::Float64
    m0::Float64
end
L1670Flatte(pars::T; kw...) where {T} = L1670Flatte(; pars, kw...)

k(m, ma, mb) = breakup(m^2, ma^2, mb^2)
function (d::L1670Flatte)(σ::Float64)
    mf, _ = d.pars
    # 
    # the value of the coupling comes from belle width / k to match G1 (const) for nominal mass
    Γ0 = 0.0272
    g, ma2, mb2 = 0.258, 1.115683, 0.547862

    iϵ = 1e-8im
    m = sqrt(σ)
    # the line does not look natural, but it's intentional
    D = mf - m - 0.5im * (Γ0 + g * k(m + iϵ, ma2, mb2))
    return 1 / D
end



function updatepars(BW::Lineshape, pars)
    fiels = fieldnames(typeof(BW))
    values = [getproperty(BW, f) for f in fiels]
    return typeof(BW)(; NamedTuple{fiels}(values)..., pars)
end



@recipe function f(BW::Lineshape)
    xv = range((BW.m1 + BW.m2)^2, (BW.m0 - BW.mk)^2, length = 300)
    intensity(σ) = abs2(BW(σ)) *
                   breakup(σ, BW.m1^2, BW.m2^2) *
                   breakup(BW.m0^2, σ, BW.mk^2) / sqrt(σ)
    yv = intensity.(xv)
    (xv, yv ./ sum(yv) .* length(yv))
end


