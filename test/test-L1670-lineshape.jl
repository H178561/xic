using Plots
using QuadGK
using Parameters

using Lc2ppiKSemileptonicModelLHCb

theme(:boxed)

struct FlatteBelle1670
    mf::Float64
    Γf::Float64
    ma2::Float64
    mb2::Float64
    g2::Float64
end

k(m, ma, mb) = Lc2ppiKSemileptonicModelLHCb.breakup(m^2, ma^2, mb^2)
function (d::FlatteBelle1670)(σ::Float64)
    iϵ = 1e-6im
    m = sqrt(σ)
    D = m - d.mf + 0.5im * (d.Γf + d.g2 * k(m + iϵ, d.ma2, d.mb2))
    return 1 / D
end

d = FlatteBelle1670(1.6744, 0.0272, 1.115683, 0.547862, 0.258)


m_range = range(1.55, 1.8, length=1000)
plot(m_range, m -> m * abs2(d(m^2));
    xlab="m [GeV]", ylab="|A|²", label="Flatte Belle L(1670)")

# 
struct FlatteNR
    mf::Float64
    g1::Float64
    ma1::Float64
    mb1::Float64
    g2::Float64
    ma2::Float64
    mb2::Float64
end
function (d::FlatteNR)(σ::Float64)
    iϵ = 1e-6im
    m = sqrt(σ)
    D = m - d.mf + 0.5im * (d.g1 * k(m + iϵ, d.ma1, d.mb1) + d.g2 * k(m + iϵ, d.ma2, d.mb2))
    return 1 / D
end

dNR = let
    m = 1.6744
    Γ1 = 0.0272
    ma1, mb1 = 0.938272046, 0.493677
    k1 = k(m, ma1, mb1)
    g1 = Γ1 / k1
    @show g1
    #
    FlatteNR(m,
        g1, ma1, mb1,
        0.258, 1.115683, 0.547862)
end


plot!(m_range, m -> m * abs2(dNR(m^2));
    xlab="m [GeV]", ylab="|A|²", label="Flatte NR L(1670)")
