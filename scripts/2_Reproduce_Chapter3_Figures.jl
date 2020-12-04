using DrWatson
quickactivate(@__DIR__, "Chemostat_Kayser2005")

using Chemostat_Dynamics
using Chemostat_Dynamics.Polytopes
using Chemostat_Dynamics.Utilities
using Chemostat_Dynamics.MonteCarlo
using Chemostat_Dynamics.MaxEnt
using Plots
using Statistics
using Base.Threads

## ------------------------------------------------------------------
ps = @time let 

    # Time
    niters = Int(50)
    
    # Params
    ϵ = 0.0
    ϕ = 1.0
    # Initialization
    Δpi = ΔPolytope()
    n = nbins(Δpi)
    Xi = 1e-5
    Di = 0.0001
    setX!((vatp, vg) -> Xi/n, Δpi)
    sgi = Δpi.p.cg
    sli = 0.0
    ps = []
    for t in 1:niters
        Xi = sum(values(Δpi.Xs)) 
        @info "iter: $t" Xi sgi sli
        
        # collect
        nϵ_z_Σvg_X_div_vg_nbins = Dict{Float64, Float64}()
        ϵ_Σvgvatp_zX_div_nbins = 0.0
        for vatp in vatp_range(Δpi)
            Σvg_Xvatpvg = sum(getX(Δpi, vatp, vg) for vg in vg_range(Δpi, vatp))
            nϵ_z_Σvg_X_div_vg_nbins[vatp] = (1 - ϵ) * z(vatp, Δpi) * Σvg_Xvatpvg / vg_nbins(Δpi, vatp)
            ϵ_Σvgvatp_zX_div_nbins += nϵ_z_Σvg_X_div_vg_nbins[vatp]
        end
        ϵ_Σvgvatp_zX_div_nbins /= nbins(Δpi)
        ϵ_Σvgvatp_zX_div_nbins *= ϵ

        # update
        Δsg, Δsl = 0.0, 0.0
        for vatp in vatp_range(Δpi)
            ΔXvatpvg = nϵ_z_Σvg_X_div_vg_nbins[vatp]
            for vg in vg_range(Δpi, vatp)
                X = getX(Δpi, vatp, vg)
                ΔXvatpvg += ϵ_Σvgvatp_zX_div_nbins - ϕ*Di*X
                setX!(Δpi, vatp, vg, max(0.0, X + ΔXvatpvg))
                Δsg -= vg*X
                Δsl -=  vg*X
            end
        end
        Δsg += Di*(Δpi.p.cg - sgi)
        Δsl += Di*(Δpi.p.cl - sli)
        sgi = max(sgi - Δsg, 0.0)
        sli = max(sli - Δsl, 0.0)

        # Init polytope
        push!(ps, Δpi)
        Δpi = addXs!(ΔPolytope(;cg = sgi, cl = sli, X0 = 0, D = Di), Δpi)
    end
    ps
end;

## ------------------------------------------------------------------
let
    p = plot()
    for pol in ps
        plot_polytope!(p, pol)
    end
    p
end

## ------------------------------------------------------------------
# function sg(p::XPolytope, sg)
#     vgX = sum(vg * p.Xs[(vatp, vg)] for vatp in p.vatp_range, vg in p.vg_range)
#     return max(0, -vgX + p.pol.D * (p.pol.cg - sg))
# end

# function sl(p::Polytope, sl)
#     vlX = sum(vl(vatp, vg, p) * p.Xs[(vatp, vg)] for vatp in p.vatp_range, vg in p.vg_range)
#     return max(0, -vlX + p.D*(p.pol.cl - sl))
# end