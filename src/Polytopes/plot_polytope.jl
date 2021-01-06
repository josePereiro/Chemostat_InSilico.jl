function plot_polytope!(plt, p::Polytope)
    
    plot!(plt, xlabel = "vatp", ylabel = "vg")
    #locals
    kwargs = (;color = :black, lw = 3, label = "")
    lb, ub = vatpL(p), vatpU(p)
    plot!(plt, (vg) -> vgU(vg, p), lb, ub; kwargs...)
    plot!(plt, (vg) -> vgL(vg, p), lb, ub; kwargs...)

    # globals
    kwargs = (;color = :black, lw = 4, ls = :dot, label = "")
    hline!(plt, [vgU(p)]; kwargs...)
    hline!(plt, [vgL(p)]; kwargs...)
    vline!(plt, [vatpL(p)]; kwargs...)
    vline!(plt, [vatpU(p)]; kwargs...)

    # margins
    xmargin = abs(vatpU(p) - vatpL(p)) * 0.1
    ymargin = abs(vgU(p) - vgL(p)) * 0.1
    plot!(xlim = [vatpL(p) - xmargin, vatpU(p) + xmargin])
    plot!(ylim = [vgL(p) - ymargin, vgU(p) + ymargin])
    plt
end
plot_polytope(p::Polytope) = plot_polytope!(plot(), p)