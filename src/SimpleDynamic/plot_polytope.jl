function plot_polytop(p)
    plt = plot(xlabel = "vatp", ylabel = "vg")
    #locals
    kwargs = (;color = :black, lw = 3, label = "")
    lb, ub = vatp_global_min(p), vatp_global_max(p)
    plot!(plt, (vg) -> vg_local_max(vg, p), lb, ub; kwargs...)
    plot!(plt, (vg) -> vg_local_min(vg, p), lb, ub; kwargs...)

    # globals
    kwargs = (;color = :black, lw = 4, ls = :dot, label = "")
    hline!(plt, [vg_global_max(p)]; kwargs...)
    hline!(plt, [vg_global_min(p)]; kwargs...)
    vline!(plt, [vatp_global_min(p)]; kwargs...)
    vline!(plt, [vatp_global_max(p)]; kwargs...)

    # margins
    xmargin = abs(vatp_global_max(p) - vatp_global_min(p)) * 0.1
    ymargin = abs(vg_global_max(p) - vg_global_min(p)) * 0.1
    plot!(xlim = [vatp_global_min(p) - xmargin, vatp_global_max(p) + xmargin])
    plot!(ylim = [vg_global_min(p) - ymargin, vg_global_max(p) + ymargin])
    plt
end