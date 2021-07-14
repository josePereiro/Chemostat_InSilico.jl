using Plots: length
@time begin 
   import Chemostat_InSilico
   const Dyn = Chemostat_InSilico.Dynamic
   import Chemostat_InSilico.Dynamic: 
      SimD2, run_simD2!, check_stst, hist, 
      n_ave_conv!, tail_ave,
      plotsdir, lglob, sglob, Container, vec!, Vi, 
      normalize!, save_sim,
      set_status, get_status,
      load_batch, save_batch,
      UNDONE_SIM_STATUS, DEAD_SIM_STATUS,
      EXPLODED_SIM_STATUS, NITERS_SIM_STATUS, 
      STST_SIM_STATUS, sim_file, collect_ts
       
   import UtilsJL
   const PltU = UtilsJL.PlotsUtils
   const Ass = UtilsJL.ProjAssistant
   const GU = UtilsJL.GeneralUtils

   using Plots
   using Base.Threads
   using ProgressMeter
   using BenchmarkTools
   using ExtractMacro

   import DataFileNames
   const DFN = DataFileNames

   Ass.set_verbose(false)
end

## ------------------------------------------------------
# batch = (; push_frec, Xts, sgts, z_avts, ug_avts, cgD_Xts, Pzts, Pugts)
# ------------------------------------------------------
Ds, ϵs, cgs, simid = Dyn.lglob(:Ds, :ϵs, :cgs, :SimD2Id)

## ------------------------------------------------------

# ## ------------------------------------------------------
# let

#    tsid = [:Xts, :sgts, :z_avts, :ug_avts, :cgD_Xts, :Pzts, :Pugts]
#    tsid = [:Xts, :sgts, :z_avts, :ug_avts, :cgD_Xts] |> rand
   
#    Ds, ϵs, cgs, simid = Dyn.lglob(:Ds, :ϵs, :cgs, :SimD2Id)
#    iter = Iterators.product(Ds, ϵs, cgs)
#    # for (simi, (D, ϵ, cg)) in enumerate(iter)
#    D = 0.5
#    ϵ = 0.8585714285714285
#    cg = Inf
#       (t, s) = collect_ts(tsid, simid, (;D, ϵ, cg))

#       p = plot(;xlabel = "time", ylabel = string(tsid))
#       plot!(p, t, s; label = "", color = :blue)
#       return p
#    # end
# end



# ## ------------------------------------------------------
# function plot_sim_summary()

   
# end
# ## ------------------------------------------------------

#     # # plots
#     # # post process time series
#     # n = 100
#     # Xts = n_ave_conv!(vec!(Xts), n)
#     # z_avts = n_ave_conv!(vec!(z_avts), n)
#     # ug_avts = n_ave_conv!(vec!(ug_avts), n)
#     # uo_avts = n_ave_conv!(vec!(uo_avts), n)
#     # cgD_Xts = n_ave_conv!(vec!(cgD_Xts), n)
#     # sgts = n_ave_conv!(vec!(sgts), n)

#     # Pzts = vec!(Pzts)
#     # Pugts = vec!(Pugts)
#     # Puots = vec!(Puots)

#     # lock(lk) do

#     #     # plots
#     #     prms = (;title = string(status), label = "", xrotation = 45)

#     #     # z
#     #     pz = bar(hist(Vz, P); xlabel = "z", ylabel = "proj", prms...)
#     #     vline!(pz, [S.z_av]; color = :blue, prms...)
#     #     vline!(pz, [S.D]; color = :red, prms...)
        
#     #     # ug
#     #     pug = bar(hist(Vug, P); xlabel = "ug", ylabel = "proj", prms...)
#     #     vline!(pug, [S.ug_av]; color = :blue, prms...)
#     #     vline!(pug, [S.cgD_X]; color = :red, prms...)

#     #     # uo
#     #     puo = bar(hist(Vuo, P); xlabel = "uo", ylabel = "proj", prms...)
#     #     vline!(puo, [S.uo_av]; color = :blue, prms...)

#     #     # ug
#     #     pug_av = plot(ug_avts; color = :blue, ylabel = "ug_av & cgD_Xts", prms...)
#     #     plot!(pug_av, cgD_Xts; color = :red, prms...)
        
#     #     # z
#     #     pz_av = plot(z_avts; ylabel = "z_av", prms...)
#     #     pz_av = hline!(pz_av, [S.D]; color = :red, prms...)
        
#     #     # X
#     #     pX = plot(Xts; ylabel = "X", prms...)
#     #     psg = plot(sgts; ylabel = "sg", prms...)
        
#     #     PltU.sfig(
#     #         [pz, pug, puo, pug_av, pX, pz_av, psg],
#     #         [plotsdir()], "sim3D", simparams, ".png"
#     #     )

#     #     dat = (;S, status, ug_avts, cgD_Xts, z_avts, Xts, sgts, Pzts)
#     #     Dyn.sdat(dat, simdfile; verbose = false)

#     # end # lock(lk) do

#     ## ------------------------------------------------------
#  # # gifs
#     # gifs_ = map(zip(Pzts, Pugts)) do hists
#     #     Pz_, Pug_ = hists
#     #     pz_ = bar(Pz_; xlabel = "z", ylabel = "proj", prms...)
#     #     pug_ = bar(Pug_; xlabel = "ug", ylabel = "proj", prms...)
#     #     PltU.sfig([pz_, pug_], tempname(), ".png")
#     # end
#     # @show typeof(gifs_)
#     # @show length(gifs_)
    