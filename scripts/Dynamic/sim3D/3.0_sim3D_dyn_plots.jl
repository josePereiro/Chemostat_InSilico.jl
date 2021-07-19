using Plots: length
@time begin 
   import Chemostat_InSilico
   const Dyn = Chemostat_InSilico.Dynamic
   import Chemostat_InSilico.Dynamic: 
      SimD3, run_simD3!, check_stst, hist, 
      n_ave_conv!, tail_ave,
      plotsdir, lglob, sglob, Container, vec!, Vi, 
      normalizeP!, 
      save_simdat, load_simdat, simdat_file,
      set_status, get_status,
      load_batch, save_batch,
      UNDONE_SIM_STATUS, DEAD_SIM_STATUS,
      EXPLODED_SIM_STATUS, NITERS_SIM_STATUS, 
      STST_SIM_STATUS, collect_ts
       
   import UtilsJL
   const PltU = UtilsJL.PlotsUtils
   const Ass = UtilsJL.ProjAssistant
   const GU = UtilsJL.GeneralUtils

   using Plots
   import GR
   !isinteractive() && GR.inline("png")

   using Base.Threads
   using ProgressMeter
   using BenchmarkTools
   using ExtractMacro

   import DataFileNames
   const DFN = DataFileNames

   Ass.set_verbose(false)
end

# ------------------------------------------------------
# batch = (; push_frec, Xts, sgts, z_avts, ug_avts, uo_avts, cgD_Xts, Pzts, Pugts, Puots)
Ds, ϵs, cgs, simid = Dyn.lglob(:Ds, :ϵs, :cgs, :SimD3Id)

## ------------------------------------------------------
include("3.1_sim3D_evolution.jl")

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
    