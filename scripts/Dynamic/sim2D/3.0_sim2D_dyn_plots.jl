@time begin 
   import Chemostat_InSilico
   const Dyn = Chemostat_InSilico.Dynamic
   import Chemostat_InSilico.Dynamic: 
      SimD2, run_simD2!, check_stst, hist, 
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
# batch = (; push_frec, Xts, sgts, z_avts, ug_avts, cgD_Xts, Pzts, Pugts)
Ds, ϵs, cgs, simid = Dyn.lglob(:Ds, :ϵs, :cgs, :SimD2Id)

## ------------------------------------------------------
# include("3.1_sim2D_evolution.jl")