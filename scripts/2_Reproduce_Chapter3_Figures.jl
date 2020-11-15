using DrWatson
quickactivate(@__DIR__, "Chemostat_Kayser2005")

using Chemostat_Dynamics
using Chemostat_Dynamics.Polytopes
using Chemostat_Dynamics.Utilities
using Chemostat_Dynamics.MonteCarlo
using Chemostat_Dynamics.MaxEnt
using Plots
using Statistics

## ------------------------------------------------------------------