# Purpose: Create hitfiles for a given run and channel. They contain calibrated energies and quality cut results.
# Useful for further analysis of calibrated energy spectrum or plotting. 
using LegendDataManagement
using LegendDataManagement: readlprops, writelprops
using LegendDataManagement.LDMUtils
using LegendSpecFits
using LegendSpecFits: get_friedman_diaconis_bin_width 
using LegendHDF5IO
using RadiationSpectra
using PropDicts
using Unitful
using TypedTables
using Statistics, StatsBase
using IntervalSets
using Makie, LegendMakie, CairoMakie 
using Unitful, Measures
using Measurements: value as mvalue

# include relevant functions 
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")
include("$(@__DIR__)/$relPath/utils/utils_plot.jl")
include("$(@__DIR__)/$relPath/processing_funcs/process_hit.jl")

# inputs
reprocess = true
asic = LegendData(:ppc01)
period = DataPeriod(3)
run = DataRun(39)
channel = ChannelId(1)
category = :cal 
e_types = [:e_trap]#, :e_cusp]

# do processing 
process_hit(asic, period, run, category, channel; reprocess = reprocess, e_types = e_types)

asic.par.rpars.ecal[period, run, channel].e_trap.cal.func

# read hit files and plot -> sanity check 
hit_par = Table(read_ldata(asic, :jlhit, category, period, run))

# sanity plot 
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Calibrated energy", ylabel = "Counts", yscale = log10, limits = ((nothing, nothing), (1, nothing)))
bins = 300:5:1500
hist!(ax, ustrip.(filter(isfinite, hit_par.e_trap)), bins = bins,  color = :silver, label = "Before qc")
hist!(ax, ustrip.(filter(isfinite, hit_par.e_trap[hit_par.qc])), bins = bins, color = :violet, label = "After qc")
axislegend(ax, position = :lt)
fig