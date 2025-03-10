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
using Plots 
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
run = DataRun(31)
channel = ChannelId(1)
category = :cal 
e_types = [:e_trap]#, :e_cusp]

# do processing 
process_hit(asic, period, run, category, channel; reprocess = reprocess, e_types = e_types)

# read hit files and plot -> sanity check 
hit_par = Table(read_ldata(asic, :jlhit, category, period, run))

# 
Plots_theme()
bins = 1000:0.2:2500
Plots.stephist(hit_par.e_trap, bins = bins, fill = true, color = :silver,label = "Before qc", ylims = [1, 100], yscale = :log10)
Plots.stephist!(hit_par.e_trap[hit_par.qc], bins = bins, fill = true, alpha = 0.5, color = :violet, label = "After qc", xlabel = "Calibrated energy", ylabel = "Counts", legend = :topleft)

