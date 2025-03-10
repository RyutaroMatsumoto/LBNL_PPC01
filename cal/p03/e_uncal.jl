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




relPath_processor = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath_processor/utils/utils_aux.jl")
include("$(@__DIR__)/$relPath_processor/utils/utils_plot.jl")


reprocess = true 
asic = LegendData(:ppc01)
period = DataPeriod(3)
run = DataRun(13)
channel = ChannelId(1)
category = :cal 
e_types = [:e_trap]#, :e_cusp]

qc = asic.par.rpars.qc[period][run][channel]

dsp_pars = Table(read_ldata(asic, :jldsp, category, period, run, channel);)
e_trap_clean = filter(!isnan, dsp_pars.e_trap) # Remove NaN values
stephist(e_trap_clean,bins = 10000)
vline!([quantile(e_trap_clean, 0.875)])