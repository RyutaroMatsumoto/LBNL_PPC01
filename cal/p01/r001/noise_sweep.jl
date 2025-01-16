using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using LegendSpecFits
using LegendDSP
using LegendDSP: get_fltpars
using RadiationDetectorSignals
using RadiationDetectorDSP
using Measurements: value as mvalue
using PropDicts
using StatsBase, IntervalSets
using Unitful
using TypedTables
using Plots
using Measures
using Optim
using BSplineKit
using Printf

# include relevant functions 
include("$(@__DIR__)/../../../../../src/filteropt_rt_optimization_blnoise.jl")
include("$(@__DIR__)/../../../../../utils/utils_plot.jl")
include("$(@__DIR__)/../../../../../utils/utils_aux.jl")
include("$(@__DIR__)/../../../../../processing_funcs/process_filteropt.jl")

ENV["LEGEND_DATA_CONFIG"] = "/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json"

# inputs
asic = LegendData(:ppc01)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)

filter_type = :cusp 

# read waveforms 
peak = :Co60a
data_peak  = read_ldata((peak), asic, :jlpeaks, category, period, run, channel)
wvfs = data_peak.waveform

# read filter optimization pars
fltopt_pars = asic.par.rpars.fltopt[period,run,channel][filter_type]

# read mean decay times 
τ_pz = mvalue(get_values(asic.par.rpars.pz[period, run, channel]).τ)

# collect inputs for noise sweep and overwrite default settings in dsp_config
ft = mvalue(fltopt_pars.ft)
rt = 1.5u"µs" :0.5u"µs" :6.0u"µs" 
dsp_config = asic.metadata.config.dsp.dsp_config.default
dsp_config[Symbol("e_grid_$(filter_type)")].rt.start = first(rt)
dsp_config[Symbol("e_grid_$(filter_type)")].rt.step = step(rt)
dsp_config[Symbol("e_grid_$(filter_type)")].rt.stop = last(rt)

# do noise sweep
result, report  = filteropt_rt_optimization_blnoise(filter_type, wvfs, DSPConfig(dsp_config), τ_pz; ft = ft) 


using Plots
scatter(ustrip.(collect(report.enc_grid_rt)), report.enc)

