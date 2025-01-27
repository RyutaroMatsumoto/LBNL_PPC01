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
# set data configuration (where to find data; and where to save results)
if gethostname() == "Lisas-MacBook-Pro.local"
    ENV["LEGEND_DATA_CONFIG"] = "/Users/lisa/Documents/Workspace/LEGEND/LBL_ASIC/ASIC_data/ppc01/config.json"
else # on NERSC 
    ENV["LEGEND_DATA_CONFIG"] = "/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json"
end 

# include relevant functions 
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/src/filteropt_rt_optimization_blnoise.jl")
include("$(@__DIR__)/$relPath/utils/utils_plot.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")
include("$(@__DIR__)/$relPath/processing_funcs/process_filteropt.jl")

# inputs 
asic = LegendData(:ppc01)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)
filter_types = [:trap, :cusp]

Plots_theme()
# load configs and modify if needed 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
dsp_config = DSPConfig(dataprod_config(asic).dsp(filekeys[1]).default)
pz_config = dataprod_config(asic).dsp(filekeys[1]).pz.default
peak =  Symbol(pz_config.peak)
τ_pz = mvalue(get_values(asic.par.rpars.pz[period, run, channel]).τ)

# do optimization 
process_filteropt(asic, period, run, category, channel, dsp_config, τ_pz, peak; reprocess = true, rt_opt_mode = :bl_noise, filter_types = filter_types)

# read filter optimization pars
fltopt_pars = asic.par.rpars.fltopt[period,run,channel]
