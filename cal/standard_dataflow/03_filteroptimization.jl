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
using Plots: scatter
using Measures
using Optim
using BSplineKit
using Printf
<<<<<<< HEAD:cal/standard_dataflow/03_filteroptimization.jl
=======
# set data configuration (where to find data; and where to save results)
# if gethostname() == "Lisas-MacBook-Pro.local"
#     ENV["LEGEND_DATA_CONFIG"] = "/Users/lisa/Documents/Workspace/LEGEND/LBL_ASIC/ASIC_data/ppc01/config.json"
# else # on NERSC 
    ENV["LEGEND_DATA_CONFIG"] = "/global/cfs/projectdirs/m2676/data/teststands/lbnl/PPC01_ryutaro/config.json"
# end 
>>>>>>> bd32341c0f9bde80f2aff9a94eedd2234c34bb9c:cal/p01/03_filteroptimization.jl

# include relevant functions 
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/src/filteropt_rt_optimization_blnoise.jl")
include("$(@__DIR__)/$relPath/utils/utils_plot.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")
include("$(@__DIR__)/$relPath/processing_funcs/process_filteropt_rt.jl")

# inputs 

asic = LegendData(:ppc01)
<<<<<<< HEAD:cal/standard_dataflow/03_filteroptimization.jl
period = DataPeriod(3)
run = DataRun(35)
=======
period = DataPeriod(1)
run = DataRun(3)
>>>>>>> bd32341c0f9bde80f2aff9a94eedd2234c34bb9c:cal/p01/03_filteroptimization.jl
channel = ChannelId(1)
category = DataCategory(:cal)
filter_types = [:trap]#, :cusp]

Plots_theme()
# load configs and modify if needed 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
dsp_config = DSPConfig(dataprod_config(asic).dsp(filekeys[1]).default)
pz_config = dataprod_config(asic).dsp(filekeys[1]).pz.default
peak =  Symbol(pz_config.peak)
τ_pz = mvalue(get_values(asic.par.rpars.pz[period, run, channel]).τ)

# do optimization 
process_filteropt_rt(asic, period, run, category, channel, dsp_config, τ_pz, peak; reprocess = true, rt_opt_mode = :bl_noise, filter_types = filter_types)
println("Filter optimization parameters: ", asic.par.rpars.fltopt[period, run, channel])

# read filter optimization pars
fltopt_pars = asic.par.rpars.fltopt[period,run,channel]
