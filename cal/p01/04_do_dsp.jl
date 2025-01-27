using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using Plots 

# set data configuration (where to find data; and where to save results)
if gethostname() == "Lisas-MacBook-Pro.local"
    ENV["LEGEND_DATA_CONFIG"] = "/Users/lisa/Documents/Workspace/LEGEND/LBL_ASIC/ASIC_data/ppc01/config.json"
else # on NERSC 
    ENV["LEGEND_DATA_CONFIG"] = "/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json"
end 

# include relevant functions 
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/src/simple_dsp.jl")
include("$(@__DIR__)/$relPath/processing_funcs/process_dsp.jl")

# inputs 
asic = LegendData(:ppc01)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)

# load configs 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
dsp_config = DSPConfig(dataprod_config(asic).dsp(filekeys[1]).default)
τ_pz = mvalue(get_values(asic.par.rpars.pz[period, run, channel]).τ)
pars_filter = asic.par.rpars.fltopt[period,run,channel]

# do dsp
process_dsp(asic, period, run, category, channel, dsp_config, τ_pz, pars_filter; reprocess = false)

# read dsp pars and plot 
dsp_pars = read_ldata(asic, :jldsp, category, period, run, channel);
Table(dsp_pars)
columnnames(Table(dsp_pars))
stephist(dsp_pars.e_trap)
