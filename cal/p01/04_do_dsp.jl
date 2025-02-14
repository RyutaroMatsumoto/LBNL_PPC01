using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using Plots 

# include relevant functions 
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/src/simple_dsp.jl")
include("$(@__DIR__)/$relPath/processing_funcs/process_dsp.jl")

# inputs 
asic = LegendData(:ppc01)
period = DataPeriod(3)
run = DataRun(2)
channel = ChannelId(1)
category = DataCategory(:cal)

# load configs 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
dsp_config = DSPConfig(dataprod_config(asic).dsp(filekeys[1]).default)
τ_pz = mvalue(get_values(asic.par.rpars.pz[period, run, channel]).τ)
pars_filter = asic.par.rpars.fltopt[period,run,channel]
# pars_filter = PropDict() # if empty, use default filter parameters 
# do dsp
process_dsp(asic, period, run, category, channel, dsp_config, τ_pz, pars_filter; reprocess = true)

# read dsp pars and plot 
dsp_pars = read_ldata(asic, :jldsp, category, period, run, channel);
Table(dsp_pars)
columnnames(Table(dsp_pars))
stephist(dsp_pars.e_trap)
