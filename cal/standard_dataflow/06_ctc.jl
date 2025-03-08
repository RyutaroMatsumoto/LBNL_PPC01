using LegendDataManagement
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using PropDicts
using LegendSpecFits
using StatsBase
using Makie, LegendMakie, CairoMakie
using Unitful
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")
include("$(@__DIR__)/$relPath/processing_funcs/process_ctc.jl")

asic = LegendData(:ppc01)
period = DataPeriod(3)
run = DataRun(2)
channel = ChannelId(1)
det = _channel2detector(asic, channel)
category = :cal 

# load configs 
reprocess = true
filekey = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])[1]
ecal_config = dataprod_config(asic).energy(filekey).default
ctc_config = dataprod_config(asic).energy(filekey).ctc.default

# run processor
result = process_ctc(data, period, run, category, channel, ecal_config, ctc_config; reprocess=reprocess, juleana_logo=false)

# check if result has been properly saved
asic.par.rpars.ctc[period, run, channel].e_trap

