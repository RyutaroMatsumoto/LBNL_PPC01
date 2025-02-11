using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendDSP
using LegendDSP: get_fltpars
using LegendSpecFits: get_friedman_diaconis_bin_width
using LegendHDF5IO
using RadiationSpectra
using RadiationDetectorDSP
using Measurements: value as mvalue
using PropDicts
using StatsBase, IntervalSets
using Unitful
using TypedTables
using CairoMakie
using Measures

# include relevant functions 
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/processing_funcs/process_peak_split.jl")
include("$(@__DIR__)/$relPath/src/simple_dsp.jl")
include("$(@__DIR__)/$relPath/src/apply_qc.jl")
include("$(@__DIR__)/$relPath/utils/utils_plot.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")

# inputs
reprocess = true
asic = LegendData(:ppc01)
period = DataPeriod(1)
run = DataRun(13)
channel = ChannelId(1)
category = DataCategory(:cal)
source = :co60

# load configs 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
qc_config = dataprod_config(asic).qc(filekeys[1]).default
ecal_config = dataprod_config(asic).energy(filekeys[1]).default
dsp_config = DSPConfig(dataprod_config(asic).dsp(filekeys[1]).default)

# run processor 
plts = process_peak_split(asic, period, run, category, channel, ecal_config, dsp_config, qc_config ; reprocess = reprocess)   

# plots 
display(plts[1])
display(plts[2])
display(plts[3])

# read peakfiles 
data_peak  = read_ldata(:Co60a, asic, :jlpeaks, category, period, run, channel)
data_peak.waveform