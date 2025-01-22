# Extract decay time of exponential tail of HPGe waveforms using peak files 
# Perform fit (truncated gaussian) to get average decay time, that is later used for pole-zero correction of all waveforms
# Note: Use peak files instead of all waveforms, because they are "good" waveforms 
env_path = "$(@__DIR__)/" * relpath(split(@__DIR__, "LBNL_PPC01")[1], @__DIR__)  * "/LBNL_PPC01/"
import Pkg; Pkg.activate(env_path)
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendDSP
using LegendDSP: get_fltpars
using LegendHDF5IO
using LegendSpecFits
using RadiationSpectra
using RadiationDetectorDSP
using Measurements: value as mvalue
using PropDicts
using StatsBase, IntervalSets
using Unitful
using TypedTables
using Plots
using Measures

# set data configuration (where to find data; and where to save results)
ENV["LEGEND_DATA_CONFIG"] = "/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json"

# load functions from hpge-ana
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/processing_funcs/process_decaytime.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")

# inputs
asic = LegendData(:ppc01)
period = DataPeriod(2)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)

pz_config = asic.metadata.config.dsp.dsp_config.pz.default

pz_config.min_tau = 0.0u"µs"
pz_config.max_tau = 650.0u"µs"
pz_config.rel_cut_fit = 0.01
process_decaytime(asic, period, run, category, channel; reprocess = true, pz_config = pz_config)

# sanity check: read decay time pars
pars_pz = asic.par.rpars.pz[period, run, channel]




