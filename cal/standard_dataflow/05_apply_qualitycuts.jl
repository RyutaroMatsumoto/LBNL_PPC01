# Purpose: Apply quality cuts to the data.
# 1. load quality cut (qc) configuration. this is where cut thresholds are defined.
# 2. apply quality cuts to the data.
# 3. save quality cut results (overall cut efficiency and flag (cut/not cut) for each waveform ) to rpars 
using LegendDataManagement
using LegendDataManagement: readlprops, writelprops
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using PropDicts
using Unitful
using TypedTables

# include relevant functions 
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/src/apply_qc.jl")
include("$(@__DIR__)/$relPath/processing_funcs/process_qualitycuts.jl")

# inputs
asic = LegendData(:ppc01)
period = DataPeriod(3)
run = DataRun(1)
channel = ChannelId(1)
category = :cal 

# load qc config and apply 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
qc_config = dataprod_config(asic).qc(filekeys[1]).default
process_qualitycuts(asic, period, run, category, channel; reprocess = true, qc_config = qc_config);

# read quality cuts from pars 
qc = asic.par.rpars.qc[period][run][channel]
qc.wvf_keep.all

sum(qc.wvf_keep.all)
qc.qc_surv.all