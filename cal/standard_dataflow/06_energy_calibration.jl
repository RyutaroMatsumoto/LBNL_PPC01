# Purpose: Perform energy calibration for a given channel and period/run.
# 1. load energy calibration (ecal) configuration. this is where calibration parameters, such as gamma lines or energy windows, are defined.
# 2. perform "simple calibration" - peak search on uncalibrated data
# 3. perform actual energy calibration - fit peaks to known gamma lines and extract calibration parameters
# 4. fit calibration curve
# 5. fit resolution curve
# 6. save calibration parameters to rpars
# several sanity plots document this procedure and are automatically saved in jlplt tier 
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
using LinearAlgebra
# using Plots 
using Makie, LegendMakie, CairoMakie
using Unitful, Measures
using Measurements: value as mvalue

# include relevant functions 
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")

include("$(@__DIR__)/$relPath/utils/utils_plot.jl")
include("$(@__DIR__)/$relPath/processing_funcs/process_energy_calibration.jl")

# inputs
reprocess = true 
asic = LegendData(:ppc01)
period = DataPeriod(3)
run = DataRun(31)
channel = ChannelId(1)
category = :cal 
e_types = [:e_trap_ctc, :e_trap]#, :e_cusp]

# load configuration for calibration
filekey = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])[1]
ecal_config = PropDict(dataprod_config(asic).energy(filekey).default)

# do calibration
process_energy_calibration(asic, period, run, category, channel, ecal_config; reprocess = reprocess, e_types = e_types)

# read calibration parameters

asic.par.rpars.ecal[period, run, channel].e_trap
asic.par.rpars.ecal[period, run, channel].e_trap_ctc.fit.Co60a.fwhm
asic.par.rpars.ecal[period, run, channel].e_trap_ctc.fit.Co60b.fwhm


