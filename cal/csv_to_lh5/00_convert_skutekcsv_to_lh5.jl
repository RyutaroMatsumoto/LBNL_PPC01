# convert .ecsv files from SkuTek "FemtoDAQ Vireo" to LEGEND-style hdf5-files. 
using LegendDataManagement
using TypedTables
using Plots 
using CSV, JSON
using Unitful
using Dates

# include relevant functions
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/utils/utils_IO.jl")

# inputs
asic = LegendData(:ppc01)
period = DataPeriod(3)
run = DataRun(43)
channel = ChannelId(1)
category = :cal 

# input for function: csv folder, heading 
csv_folder = asic.tier[DataTier(:raw_csv), category , period, run]
timestep = 0.01u"µs"

skutek_csv_to_lh5(asic, period, run, category, channel, csv_folder; timestep = timestep)

# test resting this data 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
data_raw = Table(read_ldata(asic, DataTier(:raw), filekeys[1], channel))
