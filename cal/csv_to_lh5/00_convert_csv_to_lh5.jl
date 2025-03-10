# convert csv files from oscilloscope to LEGEND-style hdf5-files. 
using LegendDataManagement
using TypedTables
using Plots 

# include relevant functions
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/utils/utils_IO.jl")

<<<<<<< HEAD:cal/csv_to_lh5/00_convert_csv_to_lh5.jl
# inputs
asic = LegendData(:ppc01)
period = DataPeriod(1)
run = DataRun(14)
=======
# set data configuration (where to find data; and where to save results)
# if gethostname() == "Lisas-MacBook-Pro.local"
#    ENV["LEGEND_DATA_CONFIG"] = "/Users/lisa/Documents/Workspace/LEGEND/LBL_ASIC/ASIC_data/ppc01/config.json"
# else # on NERSC 
    # run(`hostname -d`)
    ENV["LEGEND_DATA_CONFIG"] = "/global/cfs/projectdirs/m2676/data/teststands/lbnl/PPC01_ryutaro/config.json"
# end 

# inputs
asic = LegendData(:ppc01)
period = DataPeriod(1)
run = DataRun(2)
>>>>>>> bd32341c0f9bde80f2aff9a94eedd2234c34bb9c:cal/p01/00_convert_csv_to_lh5.jl
channel = ChannelId(1)
category = :cal 
ti = 0.0u"µs"..550.0u"µs" # time interval within waveform is truncated and saved (to reduce file size)
csv_folder = asic.tier[DataTier(:raw_csv), category , period, run]

# convert csv files to lh5 files 
csv_to_lh5(asic, period, run, category, channel, csv_folder; csv_heading = 14, nChannels = 1, nwvfmax = NaN, ti = ti)

### Cross checks 
# read waveforms as sanity check 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
data = read_ldata(asic, DataTier(:raw), filekeys, channel)
Table(data)

# how many waveforms had to be removed? 
Plots.plot(1:length(data.eventnumber)-1, diff(data.eventnumber).-1, xlabel = "eventnumber", ylabel = "number of waveforms removed between events")

# how long did measurement take
time_daq_h = (data.timestamp[end] - data.timestamp[1])/(60*60)
time_per_wvf_h = (time_daq_h)/4700
time_per_100kwvf_h = time_per_wvf_h * 100e3

# confirm scope samping: number of data points 
data.waveform[1].signal
uconvert(u"ns",step(data.waveform[1].time))

# plot random waveform to check if conversion was successful
begin
    i = rand(1:length(data.waveform))
    Plots.plot(data.waveform[i] , label = "waveform $i")
end 

