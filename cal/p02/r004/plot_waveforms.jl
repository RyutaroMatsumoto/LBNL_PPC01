## analysis of PPC01 data - period 02 (GFET)

# activate environment and load modules. 
env_path = "$(@__DIR__)/" * relpath(split(@__DIR__, "LBNL_PPC01")[1], @__DIR__)  * "/LBNL_PPC01/"
import Pkg; Pkg.activate(env_path)
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using HDF5, LegendHDF5IO
using PropDicts
using Unitful
using TypedTables
using CairoMakie 
using Measures

# set data configuration (where to find data; and where to save results)
ENV["LEGEND_DATA_CONFIG"] = "/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json"#"/Users/lisa/Documents/Workspace/LEGEND/LBL_ASIC/ASIC_data/ppc01/config.json"
# ENV["LEGEND_DATA_CONFIG"] = "/Users/lisa/Documents/Workspace/LEGEND/LBL_ASIC/ASIC_data/ppc01/config.json"

# load functions from hpge-ana
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__)
include("$relPath/hpge-ana/utils/utils_plot.jl")

# inputs and setup 
period = DataPeriod(2)
run = DataRun(1)
channel = ChannelId(1)
category = :cal 
asic = LegendData(:ppc01)

# read waveforms (raw tier)
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
data = read_ldata(asic, DataTier(:raw), filekeys, channel)
Table(data)

# how long did measurement take
time_daq_h = (data.timestamp[end] - data.timestamp[1])/(60*60)
time_per_wvf_h = (time_daq_h)/4700

time_per_100kwvf_h = time_per_wvf_h * 100e3

# confirm scope samping: number of data points 
data.waveform[1].signal
uconvert(u"ns",step(data.waveform[1].time))

# plot waveforms 
Makie_theme(; fs = 23, xgridvisible = true, ygridvisible = true)
plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(asic, filekey, :waveform)
begin
    i = rand(1:length(data.waveform))
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1, 1], title = label = "Waveform $i", xlabel = "Time ($(unit(data.waveform[i].time[1])))", ylabel = "Signal (V)")
    p = lines!(ax, data.waveform[i].time, data.waveform[i].signal , color = :dodgerblue, linewidth = 3, label = "waveform $i")
    # axislegend(ax, position = :rb)
    fig
    pname = plt_folder * "/waveform_$(i).png"
    save(pname, fig)
    @info "Saved waveform plot to $pname"
    fig
end 

# plot energy spectrum (DAQ energy = waveform max)
f2 = Figure(size = (600, 400))
ax2 = Axis(f2[1, 1] )
hist!(ax2, data.daqenergy, bins = 1000)
f2 
