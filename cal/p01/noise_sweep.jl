using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using LegendSpecFits
using LegendDSP
using LegendDSP: get_fltpars
using RadiationDetectorSignals
using RadiationDetectorDSP
using Measurements: value as mvalue
using PropDicts
using StatsBase, IntervalSets
using Unitful
using TypedTables
using CairoMakie, LegendPlots
using Measures
using Optim
using BSplineKit
using Printf
# set data configuration (where to find data; and where to save results)

# include relevant functions 
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/src/filteropt_rt_optimization_blnoise.jl")
include("$(@__DIR__)/$relPath/utils/utils_plot.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")
include("$(@__DIR__)/$relPath/processing_funcs/process_filteropt.jl")

# inputs 
asic = LegendData(:ppc01)
period = DataPeriod(1)
run = DataRun(10)
begin 
channel = ChannelId(1)
category = DataCategory(:cal)
filter_types = [:trap, :cusp]
rt_opt_mode = :bl_noise

# load configs and modify if needed 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
pz_config = dataprod_config(asic).dsp(filekeys[1]).pz.default
peak =  Symbol(pz_config.peak)
τ_pz = mvalue(get_values(asic.par.rpars.pz[period, run, channel]).τ)

# Load DSP config and modify rise-time grid for trap filter
# dsp_config = DSPConfig(dataprod_config(asic).dsp(filekeys[1]).default)
dsp_config_pd = dataprod_config(asic).dsp(filekeys[1]).default#.e_grid_trap.rt
dsp_config_pd.e_grid_trap.rt = PropDict( :start => 1.0u"µs", :stop => 5.0u"µs", :step => 0.25u"µs")
dsp_config = DSPConfig(dsp_config_pd)

data = read_ldata(asic, DataTier(:raw), filekeys, channel)
wvfs = data.waveform
def_ft = 2.0u"µs"

filter_type = :trap 
grid_rt   = getproperty(dsp_config, Symbol("e_grid_rt_$(filter_type)"))
result_rt, report_rt = filteropt_rt_optimization_blnoise(filter_type, wvfs, dsp_config, τ_pz; ft = def_ft)

fig = Figure()
ax = Axis(fig[1, 1],  xlabel = "Rise time ($(unit(report_rt.rt)))", ylabel = "RMS noise (V)", 
        title = "Noise sweep ($filter_type) $period-$run-$channel",
        xgridvisible = true, ygridvisible = true,)
    p = lines!(ax, ustrip.(collect(report_rt.enc_grid_rt)), report_rt.enc, 
    label = "trap. filter: ft = $(def_ft)")
    scatter!(ax, ustrip.(collect(report_rt.rt)), report_rt.min_enc, color = :orange, 
    label = @sprintf("min ENC = %.1e V",report_rt.min_enc))
axislegend(position = :lt)
fig
plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(asic, filekeys[1], :noise) * "/"
if !isdir(plt_folder)
    mkdir(plt_folder)
end
plt_name = plt_folder * _get_pltfilename(asic, filekeys[1], channel, Symbol("noise-sweep_$(filter_type)_ft$(def_ft)"))
save(plt_name, fig)
end 
fig