# VERY simple proessing without any filter optimization or classic processing ChannelId
# for quick tests 
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
include("$(@__DIR__)/$relPath/utils/utils_plot.jl")
Plots_theme()

# setup 
asic = LegendData(:ppc01)
period = DataPeriod(2)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)


# 1. load waveforms 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
data = read_ldata(asic, DataTier(:raw), filekeys, channel)
wvfs = data.waveform
plot(wvfs[rand(1:length(wvfs))], title = "Example waveform", xlabel = "Time", ylabel = "Amplitude [V]", legend = false)

# 2. get decay times for pole-zero correction: 
bl_window = 0.0u"µs" .. 20.0u"µs"
tail_window = 30.0u"µs" .. 130.0u"µs"
decay_times = dsp_decay_times(wvfs, bl_window, tail_window)
filter!(x -> 0u"µs" < x < 1000u"µs", decay_times)
stephist(decay_times, bins = 100, xlabel = "Decay time", ylabel = "Counts", title = "Decay times of waveforms", label = false)
decay_time_pz = 200.0u"µs" # choose decay time for pole-zero correction by hand 

## 3. do very simple dsp using trap filter 
t0_threshold = 0.01

# shift waveform to baseline
bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))
wvfs = shift_waveform.(wvfs, -bl_stats.mean)# substract baseline from waveforms

# pole-zero correction
deconv_flt = InvCRFilter(decay_time_pz)
wvfs_pz = deconv_flt.(wvfs)
i = rand(1:length(wvfs))
plot(wvfs[i], title = "Example waveform", xlabel = "Time (µs)", ylabel = "Amplitude (V)", label = "raw")
plot!(wvfs_pz[i], label = "pz - corrected")

# tail_stats = tailstats.(wvfs, leftendpoint(tail_window), rightendpoint(tail_window)) # tail analysis 
# # get raw wvf maximum/minimum
wvf_max = maximum.(wvfs.signal)
wvf_min = minimum.(wvfs.signal)

# get tail mean, std and slope
pz_stats = signalstats.(wvfs, leftendpoint(tail_window), rightendpoint(tail_window))

# characteristic times in waveform 
t0 = get_t0(wvfs, t0_threshold)
t10 = get_threshold(wvfs, wvf_max .* 0.1)
t50 = get_threshold(wvfs, wvf_max .* 0.5)
t90 = get_threshold(wvfs, wvf_max .* 0.9)

# trap-filter: signal estimator for precise energy reconstruction
trap_rt = 3.0u"µs"
trap_ft = 1.0u"µs"
uflt_trap_rtft = TrapezoidalChargeFilter(trap_rt, trap_ft)
signal_estimator = SignalEstimator(PolynomialDNI(3, 100u"ns"))
e_trap = signal_estimator.(uflt_trap_rtft.(wvfs), t50 .+ (trap_rt + trap_ft/2))

dsp_par =   Table(blmean = bl_stats.mean, blslope = bl_stats.slope, 
                 tailmean = pz_stats.mean, tailslope = pz_stats.slope,
                t0 = t0, t10 = t10, t50 = t50, t90 = t90,
                e_max = wvf_max, e_min = wvf_min,
                e_trap = e_trap)
# DSP done. 

# 4. apply some very rough quality cuts   
filter!(x -> 20.0u"µs" .< x.t0 .< 27.0u"µs", dsp_par)
filter!(x ->  x.e_trap .> 0.0, dsp_par)
filter!(x -> abs(x.blslope) .< 0.002 * 1/u"µs", dsp_par)

# energy calibration: very rough. 
ecal_config = asic.metadata.config.energy.energy_config.default
source = :co60
calib_type = :gamma
gamma_lines =  [ecal_config[Symbol("$(source)_lines")][1]]
left_window_sizes = [ecal_config[Symbol("$(source)_left_window_sizes")][1]]
right_window_sizes = [ecal_config[Symbol("$(source)_right_window_sizes")][1]]
gamma_names = [ecal_config[Symbol("$(source)_names")][1]]
fit_funcs = [Symbol.(ecal_config[Symbol("$(source)_fit_func")])[1]]
gamma_lines_dict = Dict(gamma_names .=> gamma_lines)
e_uncal = dsp_par.e_trap
p = stephist(e_uncal)
q = 0.5
vline!([quantile(e_uncal,q)])
cal_simple = gamma_lines[1] / quantile(e_uncal,q)


# roughly calibrated spectrum 
e_cal = e_uncal .* cal_simple
stephist(e_cal, xlabel = "Energy", nbins = 100)
e_cal_cut = filter(x-> x > 800u"keV", e_cal)
filter!(x-> x < 1600u"keV", e_cal_cut)
stephist(e_cal_cut, nbins = 100)

using Distributions
result_fit = fit(Normal, ustrip.(e_cal_cut))
fwhm_approx = result_fit.σ * 2.355u"keV"