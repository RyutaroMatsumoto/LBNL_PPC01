# This is a script that performs a very simple processing of HPGe waveforms. You will need to adjust some settings like baseline, tail windows, and t0 windows,...
# It is intended to be used for quick tests and to get a first impression of the data.
# 1. load waveforms (lh5)
# 2. get decay times for pole-zero correction, and apply pole-zero correction
# 3. do a simple DSP using a trapezoidal filter
# 4. rough quality cuts 
# 5. energy calibration by selecting a peak "by eye" and calibrating the spectrum
# 6. fit gamma peaks with truncated Gauss (not full gama peak shape. results will be approximate!)
# plot gamma peak fits and print FWHM and peak position
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendDSP
using LegendDSP: get_fltpars
using LegendHDF5IO
using LegendSpecFits
using RadiationSpectra
using RadiationDetectorDSP
using Measurements: value as mvalue, uncertainty as muncert
using PropDicts
using StatsBase, IntervalSets
using Unitful
using TypedTables
using Plots
using Measures
using Distributions
using Printf

# load functions from hpge-ana
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/processing_funcs/process_decaytime.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")
include("$(@__DIR__)/$relPath/utils/utils_plot.jl")
Plots_theme()

# setup 
asic = LegendData(:ppc01)
period = DataPeriod(3)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)
det = _channel2detector(asic, channel)

# 1. load waveforms 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
data = read_ldata(asic, DataTier(:raw), filekeys, channel)
wvfs = data.waveform
Plots.plot(wvfs[rand(1:length(wvfs))], title = "Example waveform", xlabel = "Time", ylabel = "Amplitude [V]", label = "raw waveform")

bl_window = 0.0u"µs" .. 7.0u"µs"
tail_window = 9.0u"µs" .. 40.0u"µs"
Plots.vline!([leftendpoint(bl_window), rightendpoint(bl_window)], label = "baseline", color = :red)
Plots.vline!([leftendpoint(tail_window), rightendpoint(tail_window)], label = "tail", color = :violet)

# 2. get decay times for pole-zero correction: 
decay_times = dsp_decay_times(wvfs, bl_window, tail_window)
filter!(x -> 0u"µs" < x < 1000u"µs", decay_times)
Plots.stephist(decay_times, bins = 100, xlabel = "Decay time", ylabel = "Counts", title = "Decay times of waveforms", label = false)
decay_time_pz = []
try 
    min_τ = 250.0u"µs"
    max_τ = 350.0u"µs"
    rel_cut_fit = 0.25 
    nbins = 100

    cuts_τ = cut_single_peak(decay_times, min_τ, max_τ,; n_bins=nbins, relative_cut=rel_cut_fit)
    result_τ, report = fit_single_trunc_gauss(decay_times, cuts_τ)
    p = Plots.plot(report, legend = :outertop)
    display(p)
    @info "truncated gaus fit - done"
    decay_time_pz = mvalue(result_τ.µ) # decay time for pole-zero 
catch e
    decay_time_pz = median(decay_times)
end 

## 3. do very simple dsp using trap filter 
t0_threshold = 200#0.01
# shift waveform to baseline
bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))
wvfs = shift_waveform.(wvfs, -bl_stats.mean)# substract baseline from waveforms

# pole-zero correction
deconv_flt = InvCRFilter(decay_time_pz)
wvfs_pz = deconv_flt.(wvfs)
i = rand(1:length(wvfs))
Plots.plot(wvfs[i], title = "Example waveform", xlabel = "Time (µs)", ylabel = "Amplitude (V)", label = "raw")
Plots.plot!(wvfs_pz[i], label = "pz - corrected")

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
trap_rt = 1.5u"µs"
trap_ft = 1.0u"µs"
uflt_trap_rtft = TrapezoidalChargeFilter(trap_rt, trap_ft)
signal_estimator = SignalEstimator(PolynomialDNI(3, 100u"ns"))
e_trap = signal_estimator.(uflt_trap_rtft.(wvfs), t50 .+ (trap_rt + trap_ft/2))

dsp_par =   Table(blmean = bl_stats.mean, blslope = bl_stats.slope, 
                 tailmean = pz_stats.mean, tailslope = pz_stats.slope,
                t0 = t0, t10 = t10, t50 = t50, t90 = t90,
                e_max = wvf_max, e_min = wvf_min,
                e_trap = e_trap)
# DSP done... 

# 4. apply some rough quality cuts  
t0_window = 5.0u"µs"..9.0u"µs"
# t0_window = 17.0u"µs"..20.0u"µs"
# t0_window = 20.0u"µs"..55.0u"µs" 
# t0_window = 35.0u"µs"..30.0u"µs" 
dsp_par_qc = filter(x -> leftendpoint(t0_window) .< x.t0 .< rightendpoint(t0_window), dsp_par)
filter!(x ->  0.0 .< x.e_trap, dsp_par_qc)
# filter!(x -> abs(x.blslope) .< 0.002 * 1/u"µs", dsp_par)

# look at uncalibrated spectrum: 
Plots.stephist(dsp_par.e_trap, bins = 1000, xlabel = "Energy (a.u.)", ylabel = "Counts", title = "Uncalibrated energy spectrum (e_trap)", label = false, fill = true)

# energy calibration "by eye" 
gamma_lines =  [1173.237u"keV", 1332.501u"keV"]
window_size =  [10.0u"keV", 10.0u"keV"]
e_uncal_range = 900.0..1150.0
Plots.stephist(filter(x-> leftendpoint(e_uncal_range) <= x <= rightendpoint(e_uncal_range), dsp_par.e_trap), bins = 1000, xlabel = "Energy (a.u.)", ylabel = "Counts", title = "Uncalibrated energy spectrum (e_trap)", label = false, fill = true)
peapos = 957#0.2355 #
vline!([peapos], linewidth = 1.5, color = :red2, label =  "peak guess")
cal_simple = gamma_lines[1] / peapos

# roughly calibrated spectrum 
e_cal = dsp_par.e_trap .* cal_simple
Plots.stephist(e_cal, xlabel = "Energy (keV)", bins = 0:5:3000, title = "Roughly calibrated energy spectrum (e_trap)", legend = false, fill = true)

# fit gamma peaks with truncated Gauss (not full gama peak shape. results will be approximate!)
function _fit_hist(energy, gamma_line, gamma_window)
    e_peak = filter(x-> gamma_line - gamma_window <= x < gamma_line + gamma_window, energy)

    cuts = cut_single_peak(e_peak, gamma_line - gamma_window, gamma_line + gamma_window,; n_bins = 100, relative_cut=0.2)
    result, report = fit_single_trunc_gauss(e_peak, cuts)
end

results = _fit_hist.(Ref(e_cal), gamma_lines, window_size)
result_fit = map(x-> x.result, results)
report_fit = map(x-> x.report, results)

µ = map(x-> x.µ, result_fit)
fwhm = map(x-> x.σ*2.355, result_fit)

# plot fit 
plts = [Plots.plot(report_fit[i], legend = false, title = @sprintf("\n FWHM = (%.1f ± %.1f) keV at E = %.0f keV , ", mvalue(ustrip(fwhm[i])), muncert(ustrip(fwhm[i])),mvalue(ustrip(gamma_lines[i]))), titlefontsize = 16) for i in eachindex(report_fit)] #
fig = Plots.plot(plts..., layout = (2, 1), size = (550, 710), left_margin = 10mm, thickness_scaling = 0.7, xlabel = "Energy (keV)")


# save figure 
plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(asic, filekeys[1], :spectrum_approx) * "/"
if !isdir(plt_folder)
    mkdir(plt_folder)
end
plt_name = plt_folder * _get_pltfilename(asic, filekeys[1], channel, Symbol("spectrum_approx_e_trap"))
savefig(fig, plt_name)
