using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using StatsBase
using LegendSpecFits
using Makie, LegendMakie, CairoMakie   
using PropDicts

# include relevant functions 
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
# include("$(@__DIR__)/$relPath/src/simple_dsp.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")
include("$(@__DIR__)/$relPath/processing_funcs/process_peakfits.jl")

# inputs 
asic = LegendData(:ppc01)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:bch)

for r in 13:22
    run = DataRun(r)
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])

# read dsp pars and plot 
dsp_pars = read_ldata(asic, :jldsp, category, period, run, channel);
Table(dsp_pars)

process_peakfits(asic, period, run, category, channel; reprocess = true, juleana_logo = false)
end 
# ### 
# # fit csa output: 
# data = asic
# det = _channel2detector(data, channel)
# filekey = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])[1]

# # waveforms 
# e_uncal = filter!(isfinite, dsp_pars.e_trap)
# qmin, qmax = 0.0, 1.0
# emin = quantile(e_uncal, qmin)
# emax = quantile(e_uncal, qmax)
# nbins = 50
# rel_cut_fit = 0.05
# cuts_τ = cut_single_peak(e_uncal, emin, emax,; n_bins=nbins, relative_cut=rel_cut_fit)
# result, report = fit_single_trunc_gauss(e_uncal, cuts_τ)
# p = LegendMakie.lplot(report, figsize = (600, 430), titlesize = 17, title = get_plottitle(filekey, det, "Energy Distribution"), juleana_logo = false, xlabel = "Energy (ADC)")
# # savelfig(LegendMakie.lsavefig, p, data, filekey, channel, :ecal)
# @info "Found mean energy time $(round(result.µ, digits=2)) for channel $channel / det $det"
# # result_wvf = (µ = result.μ, fit = result)


# e_pulser = filter!(isfinite, dsp_pars.pulser_e_trap)
# qmin, qmax = 0.0, 1.0
# emin = quantile(e_pulser, qmin)
# emax = quantile(e_pulser, qmax)
# nbins = 50
# rel_cut_fit = 0.05
# cuts_τ = cut_single_peak(e_pulser, emin, emax,; n_bins=nbins, relative_cut=rel_cut_fit)
# result_pulser, report_pulser = fit_single_trunc_gauss(e_pulser, cuts_τ)
# p = LegendMakie.lplot(report_pulser, figsize = (600, 430), titlesize = 17, title = get_plottitle(filekey, det, "Pulser Energy Distribution"), juleana_logo = false, xlabel = "Energy (ADC)")
# # savelfig(LegendMakie.lsavefig, p, data, filekey, channel, :ecal)
# @info "Found mean energy time $(round(result_pulser.µ, digits=2)) for channel $channel / det $det"
# # result_pulser = (µ = result_pulser.μ, fit = result)

# result_ecal = (µ = result.μ, µ_pulser = result_pulser.μ, fit = result, fit_pulser = result_pulser)
# writelprops(data.par[category].rpars.ecal[period], run, PropDict("$channel" => result_ecal))
# @info "Saved pars to disk"


# # data.par.rpars.ecal[period, run, channel]

# # display(p)