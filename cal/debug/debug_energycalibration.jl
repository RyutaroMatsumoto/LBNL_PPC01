# DEBUG script 
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
run = DataRun(2)
channel = ChannelId(1)
category = :cal 
e_types = [:e_trap]#, :e_cusp]

# load configuration for calibration
filekey = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])[1]
ecal_config = dataprod_config(asic).energy(filekey).default

# DEBUG START 
data = asic
source = :co60
calib_type = :gamma
gamma_lines =  ecal_config[Symbol("$(source)_lines")]
left_window_sizes = ecal_config[Symbol("$(source)_left_window_sizes")]
right_window_sizes = ecal_config[Symbol("$(source)_right_window_sizes")]
gamma_names = ecal_config[Symbol("$(source)_names")]
fit_funcs = Symbol.(ecal_config[Symbol("$(source)_fit_func")])
gamma_lines_dict = Dict(gamma_names .=> gamma_lines)

# read quality cuts from pars
qc = data.par.rpars.qc[period][run][channel]

# read dsp parameters
dsp_pars = Table(read_ldata(data, :jldsp, category, period, run, channel);)
fs = 12
e_type = :e_trap 
   
# select energy filter and apply qc
e_uncal = getproperty(dsp_pars, Symbol("$e_type"))[findall(qc.wvf_keep.all)]
e_uncal_func = "$e_type"

# do simple calibration and plot 
result_simple, report_simple = simple_calibration(e_uncal, gamma_lines , left_window_sizes, right_window_sizes,; 
            calib_type = calib_type, binning_peak_window=ecal_config.binning_peak_window, quantile_perc=NaN, 
            peak_quantile= ecal_config.left_peak_quantile..ecal_config.right_peak_quantile, 
            bin_quantile = ecal_config.left_bin_quantile..ecal_config.right_bin_quantile, 
            peakfinder_threshold = ecal_config.peakfinder_threshold, 
            peakfinder_σ = ecal_config.peakfinder_σ);

fig = LegendMakie.lplot(report_simple, juleana_logo = false)
[delete!(leg) for leg in fig.content if leg isa Legend]
vl = vlines!([report_simple.peak_guess * ustrip(report_simple.c)], color = :red2, label = "Peak Guess", alpha = 0.5, linewidth = 3)
axislegend(Makie.current_axis(), [vl], ["Peak Guess"], position = :lt)
Makie.xlims!(Makie.current_axis(), 300, 2000)
fig
m_cal_simple = result_simple.c
filekey = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])[1]

report_simple2 = (h_calsimple = report_simple.h_calsimple, 
h_uncal = report_simple.h_uncal,
c = report_simple.c,
fep_guess = report_simple.peak_guess,
peakhists = report_simple.peakhists,
peakstats = report_simple.peakstats)
   
plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(data, filekey, :energy_calibration) * "/"
det = _channel2detector(data, channel)
p = LegendMakie.lplot(report_simple2, title = get_plottitle(filekey, det, "Simple Calibration"; additiional_type=string(e_type)), cal = true, juleana_logo = false)
Makie.current_axis().titlesize = 16; p

# # fit peaks
@debug "Fit all peaks"
result_fit, report_fit = fit_peaks(result_simple.peakhists, result_simple.peakstats, gamma_names; 
                            e_unit=result_simple.unit, calib_type=:th228, fit_func = fit_funcs, m_cal_simple=m_cal_simple)
fig_fit = LegendMakie.lplot(report_fit, figsize = (600, 400*length(report_fit)), title = get_plottitle(filekey, det, "Peak Fits"; additiional_type=string(e_type)), juleana_logo = false)
fig_fit
pname = plt_folder * _get_pltfilename(data, filekey, channel, Symbol("peak_fits_$(e_type)"))
save(pname, fig_fit)
@info "Save peak fits plot to $(pname)"

# # calibration curve 
gamma_names_cal_fit = [p for p in gamma_names if !(p in Symbol.(ecal_config.cal_fit_excluded_peaks))]
μ_fit =  [result_fit[p].centroid for p in gamma_names_cal_fit]
pp_fit = [gamma_lines_dict[p] for p in gamma_names_cal_fit]
result_calib, report_calib = fit_calibration(ecal_config.cal_pol_order, μ_fit, pp_fit; e_expression=e_uncal_func)
@debug "Found $e_type calibration curve: $(result_calib.func)"
# plot calibration curve 
μ_notfit =  [result_fit[p].centroid for p in gamma_names if !(p in gamma_names_cal_fit)]
pp_notfit = [gamma_lines_dict[p] for p in gamma_names if !(p in gamma_names_cal_fit)]
pname_calib = plt_folder * _get_pltfilename(data, filekey, channel, Symbol("calibration_curve_$(e_type)"))
# Plots_theme(; fs = fs, grid = :on)
fig_calib = if isempty(μ_notfit)
    # p = plot(report_calib, xerrscaling=100)
    LegendMakie.lplot(report_calib, xerrscaling=100, title = get_plottitle(filekey, det, "Calibration Curve"; additiional_type=string(e_type)), juleana_logo = false)
else
    # p = plot(report_calib, xerrscaling=100, additional_pts=(μ = μ_notfit, peaks = pp_notfit))
    LegendMakie.lplot(report_calib, xerrscaling=100, additional_pts=(μ = μ_notfit, peaks = pp_notfit), title = get_plottitle(filekey, det, "Calibration Curve"; additiional_type=string(e_type)), juleana_logo = false)
end
Makie.current_axis().titlesize = 17
# plot!(fig_calib, thickness_scaling = 1.5, legend = :topleft,  ylabel = "Energy (a.u.)")
# plot!(plot_title=get_plottitle(filekey, _channel2detector(data, channel), " Calibration Curve"; additiional_type=string(e_type)), plot_titlelocation=(0.5,0.3), plot_titlefontsize=7)
# savefig(fig_calib, pname_calib)
save(pname_calib, fig_calib)
@info "Save calibration curve plot to $(pname_calib)"

# # resolution curve 
f_cal_widths(x) = report_calib.f_fit(x) .* report_calib.e_unit .- first(report_calib.par)
gamma_names_fwhm_fit = [p for p in gamma_names if !(p in Symbol.(ecal_config.fwhm_fit_excluded_peaks))]
fwhm_fit = f_cal_widths.([result_fit[p].fwhm for p in gamma_names_fwhm_fit])
pp_fit = [gamma_lines_dict[p] for p in gamma_names_fwhm_fit] 
result_fwhm, report_fwhm = fit_fwhm(ecal_config.fwhm_pol_order, pp_fit, fwhm_fit; e_type_cal=Symbol("$(e_type)_cal"), e_expression=e_uncal_func, uncertainty=true)
@debug "Found $e_type FWHM: $(round(u"keV", result_fwhm.qbb, digits=2))"
# plot resolution curve
fwhm_notfit =  f_cal_widths.([result_fit[p].fwhm for p in gamma_names if !(p in gamma_names_fwhm_fit)])
pp_notfit = [gamma_lines_dict[p] for p in gamma_names if !(p in gamma_names_fwhm_fit)]
pname_fwhm = plt_folder * _get_pltfilename(data, filekey, channel, Symbol("fwhm_$(e_type)"))
fig_fwhm = if isempty(fwhm_notfit)
    # plot(report_fwhm)
    LegendMakie.lplot(report_fwhm, title = get_plottitle(filekey, det, "FWHM"; additiional_type=string(e_type)), juleana_logo = false)
else
    # plot(report_fwhm, additional_pts=(peaks = pp_notfit, fwhm = fwhm_notfit))
    LegendMakie.lplot(report_fwhm, additional_pts=(peaks = pp_notfit, fwhm = fwhm_notfit), title = get_plottitle(filekey, det, "FWHM"; additiional_type=string(e_type)), juleana_logo = false)
end 
# plot!(fig_fwhm, thickness_scaling = 1.5, xticks = 0:500:3000, legend = :bottomright, xlabel = "Energy (keV)", right_margin = 3mm)
# plot!(plot_title=get_plottitle(filekey, _channel2detector(data, channel), "FWHM"; additiional_type=string(e_type)), plot_titlelocation=(0.5,0.3), plot_titlefontsize=8)
# savefig(fig_fwhm, pname_fwhm)
save(pname_fwhm, fig_fwhm)
@info "Save FWHM plot to $(pname_fwhm)"

# calibrate best fit values using energy calibration curve 
f_cal_pos(x) = report_calib.f_fit(x) .* report_calib.e_unit
width_pars = [:σ, :fwhm]
position_pars = [:μ, :centroid]
for (peak, result_peak) in result_fit
    result_peak = merge(result_peak, NamedTuple{Tuple(width_pars)}([f_cal_widths(result_peak[k]) for k in width_pars]...))
    result_peak = merge(result_peak, NamedTuple{Tuple(position_pars)}([f_cal_pos(result_peak[k]) for k in position_pars]...))
    result_fit[peak] = result_peak
end
result_fit[gamma_names[1]].fwhm
result_fit[gamma_names[2]].fwhm

result_energy = (
    m_cal_simple = m_cal_simple,
    fwhm = result_fwhm,
    cal = result_calib,
    fit  = result_fit,
)

