env_path = "$(@__DIR__)/" * relpath(split(@__DIR__, "LBNL_PPC01")[1], @__DIR__)  * "/LBNL_PPC01/"
import Pkg; Pkg.activate(env_path)
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
using Plots 
using Unitful, Measures
using Measurements: value as mvalue
# set data configuration (where to find data; and where to save results)
ENV["LEGEND_DATA_CONFIG"] = "/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json"

# include relevant functions 
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")
include("$(@__DIR__)/$relPath/utils/utils_plot.jl")
include("$(@__DIR__)/$relPath/processing_funcs/process_energy_calibration.jl")

# inputs
reprocess = true 
asic = LegendData(:ppc01)
period = DataPeriod(2)
run = DataRun(1)
channel = ChannelId(1)
category = :cal 
e_type = :e_trap
ecal_config = asic.metadata.config.energy.energy_config.default
data = asic 
mode = :manual 

# do calibration 
if mode == :auto
    process_energy_calibration(asic, period, run, category, channel; reprocess = reprocess, ecal_config = ecal_config, e_types = [:e_trap, :e_cusp, :e_zac])
    # read calibration parameters
    asic.par.rpars.ecal[period, run, channel].e_trap
end 

# MANUAL 
result_dict = Dict{Symbol, NamedTuple}()
# load configuration for calibration
@info "Load calibration configuration for $(source) source"
if source == :th228
    calib_type = :th228
elseif source == :co60
    calib_type = :gamma
else
    error("Unknown source $(source)")
end
gamma_lines =  [ecal_config[Symbol("$(source)_lines")][1]]
left_window_sizes = [ecal_config[Symbol("$(source)_left_window_sizes")][1]]
right_window_sizes = [ecal_config[Symbol("$(source)_right_window_sizes")][1]]
gamma_names = [ecal_config[Symbol("$(source)_names")][1]]
fit_funcs = [Symbol.(ecal_config[Symbol("$(source)_fit_func")])[1]]
gamma_lines_dict = Dict(gamma_names .=> gamma_lines)

# read quality cuts from pars
qc = data.par.rpars.qc[period][run][channel]

# read dsp parameters
dsp_pars = Table(read_ldata(data, :jldsp, category, period, run, channel);)

fs = 12
function _energy_calibration(e_type::Symbol)
    if !reprocess && haskey(data.par.rpars.ecal[period, run, channel], e_type)
        @info "Load existing calibration pars for $(e_type)"
        return NamedTuple(data.par.rpars.ecal[period, run, channel][e_type])
    end
    # select energy filter and apply qc
    e_uncal = getproperty(dsp_pars, Symbol("$e_type"))[findall(qc.wvf_keep)]
    e_uncal_func = "$e_type"

    # do simple calibration and plot 
    result_simple, report_simple = simple_calibration(e_uncal, gamma_lines , left_window_sizes, right_window_sizes,; 
                calib_type = calib_type, binning_peak_window=ecal_config.binning_peak_window, quantile_perc=NaN, 
                peak_quantile= ecal_config.left_peak_quantile..ecal_config.right_peak_quantile, 
                bin_quantile = ecal_config.left_bin_quantile..ecal_config.right_bin_quantile, 
                peakfinder_threshold = ecal_config.peakfinder_threshold, 
                peakfinder_σ = ecal_config.peakfinder_σ);

    m_cal_simple = result_simple.c
    Plots_theme(; fs = fs, grid = :on)
    p = plot(report_simple, right_margin=4mm, top_margin = 0mm, yformatter=:plain, thickness_scaling=1.5, cal=true,  size = (620, 400))
    filekey = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])[1]
    title!(p, get_plottitle(filekey, _channel2detector(data, channel), "Simple Calibration"; additiional_type=string(e_type)), titlefontsize = 7)
    plt_simple = savelfig(savefig, p, data, filekey, channel, Symbol("simple_calibration_$(e_type)"))
    @info "Save simple calibration plot to $(plt_simple)"

    # # fit peaks
    @debug "Fit all peaks"
    result_fit, report_fit = fit_peaks(result_simple.peakhists, result_simple.peakstats, gamma_names; 
                                e_unit=result_simple.unit, calib_type=:th228, fit_func = fit_funcs, m_cal_simple=m_cal_simple)
    Plots_theme(; fs = fs+2, grid = :off)
    p = plot(broadcast(k -> plot(report_fit[k], left_margin=2mm, right_margin=2mm, top_margin=-5mm, bottom_margin=-2mm, title=string(k), ms=2), keys(report_fit))..., layout=(length(report_fit), 1), size=(1000,710*length(report_fit)) , thickness_scaling=1.8, titlefontsize = 10, legendfontsize = 8, yguidefontsize = 9, xguidefontsize=11)
    plot!(p, plot_title=get_plottitle(filekey, _channel2detector(data, channel), "Peak Fits"; additiional_type=string(e_type)), plot_titlelocation=(0.5,0.2), plot_titlefontsize = 12)
    plt_peaks = savelfig(savefig, p, data, filekey, channel , Symbol("peak_fits_$(e_type)"))
    @info "Save peak fits plot to $(plt_peaks)"

    # calibration curve 
    gamma_names_cal_fit = [p for p in gamma_names if !(p in Symbol.(ecal_config.cal_fit_excluded_peaks))]
    μ_fit =  [result_fit[p].centroid for p in gamma_names_cal_fit]
    pp_fit = [gamma_lines_dict[p] for p in gamma_names_cal_fit]
    result_calib, report_calib = fit_calibration(ecal_config.cal_pol_order, μ_fit, pp_fit; e_expression=e_uncal_func)
    @debug "Found $e_type calibration curve: $(result_calib.func)"
    # plot calibration curve 
    Plots_theme(; fs = fs, grid = :on)
    μ_notfit =  [result_fit[p].centroid for p in gamma_names if !(p in gamma_names_cal_fit)]
    pp_notfit = [gamma_lines_dict[p] for p in gamma_names if !(p in gamma_names_cal_fit)]
    if isempty(μ_notfit)
        p = plot(report_calib, xerrscaling=100)
    else
        p = plot(report_calib, xerrscaling=100, additional_pts=(μ = μ_notfit, peaks = pp_notfit))
    end
    plot!(p, thickness_scaling = 1.5, legend = :topleft,  ylabel = "Energy (a.u.)")
    plot!(plot_title=get_plottitle(filekey, _channel2detector(data, channel), " Calibration Curve"; additiional_type=string(e_type)), plot_titlelocation=(0.5,0.3), plot_titlefontsize=7)
    pcal = savelfig(savefig, p, data, filekey, channel, Symbol("calibration_curve_$(e_type)"))
    @info "Save calibration curve plot to $(pcal)"

    # resolution curve 
    f_cal_widths(x) = report_calib.f_fit(x) .* report_calib.e_unit .- first(report_calib.par)
    # gamma_names_fwhm_fit = [p for p in gamma_names if !(p in Symbol.(ecal_config.fwhm_fit_excluded_peaks))]
    # fwhm_fit = f_cal_widths.([result_fit[p].fwhm for p in gamma_names_fwhm_fit])
    # pp_fit = [gamma_lines_dict[p] for p in gamma_names_fwhm_fit] 
    # result_fwhm, report_fwhm = fit_fwhm(ecal_config.fwhm_pol_order, pp_fit, fwhm_fit; e_type_cal=Symbol("$(e_type)_cal"), e_expression=e_uncal_func, uncertainty=true)
    # @debug "Found $e_type FWHM: $(round(u"keV", result_fwhm.qbb, digits=2))"
    # # plot resolution curve
    # fwhm_notfit =  f_cal_widths.([result_fit[p].fwhm for p in gamma_names if !(p in gamma_names_fwhm_fit)])
    # pp_notfit = [gamma_lines_dict[p] for p in gamma_names if !(p in gamma_names_fwhm_fit)]
    # if isempty(fwhm_notfit)
    #     p = plot(report_fwhm)
    # else
    #     p = plot(report_fwhm, additional_pts=(peaks = pp_notfit, fwhm = fwhm_notfit))
    # end 
    # plot!(p, thickness_scaling = 1.5, xticks = 0:500:3000, legend = :bottomright, xlabel = "Energy (keV)", right_margin = 3mm)
    # plot!(plot_title=get_plottitle(filekey, _channel2detector(data, channel), "FWHM"; additiional_type=string(e_type)), plot_titlelocation=(0.5,0.3), plot_titlefontsize=8)
    # savelfig(savefig, p, data, filekey, channel, Symbol("fwhm_$(e_type)"))

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
    # result_fit[gamma_names[2]].fwhm

    result_energy = (
        m_cal_simple = m_cal_simple,
        fwhm = NaN,#result_fwhm,
        cal = result_calib,
        fit  = result_fit,
)
end

e_types = [:e_trap, :e_cusp, :e_zac]
for e_type in e_types
    result_dict[e_type] = _energy_calibration(e_type);
end 

result = PropDict(Dict("$channel" => result_dict))
writelprops(data.par.rpars.ecal[period], run, result)
@info "Saved pars to disk"


