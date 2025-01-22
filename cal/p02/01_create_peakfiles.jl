
# activate environment and load modules. 
env_path = "$(@__DIR__)/" * relpath(split(@__DIR__, "LBNL_PPC01")[1], @__DIR__)  * "/LBNL_PPC01/"
import Pkg; Pkg.activate(env_path)
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendDSP
using LegendDSP: get_fltpars
using LegendSpecFits: get_friedman_diaconis_bin_width
using LegendHDF5IO
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

# include relevant functions 
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/processing_funcs/process_peak_split.jl")
include("$(@__DIR__)/$relPath/src/simple_dsp.jl")
include("$(@__DIR__)/$relPath/src/apply_qc.jl")
include("$(@__DIR__)/$relPath/utils/utils_plot.jl")
include("$(@__DIR__)/$relPath/utils/utils_aux.jl")
# inputs
reprocess = true
asic = LegendData(:ppc01)
period = DataPeriod(2)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)
source = :co60

# load configs 
qc_config = asic.metadata.config.qc.qc_config.default
ecal_config = asic.metadata.config.energy.energy_config.default
dsp_config = DSPConfig(asic.metadata.config.dsp.dsp_config.default)

# modify peakfinder threshold and sigma if needed
# TODO write overwrites !!
ecal_config.peakfinder_σ = 3.0
ecal_config.peakfinder_threshold = 20

# AUTOMATIC PEAK SPLIITING DID NOT WORK CORRECTLY. 2nd peak is strongly suppressed. 
# thats why we do it the manual way, with tuned parameters using only Co60a peak 
mode = :manual 
if mode == :auto
    h_uncals, peakpos = process_peak_split(asic, period, run, category, channel ; source = source, reprocess = true,
                        ecal_config = ecal_config,
                        dsp_config = dsp_config )   
    # make sanity plot. did peakfinder work as expected?
    Makie_theme(; fs = 30, xgridvisible = false, ygridvisible = false)
    nrow = ceil(Int,length(h_uncals)/2)
    f1 = Figure(size = (600*nrow, 2*400))
    row = 1
    for i in eachindex(h_uncals)
        if mod(i, 2) == 1 # odd
            local ax = Axis(f1[row, 1] , xlabel = "Energy")
        elseif mod(i,2) == 0 #even
            local ax = Axis(f1[row, 2] , xlabel = "Energy")
            row = row + 1
        end
        plot!(ax, h_uncals[i])
        vlines!(peakpos[i], color = :red2, linewidth = 3, label = "peakpos")
        hlines!([ecal_config.peakfinder_threshold], color = :green, label = " peakfinder threshold")
        axislegend(ax, position = :rt)
    end
    f1
end 


# generate peakfiles manually 
# load filekeys 
data = asic
filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
peak_folder =  asic.tier[DataTier(:jlpeaks), category , period, run] * "/"
if !ispath(peak_folder)
mkpath(peak_folder)
@info "create path: $peak_folder"
end
peak_files = peak_folder .* string.(filekeys) .* "-tier_jlpeaks.lh5"
if reprocess == false
    if sum(isfile.(peak_files)) == length(filekeys)
        @info "all peak files are already processed. You're already done!"
        println("to read peak files: use read_ldata(data, :jlpeaks, category , period, run, channel)")
        return 
    else
        @info "$(sum(isfile.(peak_files))) peak files are already processed - skip"
    end
    filekeys = filekeys[.!isfile.(peak_files)]
end

for f in eachindex(filekeys) 
    filekey = filekeys[f]
    peak_file = peak_files[f]
    data_ch = read_ldata(data, DataTier(:raw), filekey, channel)
    wvfs = data_ch.waveform
    e_uncal = data_ch.daqenergy
    eventnumber = data_ch.eventnumber
    timestamp = data_ch.timestamp

    gamma_lines =  ecal_config.co60_lines[1]
    gamma_names =  ecal_config.co60_names[1]
    left_window_sizes = ecal_config.co60_left_window_sizes[1]
    right_window_sizes = ecal_config.co60_right_window_sizes[1]

    bin_width = get_friedman_diaconis_bin_width(filter(in(quantile(e_uncal, 0.1)..quantile(e_uncal, 0.9)), e_uncal))
    nbins = ceil(Int, (maximum(e_uncal)-minimum(e_uncal))/bin_width)
    h_uncal = fit(Histogram, e_uncal, nbins=nbins)
    _, peakpos = RadiationSpectra.peakfinder(h_uncal, σ=ecal_config.peakfinder_σ, backgroundRemove=true, threshold = ecal_config.peakfinder_threshold)
   
    # sanity plot: 
    Plots_theme(; fs = 20)
    begin 
        p = plot(h_uncal, color = :dodgerblue, ylims = (0, :auto), label = false,  xlabel = "Energy (a.u.)", linewidth = 0)
        vline!(peakpos, color = :red2, label = "peakfinder result", linewidth = 3)
        hline!([ecal_config.peakfinder_threshold], color = :orange, linewidth = 2.5, linestyle = :dash, label = "peakfinder threshold")
    end
    title!(get_plottitle(filekey, _channel2detector(data, channel), "peak split"); titlefontsize = 15)
    plt_simple = savelfig(savefig, p, data, filekey, channel, Symbol("peak_split"))
    dest = join(split(plt_simple,"/")[1:end-1] .* "/") * replace(split(plt_simple,"peak_split/")[end], ".png" => "_$(split("$filekey", "-")[end]).png")
    mv(plt_simple, dest, force = true)
    cal_simple = gamma_lines./peakpos[1]
    e_simplecal = e_uncal .* cal_simple

    # save
    if !ispath(peak_folder)
    mkpath(peak_folder)
    @info "create path: $peak_folder"
    end

    fid = lh5open(peak_file, "w") 
    for i in eachindex(gamma_lines)
        if length(gamma_lines) == 1
            gamma_name = gamma_names
        else
            gamma_name = gamma_names[i]
        end
        peakIdx = findall((gamma_lines[i] - left_window_sizes[i]) .<= e_simplecal .< (gamma_lines[i] + right_window_sizes[i]))
        # do simple dsp. only for peaks. 
        dsp_par = simple_dsp_qc(Table(waveform = wvfs[peakIdx]), dsp_config)
        qc_cuts, _ = apply_qc(dsp_par, qc_config)
        qc_flag = qc_cuts.wvf_keep.all
        qc_idx = findall(qc_flag)
        qc_surv = qc_cuts.qc_surv.all
        fid["$channel/jlpeaks/$(gamma_name)/waveform"] = wvfs[peakIdx][qc_idx]
        fid["$channel/jlpeaks/$(gamma_name)/daqenergy"] = e_uncal[peakIdx][qc_idx]
        fid["$channel/jlpeaks/$(gamma_name)/e_simplecal"] = e_simplecal[peakIdx][qc_idx]
        fid["$channel/jlpeaks/$(gamma_name)/eventnumber"] = eventnumber[peakIdx][qc_idx]
        fid["$channel/jlpeaks/$(gamma_name)/timestamp"] = timestamp[peakIdx][qc_idx]
        fid["$channel/jlpeaks/$(gamma_name)/gamma_line"] = fill(gamma_lines[i], length(wvfs[peakIdx]))
        fid["$channel/jlpeaks/$(gamma_name)/qc_flag"] = qc_flag[qc_idx]# fill(NaN, length(wvfs[peakIdx]))qc_flag[qc_idx]
        fid["$channel/jlpeaks/$(gamma_name)/qc_surv"] = fill(qc_surv, length(qc_idx)) #qc_surv#fill(NaN, length(wvfs[peakIdx]))
        for par in columnnames(dsp_par)
            fid["$channel/jlpeaks/$(gamma_name)/$par"]  = getproperty(dsp_par, par)
        end
    end
    println("peak file processing done for $filekey")
    close(fid)
end 

# # sanity check: open peak files and plot 
peakA = read_ldata(:Co60a, asic, :jlpeaks, category, period, run, channel)
# peakB = read_ldata(:Co60b, asic, :jlpeaks, category, period, run, channel)
