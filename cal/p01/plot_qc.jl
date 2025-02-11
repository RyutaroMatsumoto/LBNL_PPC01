# to dsp par and their accepted qc cuts 
using LegendDataManagement
using LegendDataManagement.LDMUtils: get_pltfolder
using LegendHDF5IO
using Plots
using Measures
using PropDicts, TypedTables
using StatsBase
using Measurements: value as mvalue
using LegendDSP, RadiationDetectorDSP, IntervalSets
using Unitful, Measures, Printf

# include relevant functions 
relPath = relpath(split(@__DIR__, "hpge-ana")[1], @__DIR__) * "/hpge-ana/"
include("$(@__DIR__)/$relPath/utils/utils_plot.jl")

# inputs
asic = LegendData(:ppc01)
period = DataPeriod(1)
run = DataRun(6)
channel = ChannelId(1)
category = :cal 

# plot path
pltpath = get_pltfolder(asic, period, run, category, :qc) * "/"
if !ispath(pltpath)
    mkpath(pltpath)
end

# load dsp parameters and qc parameter 
dsp_par = Table(read_ldata(asic, :jldsp, category, period, run, channel))

# load quality cuts 
qc_cuts = asic.par.rpars.qc[period][run][channel]
qc_pars = filter!(x -> (x != :all) && (x != :finite), Symbol.(keys(qc_cuts.qc_surv)))

# plot histograms of dsp parameters with and without qc cuts
Plots_theme()
p = Vector{Plots.Plot}(undef, length(qc_pars))
for i in eachindex(p)
    par = ustrip.(getproperty(dsp_par, qc_pars[i])[qc_cuts.wvf_keep.finite])
    ustr = "$(unit(getproperty(dsp_par, qc_pars[i])[1]))"
    if ustr == "NoUnits" ustr = "a.u." end
    if qc_pars[i] == :blslope || qc_pars[i] == :tailslope
        filter!(x -> quantile(par, 0.025) < x < quantile(par, 0.975), par)
        exp_i = floor(log10(quantile(par, 0.975)) )
        par = par .* 10^-exp_i
        xl = string(qc_pars[i]) * " (10^$(round(Int, exp_i)) $ustr)"
    else
        xl = string(qc_pars[i]) * " ($ustr)"
        exp_i = 0.0
    end

    h = fit(Histogram, par, nbins = 200)
    bin_centers = collect(h.edges[1]) #.+ step(h.edges[1])/2
    p[i] = histogram(par, #getproperty(dsp_par, qc_pars[i]) 
            fill = true, color = :silver,
            xlabel = xl, 
            ylabel = "Counts", 
            bins = bin_centers, 
            label = "all",
            linecolor = :transparent)

    qc_frac_rejected = round(100 - prod(qc_cuts.qc_surv[qc_pars[i]][k] for k in keys(qc_cuts.qc_surv[qc_pars[i]]))*100, digits = 1)
    histogram!(10^-exp_i.*ustrip.(getproperty(dsp_par, qc_pars[i])[.!qc_cuts.wvf_keep[qc_pars[i]]]), 
                fill = true, color = :orange, 
                alpha = 0.8, linecolor = :orange,
                bins = bin_centers,
                label = "rejected $qc_frac_rejected %"
                ) 
    
end
ptot = plot(p..., layout = (2,2), size = (1000, 700), legend = :best, margin = 5mm)
savefig(ptot, pltpath * "dsp_pars_qc_cuts.png")

# plot some NOT cut waveforms 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
dsp_config = DSPConfig(dataprod_config(asic).dsp(filekeys[1]).default)
nplts = 25
wvf = read_ldata(asic, DataTier(:raw), category, period, run, channel)
waveforms_keep = wvf.waveform[findall(qc_cuts.wvf_keep.all)][1:nplts*10]
bl_stats = signalstats.(waveforms_keep, leftendpoint(dsp_config.bl_window), rightendpoint(dsp_config.bl_window))
waveforms_keep = shift_waveform.(waveforms_keep, -bl_stats.mean)
deconv_flt = InvCRFilter(mvalue(asic.par.rpars.pz[period, run, channel].τ))
wvfs_pz = deconv_flt.(waveforms_keep)
pltpath_kept = pltpath * "wvfs/kept/"
if !ispath(pltpath_kept)
    mkpath(pltpath_kept)
end
for j = 1:nplts
    pw = Vector{Plots.Plot}(undef, 10)
    for i= 1:10
        wvf_idx = 10*(j-1)+i
        @info wvf_idx
        pw[i] = plot(waveforms_keep[wvf_idx], label = false, legend = false)
        plot!(wvfs_pz[wvf_idx], label = false, color = :red2)
    end
    plot(pw..., layout = (5,2), plot_title = "Waveforms: QC passed ",
            size = (600, 900), 
            left_margin = 5mm, 
            margins = 2mm)
    savefig(pltpath_kept * "waveforms_notcut_$j.png")
end

# plot some cut waveforms 
nplts = 25
wvf = read_ldata(asic, DataTier(:raw), category, period, run, channel)
waveforms_cut = wvf.waveform[findall(.!qc_cuts.wvf_keep.all)][1:nplts*10]
bl_stats = signalstats.(waveforms_cut, leftendpoint(dsp_config.bl_window), rightendpoint(dsp_config.bl_window))
waveforms_cut = shift_waveform.(waveforms_cut, -bl_stats.mean)
deconv_flt = InvCRFilter(mvalue(asic.par.rpars.pz[period, run, channel].τ))
wvfs_pz = deconv_flt.(waveforms_cut)
pltpath_cut = pltpath * "wvfs/cut/"
if !ispath(pltpath_cut)
    mkpath(pltpath_cut)
end
for j = 1:nplts
    pw = Vector{Plots.Plot}(undef, 10)
    for i= 1:10
        wvf_idx = 10*(j-1)+i
        @info wvf_idx
        pw[i] = plot(waveforms_keep[wvf_idx], label = false, legend = false)
        plot!(wvfs_pz[wvf_idx], label = false, color = :red2)
    end
    plot(pw..., layout = (5,2), plot_title = "Waveforms QC cut ",
            size = (600, 900), 
            left_margin = 5mm, 
            margins = 2mm)
    savefig(pltpath_cut * "waveforms_cut_$j.png")
end


